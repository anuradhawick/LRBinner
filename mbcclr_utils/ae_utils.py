import numpy as np
import torch
from torch import nn, optim
from torch.utils.data import DataLoader
from torch.utils.data.dataset import TensorDataset
from sklearn.preprocessing import MinMaxScaler
import logging
from tqdm import tqdm
from collections import defaultdict
import os
import json

logger = logging.getLogger('LRBinner')

h_params = open(os.path.dirname(__file__) + '/hyper_params.json')
h_params = json.load(h_params)


def make_data_loader(covs, profs, *, batch_size=1024, drop_last=True, shuffle=True, cuda=False):
    # Scaling profs
    profs = MinMaxScaler().fit_transform(profs)
    covs = MinMaxScaler().fit_transform(covs)

    covs = torch.from_numpy(covs).float()
    profs = torch.from_numpy(profs).float()
    idx = torch.arange(len(covs))

    n_workers = 4 if cuda else 1

    dataset = TensorDataset(covs, profs, idx)
    return DataLoader(dataset=dataset, batch_size=batch_size, drop_last=drop_last,
                      shuffle=shuffle, pin_memory=cuda, num_workers=n_workers)


class VAE(nn.Module):
    def __init__(self, cov_size, prof_size, *, latent_dims=8, hidden_layers=[128, 128], constraints=None, device='cpu'):
        super(VAE, self).__init__()

        self.cov_size = cov_size
        self.prof_size = prof_size
        self.hidden_layers = hidden_layers
        self.latent_dims = latent_dims
        self.dropout = 0.1
        # TODO
        # self.beta = 5

        self.encoderlayers = nn.ModuleList()
        self.encodernorms = nn.ModuleList()
        self.decoderlayers = nn.ModuleList()
        self.decodernorms = nn.ModuleList()

        # Encoding layers
        for nin, nout in zip([self.cov_size + self.prof_size] + self.hidden_layers, self.hidden_layers):
            self.encoderlayers.append(nn.Linear(nin, nout))
            self.encodernorms.append(nn.BatchNorm1d(nout))

        # Latent layers
        self.mu = nn.Linear(self.hidden_layers[-1], self.latent_dims)
        self.logsigma = nn.Linear(self.hidden_layers[-1], self.latent_dims)

        # Decoding layers
        for nin, nout in zip([self.latent_dims] + self.hidden_layers[::-1], self.hidden_layers[::-1]):
            self.decoderlayers.append(nn.Linear(nin, nout))
            self.decodernorms.append(nn.BatchNorm1d(nout))

        # Final layer
        self.outputlayer = nn.Linear(
            self.hidden_layers[0], self.cov_size + self.prof_size)

        # Activation functions
        self.relu = nn.LeakyReLU()
        self.softplus = nn.Softplus()
        self.softmax = nn.Softmax(dim=1)
        self.dropoutlayer = nn.Dropout(p=self.dropout)

        self.constraints = constraints
        self.device = device

        if self.constraints:
            self.constrained_idx = set()
            ml = np.array(constraints['ml'])
            self.ml_index = {}
            self.mnl_index = {}

            for [k, v] in ml:
                self.constrained_idx.add(k)
                self.constrained_idx.add(v)
                self.ml_index[(k, v)] = 1

            mnl = np.array(constraints['mnl'])

            for [k, v] in mnl:
                self.constrained_idx.add(k)
                self.constrained_idx.add(v)
                self.mnl_index[(k, v)] = 1

        self.to(self.device)


    def _search_index(self, batch_indices):
        must_link_pairs = []
        must_not_link_pairs = []
        batch_indices = batch_indices.tolist()
        valid_idx = list(sorted(self.constrained_idx.intersection(batch_indices)))
        idx_local_map = {i:n for n, i in enumerate(batch_indices)}

        must_link_pairs = []
        must_not_link_pairs = []

        for i in valid_idx:
            for j in valid_idx:
                if j>=i:
                    break

                if (i, j) in self.ml_index or (j, i) in self.ml_index:
                    must_link_pairs.append([idx_local_map[i], idx_local_map[j]])

                if (i, j) in self.mnl_index or (j, i) in self.mnl_index:
                    must_not_link_pairs.append([idx_local_map[i], idx_local_map[j]])
        
        must_link_pairs = np.array(must_link_pairs)
        must_not_link_pairs = np.array(must_not_link_pairs)

        return must_link_pairs, must_not_link_pairs
    

    def _encode(self, tensor):
        tensors = list()

        for encoderlayer, encodernorm in zip(self.encoderlayers, self.encodernorms):
            tensor = encodernorm(self.dropoutlayer(
                self.relu(encoderlayer(tensor))))
            tensors.append(tensor)

        # Latent layers
        mu = self.mu(tensor)
        logsigma = self.softplus(self.logsigma(tensor))

        return mu, logsigma

    def encode(self, data_loader):
        self.eval()

        covs, profs, indices = data_loader.dataset.tensors
        length = len(covs)

        # from vamb
        latent = np.empty((length, self.latent_dims), dtype=np.float32)

        row = 0
        with torch.no_grad():
            for covs, profs, indices in data_loader:
                covs = covs.to(self.device)
                profs = profs.to(self.device)

                mu, logsigma = self.forward_predict(covs, profs)
                latent[row: row + len(mu)] = mu.to("cpu")
                row += len(mu)

        assert row == length
        return latent

    def reparameterize(self, mu, logsigma):
        epsilon = torch.randn(mu.size(0), mu.size(1), device=self.device)
        epsilon.requires_grad = True
        latent = mu + epsilon * torch.exp(logsigma/2)

        return latent

    def _decode(self, tensor):
        tensors = list()

        for decoderlayer, decodernorm in zip(self.decoderlayers, self.decodernorms):
            tensor = decodernorm(self.dropoutlayer(
                self.relu(decoderlayer(tensor))))
            tensors.append(tensor)

        reconstruction = self.outputlayer(tensor)
        covs_out = reconstruction.narrow(1, 0, self.cov_size)
        profs_out = reconstruction.narrow(1, self.cov_size, self.prof_size)

        return covs_out, profs_out

    def forward(self, covs, profs):
        tensor = torch.cat((covs, profs), 1)
        mu, logsigma = self._encode(tensor)
        latent = self.reparameterize(mu, logsigma)

        covs_out, profs_out = self._decode(latent)

        return covs_out, profs_out, mu, logsigma

    def forward_predict(self, covs, profs):
        tensor = torch.cat((covs, profs), 1)
        mu, logsigma = self._encode(tensor)

        return mu, logsigma

    def trainepoch(self, data_loader, epoch, optimizer, batchsteps, logfile):
        self.train()

        epoch_loss = 0
        epoch_kldloss = 0
        epoch_mseloss = 0
        epoch_celoss = 0

        # Doubling the batch size at each batchstep
        if epoch in batchsteps:
            data_loader = DataLoader(dataset=data_loader.dataset,
                                     batch_size=data_loader.batch_size * 2,
                                     shuffle=True,
                                     drop_last=True,
                                     num_workers=data_loader.num_workers,
                                     pin_memory=data_loader.pin_memory
                                     )

        for n, (covs_in, profs_in, indices) in enumerate(data_loader):
            covs_in = covs_in.to(self.device)
            profs_in = profs_in.to(self.device)

            covs_in.requires_grad = True
            profs_in.requires_grad = True

            optimizer.zero_grad()

            covs_out, profs_out, mu, logsigma = self(covs_in, profs_in)

            loss, e_cov, e_comp, kld = self.calc_loss(
                covs_in, covs_out, profs_in, profs_out, mu, logsigma, indices)

            loss.backward()
            optimizer.step()

            epoch_loss += loss.data.item()
            epoch_kldloss += kld.data.item()
            epoch_mseloss += e_comp.data.item()
            epoch_celoss += e_cov.data.item()

        logger.debug(f'Epoch: {epoch + 1:4} Loss: {epoch_loss / (1+len(data_loader)):.6f}\tEC: {epoch_celoss / (1+len(data_loader)):.7f}\tEP: {epoch_mseloss / (1+len(data_loader)):.6f}\tKLD: {epoch_kldloss / (1+len(data_loader)):.4f}\tBatchsize: {data_loader.batch_size}')

        return data_loader

    def calc_loss(self, covs_in, covs_out, profs_in, profs_out, mu, logsigma, indices):
        loss_ml = 0
        loss_mnl = 0

        if self.constraints is not None:
            must_link_pairs, must_not_link_pairs = self._search_index(indices)
            
            if len(must_link_pairs) > 0:
                loss_ml = (mu[must_link_pairs.T[0]] - mu[must_link_pairs.T[1]]).pow(2).sum(dim=1).mean()
            if len(must_link_pairs) > 0:
                loss_mnl = torch.max(torch.tensor(0).to(self.device), 10 - (mu[must_not_link_pairs.T[0]] - mu[must_not_link_pairs.T[1]]).pow(2).sum(dim=1).mean())

        e_cov = (covs_out - covs_in).pow(2).sum(dim=1).mean()
        # TODO
        # e_cov_weight = 0.1
        e_cov_weight = h_params[str(self.prof_size)]["e_cov_weight"]

        e_comp = (profs_out - profs_in).pow(2).sum(dim=1).mean()
        # e_comp_weight = 1
        e_comp_weight = h_params[str(self.prof_size)]["e_comp_weight"]

        kld = -0.5 * (1 + logsigma - mu.pow(2) -
                      logsigma.exp()).sum(dim=1).mean()
        # kld_weight = 1 / (self.prof_size * self.beta)
        kld_weight = h_params[str(self.prof_size)]["kld_weight"]

        loss = e_cov * e_cov_weight + e_comp * e_comp_weight + kld * kld_weight + loss_ml + loss_mnl

        return loss, e_cov, e_comp, kld

    def trainmodel(self, dataloader, *, nepochs=500, lrate=1e-3,
                   batchsteps=[25, 75, 150, 300], logfile=None, save_path=None):
        batchsteps_set = set(batchsteps)
        optimizer = optim.Adam(self.parameters(), lr=lrate)

        for epoch in tqdm(range(nepochs), total=nepochs, desc="Training VAE"):
            dataloader = self.trainepoch(
                dataloader, epoch, optimizer, batchsteps_set, logfile)
        self.save(save_path)

    def save(self, save_path):
        state = {'cov_size': self.cov_size,
                 'prof_size': self.prof_size,
                 'dropout': self.dropout,
                 'hidden_layers': self.hidden_layers,
                 'latent_dims': self.latent_dims,
                 'state': self.state_dict(),
                 }

        torch.save(state, save_path)


def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


def vae_encode(output, latent_dims, hidden_layers, epochs, constraints, cuda):
    comp_profiles = np.load(f"{output}/profiles/com_profs.npy")
    cov_profiles = np.load(f"{output}/profiles/cov_profs.npy")

    device = "cpu"

    if cuda:
        device = "cuda"

    vae = VAE(cov_profiles.shape[1], comp_profiles.shape[1],
              latent_dims=latent_dims, 
              hidden_layers=hidden_layers, 
              constraints=constraints, 
              device=device)

    logger.debug(f"Model param count = {count_parameters(vae)}") 
    logger.debug(vae)

    dloader = make_data_loader(cov_profiles, comp_profiles, cuda=cuda)
    vae.trainmodel(
        dloader, save_path=f"{output}/model.pt", nepochs=epochs, batchsteps=[50, 100, 150])

    dloader_encode = make_data_loader(
        cov_profiles, comp_profiles, drop_last=False, shuffle=False, cuda=cuda)
    latent = vae.encode(dloader_encode)

    np.save(f"{output}/latent", latent)
