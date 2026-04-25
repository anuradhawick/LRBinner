"""Tests for ae_utils (make_data_loader, VAE)."""

import pytest
import numpy as np
import torch

from lrbinner.mbcclr_utils.ae_utils import make_data_loader, VAE


N = 200       # number of samples — must be > batch_size (1024 default) so use small batch
COV_DIM = 4
PROF_DIM = 32
LATENT = 8


def _random_arrays(n=N, cov_dim=COV_DIM, prof_dim=PROF_DIM):
    rng = np.random.default_rng(42)
    covs = rng.random((n, cov_dim)).astype(np.float32)
    profs = rng.random((n, prof_dim)).astype(np.float32)
    return covs, profs


class TestMakeDataLoader:
    def test_returns_dataloader(self):
        covs, profs = _random_arrays()
        dl = make_data_loader(covs, profs, batch_size=32, drop_last=False, shuffle=False)
        assert dl is not None

    def test_batch_shapes(self):
        covs, profs = _random_arrays()
        dl = make_data_loader(covs, profs, batch_size=32, drop_last=False, shuffle=False)
        cov_batch, prof_batch, idx_batch = next(iter(dl))
        assert cov_batch.shape[1] == COV_DIM
        assert prof_batch.shape[1] == PROF_DIM

    def test_values_in_unit_interval(self):
        """After MinMaxScaler the tensor values should lie in [0, 1]."""
        covs, profs = _random_arrays()
        dl = make_data_loader(covs, profs, batch_size=N, drop_last=False, shuffle=False)
        cov_batch, prof_batch, _ = next(iter(dl))
        assert cov_batch.min().item() >= 0.0
        assert cov_batch.max().item() <= 1.0 + 1e-6
        assert prof_batch.min().item() >= 0.0
        assert prof_batch.max().item() <= 1.0 + 1e-6

    def test_index_tensor_length(self):
        covs, profs = _random_arrays()
        dl = make_data_loader(covs, profs, batch_size=N, drop_last=False, shuffle=False)
        _, _, idx = next(iter(dl))
        assert len(idx) == N


class TestVAE:
    def _make_vae(self):
        return VAE(cov_size=COV_DIM, prof_size=PROF_DIM, latent_dims=LATENT, hidden_layers=[64, 64])

    def test_instantiation(self):
        vae = self._make_vae()
        assert vae is not None

    def test_forward_output_shapes(self):
        vae = self._make_vae()
        vae.eval()
        bs = 16
        covs = torch.randn(bs, COV_DIM)
        profs = torch.randn(bs, PROF_DIM)
        with torch.no_grad():
            cov_out, prof_out, mu, logsigma = vae(covs, profs)
        assert cov_out.shape == (bs, COV_DIM)
        assert prof_out.shape == (bs, PROF_DIM)
        assert mu.shape == (bs, LATENT)
        assert logsigma.shape == (bs, LATENT)

    def test_latent_shape_via_encode(self):
        covs, profs = _random_arrays()
        dl = make_data_loader(covs, profs, batch_size=32, drop_last=False, shuffle=False)
        vae = self._make_vae()
        vae.eval()
        latent = vae.encode(dl)
        assert latent.shape == (N, LATENT)

    def test_gpu_skip_if_unavailable(self):
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
        vae = VAE(cov_size=COV_DIM, prof_size=PROF_DIM, latent_dims=LATENT,
                  hidden_layers=[64, 64], device="cuda")
        covs = torch.randn(16, COV_DIM).cuda()
        profs = torch.randn(16, PROF_DIM).cuda()
        cov_out, prof_out, mu, logsigma = vae(covs, profs)
        assert mu.shape == (16, LATENT)
