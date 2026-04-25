"""Tests for cluster_utils (normalize, calc_distances, smaller_indices)."""

import math
import pytest
import numpy as np
import torch

from lrbinner.mbcclr_utils.cluster_utils import normalize, calc_distances, smaller_indices


class TestNormalize:
    def test_output_type_is_tensor(self):
        mat = np.random.rand(10, 8).astype(np.float32)
        result = normalize(mat)
        assert isinstance(result, torch.Tensor)

    def test_accepts_tensor_input(self):
        mat = torch.rand(10, 8)
        result = normalize(mat)
        assert isinstance(result, torch.Tensor)

    def test_row_norms_after_normalize(self):
        """After normalization each row should have L2 norm == 1/sqrt(2)."""
        mat = torch.rand(20, 8)
        norm_mat = normalize(mat)
        row_norms = norm_mat.norm(dim=1)
        expected = 1.0 / math.sqrt(2)
        for n in row_norms.tolist():
            assert math.isclose(n, expected, abs_tol=1e-5)

    def test_zero_row_handled(self):
        """An all-zero row must not produce NaN after normalization."""
        mat = torch.zeros(5, 8)
        mat[1] = torch.rand(8)
        result = normalize(mat)
        assert not torch.isnan(result).any()

    def test_does_not_modify_original(self):
        mat = torch.rand(10, 8)
        original = mat.clone()
        normalize(mat)
        assert torch.equal(mat, original)


class TestCalcDistances:
    def test_self_distance_is_zero(self):
        mat = torch.rand(10, 8)
        mat = normalize(mat)
        dists = calc_distances(mat, 0)
        assert math.isclose(dists[0].item(), 0.0, abs_tol=1e-6)

    def test_distances_non_negative(self):
        mat = normalize(torch.rand(20, 8))
        dists = calc_distances(mat, 0)
        assert (dists >= 0).all()

    def test_distances_at_most_one(self):
        """Cosine distance is bounded by [0, 1] for normalised positive vectors."""
        mat = normalize(torch.rand(20, 8))
        dists = calc_distances(mat, 0)
        assert (dists <= 1.0 + 1e-6).all()

    def test_length_equals_matrix_rows(self):
        mat = normalize(torch.rand(15, 8))
        dists = calc_distances(mat, 0)
        assert len(dists) == 15


class TestSmallerIndices:
    def test_returns_tensor(self):
        dists = torch.tensor([0.1, 0.5, 0.2, 0.8])
        result = smaller_indices(dists, 0.3)
        assert isinstance(result, torch.Tensor)

    def test_correct_indices_returned(self):
        dists = torch.tensor([0.1, 0.5, 0.2, 0.8])
        result = smaller_indices(dists, 0.3)
        assert set(result.tolist()) == {0, 2}

    def test_empty_when_threshold_too_low(self):
        dists = torch.tensor([0.5, 0.6, 0.7])
        result = smaller_indices(dists, 0.1)
        assert len(result) == 0

    def test_all_when_threshold_high(self):
        dists = torch.tensor([0.1, 0.2, 0.3])
        result = smaller_indices(dists, 1.0)
        assert len(result) == 3
