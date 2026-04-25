"""Tests for the lrbinner package metadata and top-level imports."""

import lrbinner


def test_version():
    assert lrbinner.__version__ == "2.1.0"


def test_version_is_string():
    assert isinstance(lrbinner.__version__, str)


def test_version_format():
    parts = lrbinner.__version__.split(".")
    assert len(parts) == 3
    assert all(p.isdigit() for p in parts)


def test_submodule_imports():
    from lrbinner.mbcclr_utils import runners_utils
    from lrbinner.mbcclr_utils import ae_utils
    from lrbinner.mbcclr_utils import cluster_utils
    from lrbinner.metacoag_utils import marker_gene_utils
