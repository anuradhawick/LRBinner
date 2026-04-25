"""Tests for split_contigs in runners_utils."""

import os
import tempfile
import textwrap
import pytest
from Bio import SeqIO

from lrbinner.mbcclr_utils.runners_utils import split_contigs


def _setup_output(tmp):
    """Create the expected output directory structure."""
    frags = os.path.join(tmp, "fragments")
    os.makedirs(frags, exist_ok=True)
    return tmp


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n{seq}\n")


class TestSplitContigs:
    def test_short_contig_not_split(self):
        """A contig < 5000 bp should appear as one fragment."""
        with tempfile.TemporaryDirectory() as tmp:
            _setup_output(tmp)
            fa = os.path.join(tmp, "contigs.fasta")
            _write_fasta(fa, {"ctg1": "ACGT" * 100})  # 400 bp
            groups, parent = split_contigs(fa, tmp)
            assert "ctg1" in groups
            assert len(groups["ctg1"]) == 1

    def test_long_contig_is_split(self):
        """A contig >= 5000 bp should be split into multiple fragments."""
        with tempfile.TemporaryDirectory() as tmp:
            _setup_output(tmp)
            fa = os.path.join(tmp, "contigs.fasta")
            _write_fasta(fa, {"ctg1": "ACGT" * 2000})  # 8000 bp
            groups, parent = split_contigs(fa, tmp)
            assert len(groups["ctg1"]) > 1

    def test_output_fasta_created(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_output(tmp)
            fa = os.path.join(tmp, "contigs.fasta")
            _write_fasta(fa, {"ctg1": "ACGT" * 100})
            split_contigs(fa, tmp)
            assert os.path.isfile(os.path.join(tmp, "fragments", "contigs.fasta"))

    def test_parent_mapping(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_output(tmp)
            fa = os.path.join(tmp, "contigs.fasta")
            _write_fasta(fa, {"ctg1": "ACGT" * 100, "ctg2": "TTTT" * 100})
            groups, parent = split_contigs(fa, tmp)
            for idx in groups["ctg1"]:
                assert parent[idx] == "ctg1"
            for idx in groups["ctg2"]:
                assert parent[idx] == "ctg2"

    def test_unique_fragment_indices(self):
        with tempfile.TemporaryDirectory() as tmp:
            _setup_output(tmp)
            fa = os.path.join(tmp, "contigs.fasta")
            _write_fasta(fa, {
                "ctg1": "ACGT" * 2000,
                "ctg2": "TTTT" * 2000,
            })
            groups, parent = split_contigs(fa, tmp)
            all_indices = groups["ctg1"] + groups["ctg2"]
            assert len(all_indices) == len(set(all_indices)), "Fragment indices must be unique"
