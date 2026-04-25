"""Tests for k-mer vector computation via pykmertools."""

import os
import tempfile
import math
import pytest
import pykmertools as kt

from lrbinner.mbcclr_utils.runners_utils import run_kmers


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_fasta(path, records):
    """Write a dict {name: seq} to a FASTA file."""
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n{seq}\n")


# Expected canonical-k-mer vector lengths from pykmertools
EXPECTED_LEN = {3: 32, 4: 136, 5: 512}


# ---------------------------------------------------------------------------
# OligoComputer unit tests (pykmertools directly)
# ---------------------------------------------------------------------------

class TestOligoComputer:
    def test_vector_length_k3(self):
        vec = kt.OligoComputer(3).vectorise_one("ACGTACGT")
        assert len(vec) == EXPECTED_LEN[3]

    def test_vector_length_k4(self):
        vec = kt.OligoComputer(4).vectorise_one("ACGTACGTACGT")
        assert len(vec) == EXPECTED_LEN[4]

    def test_vector_length_k5(self):
        vec = kt.OligoComputer(5).vectorise_one("ACGTACGTACGTACGT")
        assert len(vec) == EXPECTED_LEN[5]

    def test_normalised_sum(self):
        """Frequency vector should sum to 1.0 for a non-trivial sequence."""
        vec = kt.OligoComputer(3).vectorise_one("ACGTACGT")
        assert math.isclose(sum(vec), 1.0, abs_tol=1e-6)

    def test_homopolymer_nonzero(self):
        """A pure-A homopolymer should produce a nonzero vector."""
        vec = kt.OligoComputer(3).vectorise_one("AAAAAAAAAA")
        assert sum(vec) > 0


# ---------------------------------------------------------------------------
# run_kmers integration tests
# ---------------------------------------------------------------------------

class TestRunKmers:
    def _run(self, records, k):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = os.path.join(tmp, "reads.fasta")
            _write_fasta(fasta, records)
            run_kmers(fasta, tmp, k, 1)
            out_path = os.path.join(tmp, "profiles", "com_profs")
            assert os.path.isfile(out_path), "com_profs file not created"
            with open(out_path) as fh:
                lines = [l.strip() for l in fh if l.strip()]
            return lines

    def test_line_count_matches_reads(self):
        records = {"r1": "ACGTACGT", "r2": "TTTTAAAA", "r3": "GCGCGCGC"}
        lines = self._run(records, k=3)
        assert len(lines) == 3

    def test_vector_length_k3(self):
        records = {"r1": "ACGTACGTACGT"}
        lines = self._run(records, k=3)
        assert len(lines[0].split()) == EXPECTED_LEN[3]

    def test_vector_length_k4(self):
        records = {"r1": "ACGTACGTACGT"}
        lines = self._run(records, k=4)
        assert len(lines[0].split()) == EXPECTED_LEN[4]

    def test_vector_length_k5(self):
        records = {"r1": "ACGTACGTACGTACGT"}
        lines = self._run(records, k=5)
        assert len(lines[0].split()) == EXPECTED_LEN[5]

    def test_normalised(self):
        records = {"r1": "ACGTACGTACGT"}
        lines = self._run(records, k=3)
        values = [float(v) for v in lines[0].split()]
        assert math.isclose(sum(values), 1.0, abs_tol=1e-5)

    def test_homopolymer_run(self):
        """Homopolymer sequences should produce a valid (non-error) output line."""
        records = {"r1": "AAAAAAAAAA"}
        lines = self._run(records, k=3)
        assert len(lines) == 1
        values = [float(v) for v in lines[0].split()]
        assert len(values) == EXPECTED_LEN[3]
