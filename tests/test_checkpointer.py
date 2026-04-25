"""Tests for the Checkpointer class in runners_utils."""

import os
import pickle
import tempfile
import pytest

from lrbinner.mbcclr_utils.runners_utils import Checkpointer


class TestCheckpointerNew:
    def test_empty_on_creation(self):
        with tempfile.TemporaryDirectory() as tmp:
            cp = Checkpointer(os.path.join(tmp, "cp.pkl"))
            assert cp.completed == {}

    def test_should_run_step_unknown_stage(self):
        with tempfile.TemporaryDirectory() as tmp:
            cp = Checkpointer(os.path.join(tmp, "cp.pkl"))
            assert cp.should_run_step("1_1", {"k": 3}) is True

    def test_log_marks_stage_complete(self):
        with tempfile.TemporaryDirectory() as tmp:
            cp = Checkpointer(os.path.join(tmp, "cp.pkl"))
            cp.log("1_1", {"k": 3})
            assert cp.should_run_step("1_1", {"k": 3}) is False

    def test_param_change_triggers_rerun(self):
        with tempfile.TemporaryDirectory() as tmp:
            cp = Checkpointer(os.path.join(tmp, "cp.pkl"))
            cp.log("1_1", {"k": 3})
            assert cp.should_run_step("1_1", {"k": 4}) is True

    def test_log_clears_later_stages(self):
        with tempfile.TemporaryDirectory() as tmp:
            cp = Checkpointer(os.path.join(tmp, "cp.pkl"))
            cp.log("1_1", {"k": 3})
            cp.log("2_1", {"x": 1})
            # re-logging stage 1 should clear stage 2
            cp.log("1_1", {"k": 3})
            assert "2_1" not in cp.completed

    def test_persistence(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "cp.pkl")
            cp = Checkpointer(path)
            cp.log("1_1", {"k": 3})
            # Load a new instance from same file
            cp2 = Checkpointer(path, _load_to_resume=True)
            assert cp2.should_run_step("1_1", {"k": 3}) is False

    def test_no_resume_ignores_existing_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "cp.pkl")
            cp = Checkpointer(path)
            cp.log("1_1", {"k": 3})
            # Without _load_to_resume the state is fresh
            cp2 = Checkpointer(path, _load_to_resume=False)
            assert cp2.completed == {}

    def test_str_representation(self):
        with tempfile.TemporaryDirectory() as tmp:
            cp = Checkpointer(os.path.join(tmp, "cp.pkl"))
            cp.log("1_1", {"k": 3})
            assert "1_1" in str(cp)
