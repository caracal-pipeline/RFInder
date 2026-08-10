"""
Regression tests for Rfinder.setArgs() (rfinder/main.py).

Scoped deliberately to changes made during the mgpls-illuminati integration/testing effort
(see RFINDER_CHANGES.md) - not a general test suite for the pre-existing codebase.
"""

import copy
import os

import yaml

from rfinder.main import DEFAULT_CONFIG, RFINDER_DIR, Rfinder

# setArgs() unconditionally touches many cfg_par keys beyond input_dir/output_dir (label,
# cleanup_enable, chunking, plot options, ...) and expects a fully-populated cfg_par as its
# baseline - exactly like real usage, where it's only ever called after an already-loaded
# default/user config, never with a bare dict. Load the real installed default once and deep-copy
# it per test so tests can't leak state into each other.
_DEFAULT_CFG = yaml.load(open(os.path.join(RFINDER_DIR, DEFAULT_CONFIG)), Loader=yaml.Loader)


def make_rfinder():
    rfi = Rfinder()
    rfi.cfg_par = copy.deepcopy(_DEFAULT_CFG)
    return rfi


def test_input_dir_sets_workdir():
    """-idir/--input-dir must update cfg_par['general']['workdir'].

    Regression test for a real bug: setArgs() checked kwargs.get('indir') (a key that never
    exists - click/scabha populate 'input_dir') instead of kwargs.get('input_dir'), so the
    flag was silently a no-op. Fixed 2026-08-10.
    """
    rfi = make_rfinder()
    rfi.setArgs({"input_dir": "/some/input/path"})
    assert rfi.cfg_par["general"]["workdir"] == "/some/input/path"


def test_output_dir_sets_outdir():
    """-odir/--output-dir must update cfg_par['general']['outdir'].

    Regression test for the same class of bug as test_input_dir_sets_workdir: setArgs()
    checked kwargs.get('outdir') (never a real key) instead of kwargs.get('output_dir').
    Fixed 2026-08-10.
    """
    rfi = make_rfinder()
    rfi.setArgs({"output_dir": "/some/output/path"})
    assert rfi.cfg_par["general"]["outdir"] == "/some/output/path"


def test_input_output_dir_combined_and_unset_safe():
    """Both flags apply independently when set together, and leave workdir/outdir untouched
    (rather than clobbering them) when absent/falsy."""
    rfi = make_rfinder()
    rfi.setArgs({"input_dir": "/in", "output_dir": "/out"})
    assert rfi.cfg_par["general"]["workdir"] == "/in"
    assert rfi.cfg_par["general"]["outdir"] == "/out"

    workdir_before = rfi.cfg_par["general"]["workdir"]
    outdir_before = rfi.cfg_par["general"]["outdir"]
    rfi.setArgs({"input_dir": None, "output_dir": None})
    assert rfi.cfg_par["general"]["workdir"] == workdir_before
    assert rfi.cfg_par["general"]["outdir"] == outdir_before
