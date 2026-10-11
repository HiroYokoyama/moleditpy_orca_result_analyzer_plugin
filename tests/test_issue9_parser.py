"""Regression coverage for issue #9 using real ORCA 5/6 outputs."""
from pathlib import Path
import pytest
from tests._parser_loader import load_standalone_parser

OrcaParser = load_standalone_parser().OrcaParser
SAMPLES = Path(__file__).parent / "sample_outputs"


def output(name):
    return (SAMPLES / name).read_text(encoding="utf-8")


@pytest.mark.parametrize("name", ["benzene-opt-ene.out", "benzene-opt-ene_5.out", "benzene-opt-eneQZ.out"])
def test_every_mo_has_all_basis_coefficients(name):
    parser = OrcaParser()
    parser.load_from_memory(output(name))
    size = sum(2 * shell["l"] + 1 for shell in parser.data["basis_set_shells"])
    assert size in (114, 522)
    assert len(parser.data["mo_coeffs"]) == size
    assert all(len(mo["coeffs"]) == size for mo in parser.data["mo_coeffs"].values())
