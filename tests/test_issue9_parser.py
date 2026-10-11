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


def test_latest_equal_sized_basis_parameters():
    parser = OrcaParser()
    parser.data.update(atoms=["H"], coords=[[0, 0, 0]])
    def basis(exponent):
        return f"BASIS SET IN INPUT FORMAT\n--------------------\nNewGTO H\nS 1\n1 {exponent} 1.0\nend;\n"
    parser.lines = (basis(1.0) + basis(2.0)).splitlines()
    parser.parse_basis_set()
    assert parser.data["basis_set_shells"][0]["exps"] == [2.0]


def test_final_basis_and_coefficients_have_matching_size():
    parser = OrcaParser()
    parser.load_from_memory(output("benzene-opt-ene.out") + "\n" + output("benzene-opt-eneQZ.out"))
    size = sum(2 * shell["l"] + 1 for shell in parser.data["basis_set_shells"])
    assert size == 522
    assert all(len(mo["coeffs"]) == size for mo in parser.data["mo_coeffs"].values())
