"""Command line: read an output file without a window and print what was found.

``moleditpy -r orca_result_analyzer parse FILE`` prints a summary, and
``--json`` the parsed data that every dialog reads, for scripts and ``jq``;
``xyz FILE`` writes the final structure, and ``xyz --all`` every step.
No QApplication is created.
"""

import argparse
import json
import os
import sys
from collections import Counter

from . import PLUGIN_VERSION
from .loading import read_orca_text
from .parser import IMAGINARY_FREQ_THRESHOLD, OrcaParser

PROGRAM = "ORCA"
_PACKAGE = __package__ or "orca_result_analyzer"
_HARTREE_TO_EV = 27.211386245988


def _prog():
    """The command that started us, so usage text names one that exists."""
    if "moleditpy.plugins.plugin_runner" in sys.modules:
        return f"moleditpy -r {_PACKAGE}"
    return f"python -m {_PACKAGE}"


def _empty(value):
    return value is None or value == [] or value == {} or value == ""


#: Keys the parser fills with a default whatever the file holds.
_ALWAYS_SET = ("termination_status", "charge", "mult")


def _found_anything(data):
    """Whether parsing recognised any part of the file; flags do not count."""
    return any(
        not _empty(value) and not isinstance(value, bool)
        for key, value in data.items()
        if key not in _ALWAYS_SET
    )


def parse_file(path):
    """Read and parse *path*, returning the parser's data dict."""
    parser = OrcaParser()
    parser.load_from_memory(read_orca_text(path), path)
    return parser.data


def _formula(atoms):
    counts = Counter(atoms)
    order = [s for s in ("C", "H") if s in counts] if "C" in counts else []
    order += sorted(s for s in counts if s not in order)
    return "".join(s + (str(counts[s]) if counts[s] > 1 else "") for s in order)


def _optimisation_converged(data):
    """ORCA marks the stationary point it reached with an opt_final step."""
    return any(step.get("type") == "opt_final" for step in data["scan_steps"])


def _frontier(orbitals):
    """HOMO and LUMO energies in Eh, alpha or restricted, or None."""
    homo = lumo = None
    for orb in orbitals or []:
        if orb.get("spin", "restricted") not in ("restricted", "alpha"):
            continue
        energy = orb.get("energy_eh", orb.get("energy"))
        if energy is None:
            continue
        if orb.get("occupation", orb.get("occ", 0)) > 0:
            homo = energy if homo is None else max(homo, energy)
        else:
            lumo = energy if lumo is None else min(lumo, energy)
    return homo, lumo


def summary_lines(path, data):
    """Human-readable lines describing what *data* holds."""
    lines = [f"File:        {path}"]
    version = data.get("version") or "(version not found)"
    if not version.startswith(PROGRAM):
        version = f"{PROGRAM} {version}"
    lines.append(f"Program:     {version}")
    lines.append(f"Status:      {data.get('termination_status') or 'unknown'}")
    atoms = data.get("atoms") or []
    if atoms:
        lines.append(
            f"Molecule:    {_formula(atoms)}, {len(atoms)} atoms, "
            f"charge {data.get('charge', 0)}, multiplicity {data.get('mult', 1)}"
        )
    if data.get("scf_energy") is not None:
        lines.append(f"SCF energy:  {data['scf_energy']:.10f} Eh")
    s2 = data.get("spin_s2")
    if isinstance(s2, dict) and s2.get("actual") is not None:
        lines.append(f"<S**2>:      {s2['actual']:.4f} (ideal {s2.get('ideal')})")
    homo, lumo = _frontier(data.get("orbital_energies"))
    if homo is not None and lumo is not None:
        lines.append(
            f"HOMO/LUMO:   {homo * _HARTREE_TO_EV:.3f} / {lumo * _HARTREE_TO_EV:.3f} eV "
            f"(gap {(lumo - homo) * _HARTREE_TO_EV:.3f} eV)"
        )
    dipole = data.get("dipole")
    if isinstance(dipole, dict) and dipole.get("magnitude_debye") is not None:
        lines.append(f"Dipole:      {dipole['magnitude_debye']:.4f} D")
    steps = data.get("scan_steps") or []
    if steps:
        kinds = Counter(step.get("type", "step") for step in steps)
        lines.append(
            "Steps:       "
            + ", ".join(f"{n} {kind}" for kind, n in sorted(kinds.items()))
        )
        if kinds.get("opt_cycle") and not data.get("is_scan"):
            state = "converged" if _optimisation_converged(data) else "not converged"
            lines.append(f"Optimised:   {state}")
    modes = [m for m in data.get("frequencies") or [] if m.get("freq")]
    if modes:
        imaginary = [m["freq"] for m in modes if m["freq"] < -IMAGINARY_FREQ_THRESHOLD]
        text = f"Frequencies: {len(modes)} modes, {len(imaginary)} imaginary"
        if imaginary:
            text += f" (lowest {min(imaginary):.1f} cm-1)"
        lines.append(text)
    thermal = data.get("thermal") or {}
    if thermal.get("gibbs") is not None:
        lines.append(
            f"Gibbs:       {thermal['gibbs']:.8f} Eh at {thermal.get('temperature', '?')} K"
        )
    excited = data.get("tddft") or data.get("excitation_energies") or []
    if excited:
        lines.append(f"Excited:     {len(excited)} states")
    if data.get("nmr_shielding"):
        lines.append(f"NMR:         {len(data['nmr_shielding'])} shieldings")
    charges = data.get("charges") or {}
    if charges:
        lines.append("Charges:     " + ", ".join(charges))
    if data.get("nbo_orbitals"):
        lines.append(
            f"NBO:         {len(data['nbo_orbitals'])} orbitals, "
            f"{len(data.get('nbo_perturbation') or [])} E(2) entries"
        )
    present = sorted(key for key, value in data.items() if not _empty(value))
    lines.append("")
    lines.append("Sections:    " + ", ".join(present))
    return lines


def _build_parser():
    parser = argparse.ArgumentParser(
        prog=_prog(),
        description=f"{PROGRAM} Result Analyzer {PLUGIN_VERSION}: analyse {PROGRAM} output "
        "without opening a window.",
    )
    commands = parser.add_subparsers(dest="command", metavar="COMMAND")
    parse = commands.add_parser(
        "parse",
        help="read an output file and print a summary or its data",
        description=f"Read an output file from {PROGRAM} and print a summary, or with "
        "--json the parsed data the analysis windows use.",
    )
    parse.add_argument("file", help=f"{PROGRAM} output file")
    parse.add_argument(
        "--json",
        action="store_true",
        help="print the parsed data as JSON instead of the summary",
    )
    parse.add_argument(
        "-s",
        "--section",
        action="append",
        metavar="NAME",
        help="with --json, only this top-level key (repeatable); "
        "the summary lists the keys present",
    )
    parse.add_argument(
        "-o",
        "--output",
        metavar="PATH",
        help="write to PATH instead of standard output",
    )
    parse.add_argument(
        "--indent", type=int, default=1, metavar="N", help="JSON indent (default 1)"
    )
    xyz = commands.add_parser(
        "xyz",
        help="write the final structure, or every step, as XYZ",
        description="Write the final structure as XYZ in Angstrom, or with --all "
        "every step as a multi-frame XYZ: the points of a scan or IRC, otherwise "
        "the optimisation cycles.",
    )
    xyz.add_argument("file", help=f"{PROGRAM} output file")
    xyz.add_argument(
        "--all",
        action="store_true",
        help="every scan or IRC point, or every optimisation cycle, one frame "
        "each with its energy in the comment line",
    )
    xyz.add_argument(
        "-o",
        "--output",
        metavar="PATH",
        help="write to PATH instead of standard output",
    )
    return parser


def _load(path):
    """Parse *path*, or report why not; returns (data, None) or (None, status)."""
    try:
        data = parse_file(path)
    except OSError as exc:
        print(f"Error: cannot read {path}: {exc}", file=sys.stderr)
        return None, 1
    if not _found_anything(data):
        print(f"Error: {path} does not look like {PROGRAM} output.", file=sys.stderr)
        return None, 1
    return data, None


def _write(args, text):
    if args.output:
        with open(args.output, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text + "\n")
    else:
        sys.stdout.write(text + "\n")


def xyz_frame(atoms, coords, comment):
    """One XYZ frame; the comment is kept to a single line."""
    lines = [str(len(atoms)), " ".join(comment.splitlines())]
    lines += [
        f"{symbol:<2} {x:15.8f} {y:15.8f} {z:15.8f}"
        for symbol, (x, y, z) in zip(atoms, coords)
    ]
    return "\n".join(lines)


def _step_comment(step):
    if step.get("type") == "scan_step" and step.get("scan_step_id") is not None:
        parts = [f"point {step['scan_step_id']}"]
    else:
        parts = [f"{step.get('type', 'step')} {step.get('step', '')}".strip()]
    if step.get("energy") is not None:
        parts.append(f"E = {step['energy']:.10f} Eh")
    if step.get("scan_coord") is not None:
        parts.append(f"coordinate = {step['scan_coord']}")
    return "  ".join(parts)


def trajectory_steps(data):
    """The steps to write: a scan's or IRC's points, otherwise every cycle."""
    steps = [
        s for s in data.get("scan_steps") or [] if s.get("atoms") and s.get("coords")
    ]
    points = [s for s in steps if s.get("type") == "scan_step"]
    return points or steps


def _xyz_command(args):
    data, status = _load(args.file)
    if data is None:
        return status
    name = os.path.basename(args.file)
    if args.all:
        steps = trajectory_steps(data)
        if not steps:
            print(
                f"Error: {name} has no optimisation, scan or IRC steps.",
                file=sys.stderr,
            )
            return 1
        frames = [
            xyz_frame(s["atoms"], s["coords"], f"{name}  {_step_comment(s)}")
            for s in steps
        ]
    else:
        if not data.get("atoms") or not data.get("coords"):
            print(f"Error: {name} has no structure.", file=sys.stderr)
            return 1
        comment = name
        if data.get("scf_energy") is not None:
            comment += f"  E = {data['scf_energy']:.10f} Eh"
        frames = [xyz_frame(data["atoms"], data["coords"], comment)]
    _write(args, "\n".join(frames))
    return 0


def _parse_command(args):
    if args.section and not args.json:
        print("Error: --section needs --json.", file=sys.stderr)
        return 2
    data, status = _load(args.file)
    if data is None:
        return status
    if args.json:
        if args.section:
            missing = [name for name in args.section if name not in data]
            if missing:
                print(
                    f"Error: no section {', '.join(missing)}; "
                    f"available: {', '.join(sorted(data))}",
                    file=sys.stderr,
                )
                return 2
            data = {name: data[name] for name in args.section}
        text = json.dumps(data, indent=args.indent, ensure_ascii=False, default=str)
    else:
        text = "\n".join(summary_lines(os.path.abspath(args.file), data))
    _write(args, text)
    return 0


def main(argv=None):
    """Run the command line with *argv* (default sys.argv[1:]); return the status."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.command == "parse":
        return _parse_command(args)
    if args.command == "xyz":
        return _xyz_command(args)
    parser.print_help()
    return 2
