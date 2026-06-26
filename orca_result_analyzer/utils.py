import os

# Known dummy / pseudo-atom labels used in ORCA and other QC output files.
_DUMMY_SYMBOLS: frozenset[str] = frozenset(
    {"*", "-", "X", "DA", "DU", "DUM", "DUMMY", "Q", "BQ", "LP"}
)


def get_default_export_path(base_path, suffix="_analyzed", extension=""):
    """
    Generates a default export path based on the input file path.
    Example: 'job.out' -> 'job_analyzed.csv'
    """
    if not base_path:
        return ""
    dirname = os.path.dirname(base_path)
    filename_base = os.path.splitext(os.path.basename(base_path))[0]
    new_filename = f"{filename_base}{suffix}{extension}"
    return os.path.join(dirname, new_filename)


def normalize_atom_symbol(raw: str) -> str:
    """Return a valid RDKit atom symbol for *raw*, mapping dummy labels to '*'.

    Handles:
    - Labels with a colon suffix  e.g. ``X:1`` → ``*``
    - Known dummy labels          e.g. ``DA``, ``DU``, ``BQ`` → ``*``
    - Unknown / non-periodic      e.g. ``Xx`` → ``*``
    - Normal elements             e.g. ``Fe``, ``C`` → returned as-is (capitalised)
    """
    sym = raw.strip()
    # Strip ORCA-style colon-suffixed labels like "X:1" or "C:2"
    if ":" in sym:
        sym = sym.split(":")[0]
    if sym.upper() in _DUMMY_SYMBOLS:
        return "*"
    sym = sym.capitalize()
    try:
        from rdkit import Chem

        if Chem.GetPeriodicTable().GetAtomicNumber(sym) <= 0:
            return "*"
    except Exception:  # noqa: BLE001
        return "*"
    return sym
