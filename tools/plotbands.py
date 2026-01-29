#!/usr/bin/env python3
import argparse
import os
import re
from typing import List, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 20,
    "axes.labelsize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 20,
    "axes.titlesize": 20
})


HEADER_RE = re.compile(r'&\s*plot\b.*?nbnd\s*=\s*(\d+)\s*,\s*nks\s*=\s*(\d+)', re.IGNORECASE)

# Define your high-symmetry points (in crystal/BZ coords)
HIGH_SYM_POINTS = {
    "G": np.array([0.0, 0.0, 0.0]),
    "X": np.array([0.5, 0.0, 0.0]),
    "S": np.array([0.5, 0.5, 0.0]),
    "Y": np.array([0.0, 0.5, 0.0]),
}

def kpoint_label(kpoint, sym_points=HIGH_SYM_POINTS, tol=1e-4):
    """
    kpoint: iterable of length 3 (kx, ky, kz)
    sym_points: dict{label: np.array([kx,ky,kz])}
    tol: numerical tolerance for matching

    Returns: label string (e.g. 'G', 'X', ...) or None if no match.
    """
    kvec = np.array(kpoint, dtype=float)
    for label, ref in sym_points.items():
        if np.linalg.norm(kvec - ref) < tol:
            return label
    return None


def parse_bands_dat_qe(path: str) -> Tuple[np.ndarray, np.ndarray, int, int]:
    """
    Parse QE bands.dat in the format:
      &plot nbnd=NN, nks=MM /
                kx ky kz
        E1 E2 ... (over multiple lines, total nbnd)
      (repeat MM times)
    Returns:
      kfrac: (MM, 3) fractional k-points
      energies: (NN, MM)
      nbnd, nks
    """
    with open(path, 'r') as f:
        lines = [ln.rstrip() for ln in f]

    # Find header with nbnd and nks
    nbnd = nks = None
    i0 = 0
    for i, ln in enumerate(lines):
        m = HEADER_RE.search(ln)
        if m:
            nbnd = int(m.group(1))
            nks = int(m.group(2))
            i0 = i + 1
            break
    if nbnd is None or nks is None:
        raise ValueError("Could not find header '&plot nbnd=..., nks=...' in file.")

    kfrac = np.zeros((nks, 3), dtype=float)
    energies = np.zeros((nbnd, nks), dtype=float)

    i = i0
    ik = 0
    while ik < nks and i < len(lines):
        # Skip empty lines
        while i < len(lines) and not lines[i].strip():
            i += 1
        if i >= len(lines):
            break

        # Expect a line with 3 floats = k-point
        parts = lines[i].split()
        if len(parts) != 3:
            raise ValueError(f"Expected 3 floats for k-point at line {i+1}, got: '{lines[i]}'")
        try:
            kx, ky, kz = map(float, parts)
        except ValueError as e:
            raise ValueError(f"Invalid k-point line at {i+1}: {lines[i]}") from e
        kfrac[ik, :] = [kx, ky, kz]
        i += 1

        # Now accumulate nbnd energies possibly over multiple lines
        acc: List[float] = []
        while len(acc) < nbnd and i < len(lines):
            ln = lines[i].strip()
            i += 1
            if not ln:
                continue
            # Split numbers in this line
            tokens = ln.split()
            # Some files may have trailing comments; keep numeric tokens only
            vals = []
            for t in tokens:
                try:
                    vals.append(float(t))
                except ValueError:
                    # stop at first non-numeric token
                    break
            if not vals:
                continue
            acc.extend(vals)

        if len(acc) < nbnd:
            raise ValueError(f"Reached EOF before reading {nbnd} energies for k-point {ik+1}.")
        energies[:, ik] = acc[:nbnd]
        ik += 1

    if ik != nks:
        raise ValueError(f"Parsed {ik} k-points, but header says nks = {nks}.")

    return kfrac, energies, nbnd, nks


def recip_from_cell(cell: np.ndarray) -> np.ndarray:
    """
    Given cell (3x3) with rows a1,a2,a3 in Angstrom, return reciprocal lattice
    matrix B with rows b1,b2,b3 in 1/Angstrom (2*pi convention).
    """
    a1, a2, a3 = cell
    vol = np.dot(a1, np.cross(a2, a3))
    if abs(vol) < 1e-12:
        raise ValueError("Cell volume is nearly zero; check CELL_PARAMETERS.")
    b1 = 2*np.pi * np.cross(a2, a3) / vol
    b2 = 2*np.pi * np.cross(a3, a1) / vol
    b3 = 2*np.pi * np.cross(a1, a2) / vol
    return np.vstack([b1, b2, b3])


def parse_cell_arg(cell_str: str) -> np.ndarray:
    """
    Parse a cell string like: "a1x,a1y,a1z; a2x,a2y,a2z; a3x,a3y,a3z"
    into a (3,3) matrix (rows = a1,a2,a3), in Angstrom.
    """
    try:
        rows = [r.strip() for r in cell_str.split(';')]
        if len(rows) != 3:
            raise ValueError
        mat = []
        for r in rows:
            comps = [float(x) for x in r.split(',')]
            if len(comps) != 3:
                raise ValueError
            mat.append(comps)
        return np.array(mat, dtype=float)
    except Exception:
        raise ValueError("Failed to parse --cell. Use 'a1x,a1y,a1z; a2x,a2y,a2z; a3x,a3y,a3z' (Å).")


def cumulative_kdist(kfrac: np.ndarray,
                     bmat: Optional[np.ndarray] = None,
                     cell: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Compute cumulative distance along k-path.
    If bmat (3x3) is provided (rows are reciprocal vectors in 1/Å),
    convert fractional k to Cartesian: k_cart = kfrac @ bmat  (since rows are b-vectors),
    otherwise use fractional coordinates directly (unitless).
    Returns s of shape (nks,).
    """
    if bmat is not None:
        k_cart = kfrac @ bmat  # (nks,3)
    else:
        k_cart = kfrac

    if cell is not None:
        tmp = cell / cell[0,0]
        ktmp = kfrac @ tmp

    diffs = np.diff(k_cart, axis=0)
    seg = np.linalg.norm(diffs, axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])

    # get labels automatically
    nk = len(kfrac)
    labels = []
    for i in range(nk):
        k_label = kpoint_label(ktmp[i])
        if k_label is not None: labels.append([s[i], k_label])

    return s, labels


def auto_breaks(kdist: np.ndarray) -> List[int]:
    """
    Locate likely segment boundaries by large jumps in step size.
    Returns indices where vertical separators should be drawn.
    """
    if len(kdist) < 3:
        return []
    steps = np.diff(kdist)
    med = np.median(steps[steps > 1e-12]) if np.any(steps > 1e-12) else 0.0
    idx = []
    for i, st in enumerate(steps, start=1):
        if med > 0 and st > 5 * med:
            idx.append(i)
    return idx


def main():
    ap = argparse.ArgumentParser(description="Plot QE band structure from bands.dat (&plot nbnd=..., nks=...) format.")
    ap.add_argument("--input", default="bands.dat", help="Path to bands.dat")
    ap.add_argument("--save", default="bandstructure.pdf", help="Output figure filename")
    ap.add_argument("--dpi", type=int, default=300, help="DPI for raster outputs")
    ap.add_argument("--fermi", type=float, default=0.0, help="Fermi level in eV (subtract from energies)")
    ap.add_argument("--emin", type=float, default=None, help="ymin (eV)")
    ap.add_argument("--emax", type=float, default=None, help="ymax (eV)")

    cell_group = ap.add_mutually_exclusive_group()
    cell_group.add_argument("--alat", type=str, default=None,
                            help="Orthorhombic cell lengths 'a,b,c' in Å (CELL_PARAMETERS diagonal).")
    cell_group.add_argument("--cell", type=str, default=None,
                            help="Full 3x3 cell: 'a1x,a1y,a1z; a2x,a2y,a2z; a3x,a3y,a3z' in Å.")

    ap.add_argument("--labels", type=str, default=None,
                    help='Comma-separated high-symmetry labels, e.g. "Γ,X,S,Y,Γ,Z"')
    ap.add_argument("--nodes", type=str, default=None,
                    help="Comma-separated node positions along k (same units as x-axis).")
    ap.add_argument("--node-indices", type=str, default=None,
                    help="Comma-separated integer indices of k-points where ticks should be placed.")

    ap.add_argument("--figsize", type=float, nargs=2, metavar=("WIDTH", "HEIGHT"), default=(6, 4),
                    help="Figure size in inches, e.g. --figsize 8 6")

    ap.add_argument("--linewidth", type=float, default=2.0, help="Line width")
    args = ap.parse_args()

    prefix = args.input.split(".")[0]
    args.save = f"{prefix}_" + args.save

    # Parse bands.dat
    kfrac, energies, nbnd, nks = parse_bands_dat_qe(args.input)

    # Reciprocal lattice (optional)
    bmat = None
    if args.alat:
        try:
            a, b, c = [float(x) for x in args.alat.split(",")]
        except Exception:
            raise SystemExit("Failed to parse --alat. Use 'a,b,c' in Å.")
        cell = np.array([[a, 0.0, 0.0],
                         [0.0, b, 0.0],
                         [0.0, 0.0, c]], dtype=float)
        bmat = recip_from_cell(cell)
    elif args.cell:
        cell = parse_cell_arg(args.cell)
        bmat = recip_from_cell(cell)
    # print("bmat = ", bmat)

    # Compute cumulative k-distance
    kdist, tick_poslabels = cumulative_kdist(kfrac, bmat=bmat, cell=cell)
    x_label = r"$k$ (1/Å)" if bmat is not None else "k (fractional units)"

    # Fermi shift
    if args.fermi != 0.0:
        energies = energies - args.fermi

    # Plot
    fig, ax = plt.subplots()
    for ib in range(nbnd):
        ax.plot(kdist, energies[ib, :], linewidth=args.linewidth)

    # ax.set_xlabel(x_label)
    ax.set_ylabel("Energy (eV)")
    ax.grid(True, which="both", linestyle="-", alpha=0.3)

    # Vertical separators and xticks
    tick_positions: Optional[List[float]] = None
    tick_labels: Optional[List[str]] = None

    if args.nodes:
        try:
            tick_positions = [float(x) for x in args.nodes.split(",")]
        except Exception:
            raise SystemExit("Failed to parse --nodes. Provide comma-separated floats.")
    elif args.node_indices:
        try:
            idx = [int(x) for x in args.node_indices.split(",")]
            tick_positions = [kdist[i] for i in idx if 0 <= i < len(kdist)]
        except Exception:
            raise SystemExit("Failed to parse --node-indices. Provide comma-separated integers.")
    elif len(tick_poslabels) > 0:
        tick_positions = [label[0] for label in tick_poslabels]
    else:
        # Auto-detect large jumps
        idx = auto_breaks(kdist)
        tick_positions = [kdist[0]] + [kdist[i] for i in idx] + [kdist[-1]]

    if args.labels:
        tick_labels = [s.strip() for s in args.labels.split(",")]
        if len(tick_labels) != len(tick_positions):
            print("Warning: number of labels != number of tick positions; labels will be ignored.")
            tick_labels = None
    elif len(tick_poslabels) > 0:
        tick_labels = [label[1] for label in tick_poslabels]

    print("tick_positions =", tick_positions)
    print("tick_labels    =", tick_labels)

    # Draw separators and xticks
    xmin, xmax = float(np.min(kdist)), float(np.max(kdist))
    if tick_positions:
        for x in sorted(set(tick_positions)):
            if x > xmin + 1e-12 and x < xmax - 1e-12:
                ax.axvline(x=x, linestyle="--", linewidth=0.8, alpha=0.5)
        if tick_labels:
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels)

    # x-limit
    ax.set_xlim(xmin, xmax)
    # y-limits
    if args.emin is not None or args.emax is not None:
        ymin = args.emin if args.emin is not None else float(np.nanmin(energies)) - 0.5
        ymax = args.emax if args.emax is not None else float(np.nanmax(energies)) + 0.5
        ax.set_ylim(ymin, ymax)

    # Zero line if Fermi-shifted
    if args.fermi != 0.0:
        ax.axhline(0.0, linestyle="--", linewidth=0.8, alpha=0.6)

    fig.tight_layout()
    ext = os.path.splitext(args.save)[1].lower()
    if ext in [".pdf", ".png", ".jpg", ".jpeg", ".tif", ".tiff"]:
        fig.savefig(args.save, dpi=args.dpi, bbox_inches="tight")
    else:
        fig.savefig(args.save, bbox_inches="tight")
    print(f"Saved plot to: {args.save}")


if __name__ == "__main__":
    main()
