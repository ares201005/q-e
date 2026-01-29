#!/usr/bin/env python3
# plot_phonon_bands.py

import re
import argparse
import numpy as np
import matplotlib.pyplot as plt


# Enable LaTeX for matplotlib
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 20,
    "axes.labelsize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 16,
    "axes.titlesize": 20
})


def pretty_label(s: str) -> str:
    mapping = {
        'G':'G','Gamma':'G','GAMMA':'G',
        'M':'M','K':'K','A':'A','L':'L','H':'H','X':'X','Y':'Y','Z':'Z'
    }
    return mapping.get(s, s)

def parse_matdyn_freq(path):
    """
    Supports:
      (A) '&plot nbnd=.., nks=.. /' format
          Each q-block has one line: 'qx qy qz s' (s = cumulative distance),
          followed by *one or more* lines totaling `nbnd` frequencies.
      (B) Legacy format: q-line then many 'freq(i) = ... [unit]' lines.
    Returns: x (nq,), q (nq,3), w (nq,nmodes), unit (str)
    """
    # keep blank-line filtering, but preserve ordering
    with open(path, 'r') as f:
        lines = [ln.rstrip() for ln in f if ln.strip()]

    float_re = re.compile(r'([+-]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)')  # also accept Fortran D

    # ---------- Format (A): &plot ----------
    if lines and lines[0].lstrip().startswith('&plot'):
        hdr = lines[0]
        m = re.search(r'nbnd\s*=\s*(\d+)\s*,\s*nks\s*=\s*(\d+)', hdr)
        if not m:
            raise ValueError("Could not read nbnd/nks from &plot header.")
        nbnd, nks = int(m.group(1)), int(m.group(2))

        q = np.zeros((nks, 3), float)
        x = np.zeros(nks, float)
        w = np.zeros((nks, nbnd), float)

        i = 1  # first line after header
        iq = 0
        while iq < nks and i < len(lines):
            # --- q-line: expect at least 4 floats (qx qy qz s)
            qvals = [float(v.replace('D', 'E').replace('d', 'e')) for v in float_re.findall(lines[i])]
            if len(qvals) < 4:
                raise ValueError(f"Expected 4 floats on q-line at block {iq+1}, got: '{lines[i]}'")
            q[iq, :] = qvals[:3]
            x[iq]    = qvals[3]
            i += 1

            # --- freq lines: may span multiple lines until we collect nbnd numbers
            freqs = []
            while len(freqs) < nbnd and i < len(lines):
                vals = [float(v.replace('D', 'E').replace('d', 'e')) for v in float_re.findall(lines[i])]
                # It’s possible (though unlikely) a malformed file jumps to next q before finishing freqs.
                # We trust header counts and just keep appending floats here.
                if vals:
                    freqs.extend(vals)
                i += 1

            if len(freqs) < nbnd:
                raise ValueError(f"Block {iq+1}: expected {nbnd} frequencies, got {len(freqs)}.")
            w[iq, :] = freqs[:nbnd]
            iq += 1

        if iq != nks:
            raise ValueError(f"Read {iq} q-points but header says nks={nks}.")
        unit = "cm$^{-1}$"  # not printed in this format
        return x, q, w, unit

    # ---------- Format (B): legacy ----------
    q_list, w_blocks, cur_w = [], [], []
    unit = "cm$^{-1}$"
    three_float = re.compile(
        r'^\s*([+-]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)\s+'
        r'([+-]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)\s+'
        r'([+-]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)\s*$'
    )
    frq_line = re.compile(r'^\s*(?:freq|omega)\s*\(\s*\d+\s*\)\s*=\s*([+-]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)(?:\s*\[(.*?)\])?', re.I)

    for line in lines:
        m_q = three_float.match(line)
        if m_q:
            if cur_w:
                w_blocks.append(cur_w); cur_w = []
            q_list.append([float(m_q.group(1).replace('D','E')),
                           float(m_q.group(2).replace('D','E')),
                           float(m_q.group(3).replace('D','E'))])
            continue
        m_f = frq_line.match(line)
        if m_f:
            val = float(m_f.group(1).replace('D','E'))
            cur_w.append(val)
            if m_f.group(2):
                u = m_f.group(2).strip()
                unit = u.replace('cm-1', r'cm$^{-1}$').replace('THz','THz')
    if cur_w: w_blocks.append(cur_w)
    if not q_list or not w_blocks:
        raise ValueError("Could not parse q-points or frequencies (unrecognized format).")
    nmodes = max(len(b) for b in w_blocks)
    if any(len(b)!=nmodes for b in w_blocks):
        raise ValueError("Inconsistent mode counts among q-points in legacy format.")
    q = np.array(q_list, float)
    w = np.array(w_blocks, float)
    dq = np.linalg.norm(q[1:] - q[:-1], axis=1)
    x = np.concatenate([[0.0], np.cumsum(dq)])
    return x, q, w, unit

def make_ticks(x, segments, labels):
    if len(labels) != len(segments)  + 1:
        raise ValueError("labels length must be segments length + 1 (one per vertex).")
    ends = np.cumsum(segments) #- 1
    tick_positions = [x[0]] + [x[i] for i in ends]
    tick_labels = [pretty_label(s) for s in labels]
    return tick_positions, tick_labels, ends

def plot_bands(x, w, tick_positions, tick_labels, ends, title=None, unit='cm$^{-1}$', outfile=None, show=True):
    plt.figure(figsize=(6.2, 4.6))
    for imode in range(w.shape[1]):
        plt.plot(x, w[:, imode], lw=3.0)
    for idx in ends[:-1]:
        plt.axvline(x=x[idx], linestyle='-', linewidth=0.6)
    plt.xticks(tick_positions, tick_labels)
    plt.xlim(x[0], x[-1])
    plt.ylim(bottom=0.0)
    plt.ylabel(f'Frequency ({unit})')
    # plt.xlabel('Wave vector')
    if title: plt.title(title)
    plt.grid(alpha=0.2, linestyle='--', linewidth=0.5)
    plt.tight_layout()
    # plt.show()
    # if show: plt.show()
    if outfile: plt.savefig(outfile, dpi=300, bbox_inches='tight')

def main():
    ap = argparse.ArgumentParser(description="Plot phonon band structure from QE matdyn flfrq.")
    ap.add_argument('--freq', required=True, help="Path to matdyn flfrq (bands) file")
    ap.add_argument('--segments', nargs='+', type=int, required=True, help="Points per path segment, e.g. 60 60 60")
    ap.add_argument('--labels',   nargs='+', required=True, help="Vertex labels; length = segments+1 (use G for Gamma)")
    ap.add_argument('--unit', default=None, help="Override y-axis unit (e.g., THz)")
    ap.add_argument('--title', default=None)
    ap.add_argument('--outfile', default=None)
    ap.add_argument('--no-show', action='store_true')
    args = ap.parse_args()

    x, q, w, unit = parse_matdyn_freq(args.freq)
    nq = w.shape[0]
    if nq != sum(args.segments) + 1:
        raise ValueError(f"nq={nq} points, but sum(segments)={sum(args.segments)}.")
    tick_positions, tick_labels, ends = make_ticks(x, args.segments, args.labels)
    plot_bands(x, w, tick_positions, tick_labels, ends,
               title=args.title, unit=(args.unit or unit),
               outfile=args.outfile, show=not args.no_show)

if __name__ == '__main__':
    main()
