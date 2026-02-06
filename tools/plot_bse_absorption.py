import argparse
import numpy as np
import matplotlib.pyplot as plt


# Enable LaTeX for matplotlib
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    #"font.serif": ["Times New Roman"],
    "font.size": 20,
    "axes.labelsize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 20,
    "axes.titlesize": 20
})

#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

def read_absorption_eh(path: str):
    """
    absorption_eh.dat columns:
      1 omega (eV)
      2 eps2
      3 eps1
      4 JDOS
    """
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    omega = data[:, 0]
    eps2  = data[:, 1]
    eps1  = data[:, 2]
    jdos  = data[:, 3]
    return omega, eps2, eps1, jdos

def read_eigenvalues(path: str):
    """
    eigenvalues.dat columns:
      eig (eV), abs(dipole)^2, Re(dipole), Im(dipole)
    """
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    e_ex  = data[:, 0]
    dip2  = data[:, 1]
    dipre = data[:, 2]
    dipim = data[:, 3]
    return e_ex, dip2, dipre, dipim

def gaussian_broaden(x_grid, centers, weights, sigma):
    """
    Sum_i weights_i * exp(-(x-centers_i)^2/(2*sigma^2))
    (not normalized to area=1 unless you want it)
    """
    if sigma <= 0:
        raise ValueError("sigma must be > 0 for broadening.")
    x = x_grid[:, None]
    c = centers[None, :]
    w = weights[None, :]
    g = np.exp(-0.5 * ((x - c) / sigma) ** 2)
    y = np.sum(w * g, axis=1)
    return y


def parse_args():
    # absorption_file = "absorption_eh.dat"
    # eigen_file      = "eigenvalues.dat"

    p = argparse.ArgumentParser(
        description="Plot BerkeleyGW absorption spectrum with exciton dipole strengths."
    )
    p.add_argument("--abs", dest="absorption_file", default="absorption_eh.dat",
                   help="Path to absorption_eh.dat (default: absorption_eh.dat)")
    p.add_argument("--eig", dest="eigen_file", default="eigenvalues.dat",
                   help="Path to eigenvalues.dat (default: eigenvalues.dat)")
    p.add_argument("--emin", type=float, default=None,
                   help="Minimum energy (eV) for plotting window")
    p.add_argument("--emax", type=float, default=None,
                   help="Maximum energy (eV) for plotting window")
    p.add_argument("--sigma", type=float, default=0.02,
                   help="Gaussian broadening sigma in eV for exciton spectrum (default: 0.05)")
    p.add_argument("--no-broaden", action="store_true",
                   help="Disable Gaussian broadening curve (sticks only)")
    p.add_argument("--no-scale", action="store_true",
                   help="Do not scale exciton |dipole|^2 to match eps2 max (plot raw values)")
    p.add_argument("--pad", type=float, default=0.1,
                   help="Energy padding (eV) for selecting excitons around window (default: 0.1)")
    p.add_argument("--save", default="absoprtion.pdf", help="Output figure filename")
    return p.parse_args()

def main():
    args = parse_args()


    # --- read files ---
    omega, eps2, eps1, jdos = read_absorption_eh(args.absorption_file)
    assert len(omega) == len(eps2) == len(eps1) == len(jdos)

    #print("omega = ", omega)
    #print("eps2  = ", eps1)
    #print("eps1  = ", eps2)
    #print("jdos  = ", jdos)

    e_ex, dip2, dipre, dipim = read_eigenvalues(args.eigen_file)
    assert len(e_ex) == len(dip2) == len(dipre) == len(dipim)
    # print("e_ex = ", e_ex)

    # --- choose plotting window ---
    file_emin, file_emax = float(np.min(omega)), float(np.max(omega))
    emin = file_emin if args.emin is None else args.emin
    emax = file_emax if args.emax is None else args.emax
    # print("plot window is ", emin, emax)
    if emin >= emax:
        raise ValueError(f"Invalid window: emin ({emin}) must be < emax ({emax}).")

    # --- stick selection / scaling options ---
    # Keep only excitons in range (a bit padded)
    pad = 0.1
    mask = (e_ex >= emin - pad) & (e_ex <= emax + pad)
    mask2 = (omega >= emin - pad) & (omega <= emax + pad)
    e_ex_plot = e_ex[mask]
    dip2_plot = dip2[mask]
    eps2_plot = eps2[mask2]
    ymax = max(eps2_plot)


    # Optional: scale dipole^2 to comparable visual magnitude
    # (purely for plotting; doesn't change physics)
    if len(dip2_plot) > 0 and np.max(dip2_plot) > 0:
        scale = np.max(eps2_plot) / np.max(dip2_plot) if np.max(eps2_plot) > 0 else 1.0
    else:
        scale = 1.0
    # scale *= 0.25
    dip2_scaled = dip2_plot * scale

    # --- optional broadening for exciton sticks ---
    do_broaden = (not args.no_broaden)
    print("do_broaden = ", do_broaden)
    sigma = float(args.sigma)

    e_grid = omega  # use same grid as absorption for easy overlay
    if do_broaden and len(e_ex_plot) > 0:
        dip_broaden = gaussian_broaden(e_grid, e_ex_plot, dip2_scaled, sigma=sigma)
    else:
        dip_broaden = None

    # --- plot ---
    fig, ax1 = plt.subplots(figsize=(8, 5))

    # Absorption: eps2
    ax1.plot(omega, eps2, linewidth=2, label=r"Im$\varepsilon_2(\omega)$", color="blue")
    ax1.set_xlabel("Energy (eV)")
    ax1.set_ylabel(r"Absorption (Im$\{\varepsilon\}$)")

    # Second y-axis for dipole^2
    ax2 = ax1.twinx()
    ax2.set_ylabel(r"Exciton intensity ($|d_i|^2$)", color="red")

    # Sticks
    if len(e_ex_plot) > 0:
        ax2.vlines(e_ex_plot, 0.0, dip2_scaled, linewidth=1.0, alpha=0.7, label=r"$|d_i|^2$", color='red')

    ### Broadened curve (optional)
    #if dip_broaden is not None:
    #    ax2.plot(e_grid, dip_broaden, linewidth=2, linestyle="--", label=f"Broadened $|d_i|^2$ ($\sigma=${sigma:.2f} eV)", color="red")

    # Limits
    ax1.set_xlim(emin, emax)
    ax1.set_ylim(0, ymax*1.05)
    ax2.set_ylim(0, ymax*1.05)

    # Combined legend
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="best", frameon=False)

    ax1.set_title("Absorption spectrum and exciton dipole strengths")
    fig.tight_layout()
    # plt.show()
    #fig.savefig(args.save, dpi=300, bbox_inches="tight")
    fig.savefig(args.save, bbox_inches="tight")

if __name__ == "__main__":
    main()

