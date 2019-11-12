#!/usr/bin/env python3

# ========================================================================
#
# Imports
#
# ========================================================================
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ========================================================================
#
# Some defaults variables
#
# ========================================================================
plt.rc("text", usetex=True)
cmap_med = [
    "#F15A60",
    "#7AC36A",
    "#5A9BD4",
    "#FAA75B",
    "#9E67AB",
    "#CE7058",
    "#D77FB4",
    "#737373",
]
cmap = [
    "#EE2E2F",
    "#008C48",
    "#185AA9",
    "#F47D23",
    "#662C91",
    "#A21D21",
    "#B43894",
    "#010202",
]
dashseq = [
    (None, None),
    [10, 5],
    [10, 4, 3, 4],
    [3, 3],
    [10, 4, 3, 4, 3, 4],
    [3, 3],
    [3, 3],
]
markertype = ["s", "d", "o", "p", "h"]


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == "__main__":

    # ========================================================================
    # Parse arguments
    parser = argparse.ArgumentParser(description="A simple plot tool")
    parser.add_argument("-s", "--show", help="Show the plots", action="store_true")
    args = parser.parse_args()

    # ========================================================================
    # NASA CFL3D output
    fname = os.path.join(os.path.abspath("nasa_data"), "coeffs_cfl3d.dat")
    df = pd.read_csv(fname)

    plt.figure(0)
    plt.semilogx(
        df["h"],
        df["cd"],
        lw=2,
        color=cmap[0],
        marker=markertype[0],
        mec=cmap[0],
        mfc=cmap[0],
        ms=10,
        label="NASA CFL3D",
    )

    plt.figure(1)
    plt.semilogx(
        df["h"],
        df["cf"],
        lw=2,
        color=cmap[0],
        marker=markertype[0],
        mec=cmap[0],
        mfc=cmap[0],
        ms=10,
        label="NASA CFL3D",
    )

    # wall cf
    fname = os.path.join(os.path.abspath("nasa_data"), "wall_cf_cfl3d.dat")
    df = pd.read_csv(fname)
    plt.figure(2)
    p = plt.plot(df["x"], df["cf"], lw=2, color=cmap[0], label="NASA CFL3D")
    p[0].set_dashes(dashseq[0])

    # profiles
    fname = os.path.join(os.path.abspath("nasa_data"), "mut_0.97_cfl3d.dat")
    df = pd.read_csv(fname)
    plt.figure(3)
    p = plt.plot(df["mut"], df["y"], lw=2, color=cmap[0], label="NASA CFL3D")
    p[0].set_dashes(dashseq[0])

    fname = os.path.join(os.path.abspath("nasa_data"), "f1f2_0.97_cfl3d.dat")
    df = pd.read_csv(fname)
    plt.figure(4)
    p = plt.plot(df["F1"], df["y"], lw=2, color=cmap[0], label="NASA CFL3D")
    p[0].set_dashes(dashseq[0])

    fname = os.path.join(os.path.abspath("nasa_data"), "sdr_tke_0.97_cfl3d.dat")
    df = pd.read_csv(fname)
    plt.figure(5)
    p = plt.plot(df["tke"], df["y"], lw=2, color=cmap[0], label="NASA CFL3D")
    p[0].set_dashes(dashseq[0])
    plt.figure(6)
    p = plt.plot(df["sdr"], df["y"], lw=2, color=cmap[0], label="NASA CFL3D")
    p[0].set_dashes(dashseq[0])

    fname = os.path.join(os.path.abspath("nasa_data"), "yp_up_0.97_cfl3d.dat")
    df = pd.read_csv(fname)
    plt.figure(7)
    p = plt.plot(df["logyp"], df["up"], lw=2, color=cmap[0], label="NASA CFL3D")
    p[0].set_dashes(dashseq[0])

    fname = os.path.join(os.path.abspath("nasa_data"), "u_0.97_cfl3d.dat")
    df = pd.read_csv(fname)
    plt.figure(8)
    p = plt.plot(df["u"], df["y"], lw=2, color=cmap[0], label="NASA CFL3D")
    p[0].set_dashes(dashseq[0])

    # ========================================================================
    # NASA FUN3D output
    fname = os.path.join(os.path.abspath("nasa_data"), "coeffs_fun3d.dat")
    df = pd.read_csv(fname)

    plt.figure(0)
    plt.semilogx(
        df["h"],
        df["cd"],
        lw=2,
        color=cmap[1],
        marker=markertype[1],
        mec=cmap[1],
        mfc=cmap[1],
        ms=10,
        label="NASA FUN3D",
    )

    plt.figure(1)
    plt.semilogx(
        df["h"],
        df["cf"],
        lw=2,
        color=cmap[1],
        marker=markertype[1],
        mec=cmap[1],
        mfc=cmap[1],
        ms=10,
        label="NASA FUN3D",
    )

    # wall cf
    fname = os.path.join(os.path.abspath("nasa_data"), "wall_cf_fun3d.dat")
    df = pd.read_csv(fname)
    plt.figure(2)
    p = plt.plot(df["x"], df["cf"], lw=2, color=cmap[1], label="NASA FUN3D")
    p[0].set_dashes(dashseq[1])

    # profiles
    fname = os.path.join(os.path.abspath("nasa_data"), "mut_0.97_fun3d.dat")
    df = pd.read_csv(fname)
    plt.figure(3)
    p = plt.plot(df["mut"], df["y"], lw=2, color=cmap[1], label="NASA FUN3D")
    p[0].set_dashes(dashseq[1])

    fname = os.path.join(os.path.abspath("nasa_data"), "sdr_tke_0.97_fun3d.dat")
    df = pd.read_csv(fname)
    plt.figure(5)
    p = plt.plot(df["tke"], df["y"], lw=2, color=cmap[1], label="NASA FUN3D")
    p[0].set_dashes(dashseq[1])
    plt.figure(6)
    p = plt.plot(df["sdr"], df["y"], lw=2, color=cmap[1], label="NASA FUN3D")
    p[0].set_dashes(dashseq[1])

    # ========================================================================
    # Nalu output
    fname = "coeffs.dat"
    df = pd.read_csv(fname)
    print(df)

    plt.figure(0)
    plt.semilogx(
        df["h"],
        df["cd"],
        ls="-",
        lw=2,
        color=cmap[2],
        marker=markertype[2],
        mec=cmap[2],
        mfc=cmap[2],
        ms=10,
        label="Nalu",
    )

    plt.figure(1)
    plt.semilogx(
        df["h"],
        df["cf"],
        ls="-",
        lw=2,
        color=cmap[2],
        marker=markertype[2],
        mec=cmap[2],
        mfc=cmap[2],
        ms=10,
        label="Nalu",
    )

    # wall cf
    fdir = "545x385/results"
    df = pd.read_csv(os.path.join(fdir, "wall_coeffs.dat"))
    plt.figure(2)
    p = plt.plot(df["x"], df["cf"], lw=2, color=cmap[2], label="Nalu")
    p[0].set_dashes(dashseq[2])

    # profiles
    profiles = pd.read_csv(os.path.join(fdir, "profiles.dat"))
    plt.figure(3)
    p = plt.plot(profiles["mut"], profiles["z"], lw=2, color=cmap[2], label="Nalu")
    p[0].set_dashes(dashseq[2])

    plt.figure(4)
    p = plt.plot(profiles["F1"], profiles["z"], lw=2, color=cmap[2], label="Nalu")
    p[0].set_dashes(dashseq[2])

    plt.figure(5)
    p = plt.plot(profiles["tke"], profiles["z"], lw=2, color=cmap[2], label="Nalu")
    p[0].set_dashes(dashseq[2])

    plt.figure(6)
    p = plt.plot(profiles["sdr"], profiles["z"], lw=2, color=cmap[2], label="Nalu")
    p[0].set_dashes(dashseq[2])

    plt.figure(7)
    p = plt.plot(
        np.log10(profiles["yp"]), profiles["uxp"], lw=2, color=cmap[2], label="Nalu"
    )
    p[0].set_dashes(dashseq[2])

    plt.figure(8)
    p = plt.plot(profiles["ux"], profiles["z"], lw=2, color=cmap[2], label="Nalu")
    p[0].set_dashes(dashseq[2])

    # ========================================================================
    # Format the plots
    plt.figure(0)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$C_d$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("cd.pdf", format="pdf", dpi=300)
    plt.savefig("cd.png", format="png", dpi=300)

    plt.figure(1)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$C_f$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("cf.pdf", format="pdf", dpi=300)
    plt.savefig("cf.png", format="png", dpi=300)

    plt.figure(2)
    ax = plt.gca()
    plt.xlabel(r"$x~[m]$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$C_f$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    plt.ylim([0.002, 0.006])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("wall_cf.pdf", format="pdf", dpi=300)
    plt.savefig("wall_cf.png", format="png", dpi=300)

    plt.figure(3)
    ax = plt.gca()
    plt.xlabel(r"$\mu_t / \mu$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y~[m]$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    plt.ylim([0.0, 0.025])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("mut.pdf", format="pdf", dpi=300)
    plt.savefig("mut.png", format="png", dpi=300)

    plt.figure(4)
    ax = plt.gca()
    plt.xlabel(r"$F_1$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y~[m]$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    plt.ylim([0.0, 0.025])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("f1.pdf", format="pdf", dpi=300)
    plt.savefig("f1.png", format="png", dpi=300)

    plt.figure(5)
    ax = plt.gca()
    plt.xlabel(r"$k / a^2$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y~[m]$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    plt.xlim([1e-10, 1e-3])
    plt.ylim([1e-7, 1e-1])
    plt.xscale("log")
    plt.yscale("log")
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("tke.pdf", format="pdf", dpi=300)
    plt.savefig("tke.png", format="png", dpi=300)

    plt.figure(6)
    ax = plt.gca()
    plt.xlabel(r"$\omega \mu / (\rho a^2)$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y~[m]$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    plt.xlim([1e-8, 1e0])
    plt.ylim([1e-7, 1e-1])
    plt.xscale("log")
    plt.yscale("log")
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("sdr.pdf", format="pdf", dpi=300)
    plt.savefig("sdr.png", format="png", dpi=300)

    plt.figure(7)
    ax = plt.gca()
    plt.xlabel(r"$\log{(y^+)}$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$u^+$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    plt.xlim([-1, 5])
    plt.ylim([0, 30])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("up.pdf", format="pdf", dpi=300)
    plt.savefig("up.png", format="png", dpi=300)

    plt.figure(8)
    ax = plt.gca()
    plt.xlabel(r"$u / u_0$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y~[m]$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    plt.xlim([0, 1])
    plt.ylim([0, 0.04])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig("u.pdf", format="pdf", dpi=300)
    plt.savefig("u.png", format="png", dpi=300)

    if args.show:
        plt.show()
