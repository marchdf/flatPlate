#!/usr/bin/env python3

# ========================================================================
#
# Imports
#
# ========================================================================
import argparse
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ========================================================================
#
# Some defaults variables
#
# ========================================================================
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Times')
cmap_med = ['#F15A60', '#7AC36A', '#5A9BD4', '#FAA75B',
            '#9E67AB', '#CE7058', '#D77FB4', '#737373']
cmap = ['#EE2E2F', '#008C48', '#185AA9', '#F47D23',
        '#662C91', '#A21D21', '#B43894', '#010202']
dashseq = [(None, None), [10, 5], [10, 4, 3, 4], [
    3, 3], [10, 4, 3, 4, 3, 4], [3, 3], [3, 3]]
markertype = ['s', 'd', 'o', 'p', 'h']


# ========================================================================
#
# Function definitions
#
# ========================================================================


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == '__main__':

    # ========================================================================
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='A simple plot tool')
    parser.add_argument(
        '-s', '--show', help='Show the plots', action='store_true')
    args = parser.parse_args()

    # ========================================================================
    # Nalu output
    fname = 'coeffs.dat'
    df = pd.read_csv(fname)
    df['h'] = np.sqrt(1. / df['N2'])

    plt.figure(0)
    plt.semilogx(df['h'], df['cd'], lw=2, color=cmap[0])
    plt.semilogx(df['h'], df['cd'], markertype[0],
                 mec=cmap[0], mfc=cmap[0], ms=10)

    plt.figure(1)
    plt.semilogx(np.sqrt(df['h']), df['cf'], lw=2, color=cmap[0])
    plt.semilogx(np.sqrt(df['h']), df['cf'], markertype[0],
                 mec=cmap[0], mfc=cmap[0], ms=10)

    # wall cf
    df = pd.read_csv('545x385/results/wall_coeffs.dat')
    plt.figure(2)
    plt.plot(df['x'], df['cf'], lw=2, color=cmap[0])

    # ========================================================================
    # NASA CFL3D output
    fname = os.path.join(os.path.abspath('nasa_data'), 'coeffs_cfl3d.dat')
    df = pd.read_csv(fname)

    plt.figure(0)
    plt.semilogx(df['h'], df['cd'], lw=2, color=cmap[1])
    plt.semilogx(df['h'], df['cd'], markertype[1],
                 mec=cmap[1], mfc=cmap[1], ms=10)

    plt.figure(1)
    plt.semilogx(np.sqrt(df['h']), df['cf'], lw=2, color=cmap[1])
    plt.semilogx(np.sqrt(df['h']), df['cf'], markertype[1],
                 mec=cmap[1], mfc=cmap[1], ms=10)

    # wall cf
    fname = os.path.join(os.path.abspath('nasa_data'), 'wall_cf_cfl3d.dat')
    df = pd.read_csv(fname)
    plt.figure(2)
    plt.plot(df['x'], df['cf'], lw=2, color=cmap[1])

    # ========================================================================
    # NASA FUN3D output
    fname = os.path.join(os.path.abspath('nasa_data'), 'coeffs_fun3d.dat')
    df = pd.read_csv(fname)

    plt.figure(0)
    plt.semilogx(df['h'], df['cd'], lw=2, color=cmap[2])
    plt.semilogx(df['h'], df['cd'], markertype[2],
                 mec=cmap[2], mfc=cmap[2], ms=10)

    plt.figure(1)
    plt.semilogx(np.sqrt(df['h']), df['cf'], lw=2, color=cmap[2])
    plt.semilogx(np.sqrt(df['h']), df['cf'], markertype[2],
                 mec=cmap[2], mfc=cmap[2], ms=10)

    # wall cf
    fname = os.path.join(os.path.abspath('nasa_data'), 'wall_cf_fun3d.dat')
    df = pd.read_csv(fname)
    plt.figure(2)
    plt.plot(df['x'], df['cf'], lw=2, color=cmap[2])

    # ========================================================================
    # Format the plots
    plt.figure(0)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$C_d$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    plt.tight_layout()
    plt.savefig('cd.pdf', format='pdf')

    plt.figure(1)
    ax = plt.gca()
    plt.xlabel(r"$h~[-]$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$C_f$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    plt.tight_layout()
    plt.savefig('cf.pdf', format='pdf')

    plt.figure(2)
    ax = plt.gca()
    plt.xlabel(r"$x~[m]$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$C_f$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    plt.ylim([0.002, 0.006])
    plt.tight_layout()
    plt.savefig('wall_cf.pdf', format='pdf')

    if args.show:
        plt.show()
