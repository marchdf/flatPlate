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
    # NASA theory
    fname = os.path.join(os.path.abspath('nasa_data'), 'retheta_cf_theory.dat')
    df = pd.read_csv(fname, comment='#')

    plt.figure(0)
    p = plt.plot(df['retheta'], df['cf'], lw=2,
                 color=cmap[0], label='K-S theory')
    p[0].set_dashes(dashseq[0])

    fname = os.path.join(os.path.abspath('nasa_data'), 'yp_up_theory.dat')
    df = pd.read_csv(fname, comment='#')
    df['log10yp'] = np.log10(df['yp'])

    plt.figure(1)
    p = plt.plot(df['log10yp'], df['up'], lw=2,
                 color=cmap[0], label='Coles theory')
    p[0].set_dashes(dashseq[0])

    # ========================================================================
    # NASA CFL3D
    fname = os.path.join(os.path.abspath('nasa_data'), 'retheta_cf_cfl3d.dat')
    df = pd.read_csv(fname, comment='#')

    plt.figure(0)
    p = plt.plot(df['retheta'], df['cf'], lw=2,
                 color=cmap[1], label='NASA CFL3D')
    p[0].set_dashes(dashseq[1])

    fname = os.path.join(os.path.abspath('nasa_data'), 'yp_up_cfl3d.dat')
    df = pd.read_csv(fname, comment='#')

    plt.figure(1)
    p = plt.plot(df['log10yp'], df['up'], lw=2,
                 color=cmap[1], label='NASA CFL3D')
    p[0].set_dashes(dashseq[1])

    # ======================================================================
    # Nalu output
    fname = os.path.join(os.path.abspath('545x385/results'), 'wall_coeffs.dat')
    df = pd.read_csv(fname)

    plt.figure(0)
    p = plt.plot(df['retheta'], df['cf'], lw=2,
                 color=cmap[2], label='Nalu')
    p[0].set_dashes(dashseq[2])

    fname = os.path.join(os.path.abspath('545x385/results'), 'yp_up.dat')
    df = pd.read_csv(fname)
    df['log10yp'] = np.log10(df['yp'])

    plt.figure(1)
    p = plt.plot(df['log10yp'], df['up'], lw=2,
                 color=cmap[2], label='Nalu')
    p[0].set_dashes(dashseq[2])

    # ======================================================================
    # Format the plots
    plt.figure(0)
    ax = plt.gca()
    plt.xlabel(r"$Re_\theta$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$C_f$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    ax.set_xlim([4000, 14000])
    ax.set_ylim([0.002, 0.004])
    legend = ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig('retheta_cf.pdf', format='pdf')
    plt.savefig('retheta_cf.png', format='png')

    plt.figure(1)
    ax = plt.gca()
    plt.xlabel(r"$\log_{10}(y^+)$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$u^+$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    legend = ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig('yp_up.pdf', format='pdf')
    plt.savefig('yp_up.png', format='png')

    if args.show:
        plt.show()
