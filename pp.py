#!/usr/bin/env python3

# ========================================================================
#
# Imports
#
# ========================================================================
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import netCDF4 as nc
import glob
import yaml

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
def parse_ic(fname):
    """Parse the Nalu yaml input file for the initial conditions"""
    with open(fname, 'r') as stream:
        try:
            dat = yaml.load(stream)
            u0 = float(dat['realms'][0]['initial_conditions']
                       [0]['value']['velocity'][0])
            rho0 = float(dat['realms'][0]['material_properties']
                         ['specifications'][0]['value'])
            mu = float(dat['realms'][0]['material_properties']
                       ['specifications'][1]['value'])

            return u0, rho0, mu

        except yaml.YAMLError as exc:
            print(exc)


# ========================================================================
def get_wall_values(enames):
    """Get wall values from Exodus file."""
    data = np.array([]).reshape(0, 3)
    lst = []

    # Loop on variables
    for ename in enames:
        dat = nc.MFDataset(ename)

        ssn = ["%s" % nc.chartostring(ss)
               for ss in dat.variables['ss_names'][:]]
        vn = ["%s" % nc.chartostring(nn)
              for nn in dat.variables['name_nod_var'][:]]
        idx_ss = ssn.index("b'bottomwall'")
        idx_tw = vn.index("b'tau_wall'")
        idx_p = vn.index("b'pressure'")

        # Get wall x coordinates
        try:
            wall_elem_idx = dat.variables[
                'elem_ss{0:d}'.format(idx_ss + 1)][:] - 1
        except:
            continue
        wall_connect1 = dat.variables['connect1'][wall_elem_idx].flatten() - 1
        wall_coordx = dat.variables['coordx'][wall_connect1]
        wall_coordz = dat.variables['coordz'][wall_connect1]

        actual_idx = np.where(wall_coordz < 1e-14)
        wall_x = wall_coordx[actual_idx]
        wall_connect1 = wall_connect1[actual_idx]

        wall_x, unique_idx = np.unique(wall_x, return_index=True)
        wall_connect1 = wall_connect1[unique_idx]

        # Get wall variables
        tau_wall = dat.variables[
            'vals_nod_var{0:d}'.format(idx_tw + 1)][-1, wall_connect1]
        pressure = dat.variables[
            'vals_nod_var{0:d}'.format(idx_p + 1)][-1, wall_connect1]

        # To dataframe
        df = pd.DataFrame(data=np.vstack((wall_x, tau_wall, pressure)).T,
                          columns=['x', 'tau_wall', 'pressure'])
        lst.append(df)

    # Save
    dfw = pd.concat(lst, ignore_index=True)
    dfw = dfw.sort_values(by=['x'])
    dfw = dfw.drop_duplicates(subset='x')
    return dfw


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == '__main__':

    # ========================================================================
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='Postprocess Nalu data')
    args = parser.parse_args()

    # ========================================================================
    # Setup
    ppdirs = ['35x25', "69x49", "137x97", "273x193", "545x385"]
    fdirs = [os.path.abspath(fdir) for fdir in ppdirs]
    rdirs = [os.path.join(fdir, 'results') for fdir in fdirs]
    ocname = os.path.join(os.path.abspath('.'), 'coeffs.dat')
    colnames = ['N2', 'h', 'cf', 'cp', 'cd', 'cl']
    dfc = pd.DataFrame(index=ppdirs, columns=colnames)

    # ========================================================================
    # Post-process
    for ppdir, fdir, rdir in zip(ppdirs, fdirs, rdirs):
        print('Post-processing directory:', ppdir)
        fname = os.path.join(rdir, 'flatPlate.dat')
        yname = os.path.join(fdir, 'flatPlate.i')
        enames = glob.glob(os.path.join(rdir, 'flatPlate.e*'))
        owname = os.path.join(rdir, 'wall_coeffs.dat')

        # Derived quantities
        L = 2.0
        W = 1.0
        area = L * W
        u0, rho0, mu = parse_ic(yname)
        dynPres = rho0 * 0.5 * u0 * u0

        # Get wall values
        dfw = get_wall_values(enames)

        # Get integrated wall values
        df = pd.read_csv(fname, delim_whitespace=True)

        # Calculate coefficients
        df['cl'] = (df['Fpy'] + df['Fvy']) / (dynPres * area)
        df['cd'] = (df['Fpx'] + df['Fvx']) / (dynPres * area)
        dfw['cf'] = dfw['tau_wall'] / dynPres
        dfw['cp'] = dfw['pressure'] / dynPres
        xslice = 0.97008
        cf_slice = np.interp(xslice, dfw['x'], dfw['cf'])
        cp_slice = np.interp(xslice, dfw['x'], dfw['cp'])

        # Save them for later
        res = [float(s) for s in ppdir.split('x')]
        N2 = (res[0] - 1) * (res[1] - 1)
        h = np.sqrt(1. / N2)
        dfc.loc[ppdir] = [N2, h, cf_slice, cp_slice,
                          df['cd'].iloc[-1], df['cl'].iloc[-1]]
        dfw.to_csv(owname, index=False)

    # Save the coefficients in a convenient table
    dfc.to_csv(ocname)