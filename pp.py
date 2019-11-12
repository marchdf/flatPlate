#!/usr/bin/env python3

# ========================================================================
#
# Imports
#
# ========================================================================
import argparse
import re
import os
import numpy as np
import pandas as pd
import netCDF4 as nc
import glob
import yaml
import scipy.integrate as spi


# ========================================================================
#
# Function definitions
#
# ========================================================================
def parse_ic(fname):
    """Parse the Nalu yaml input file for the initial conditions"""
    with open(fname, "r") as stream:
        try:
            dat = yaml.load(stream)
            u0 = float(
                dat["realms"][0]["initial_conditions"][0]["value"]["velocity"][0]
            )
            rho0 = float(
                dat["realms"][0]["material_properties"]["specifications"][0]["value"]
            )
            mu = float(
                dat["realms"][0]["material_properties"]["specifications"][1]["value"]
            )

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
        dat = nc.Dataset(ename)

        ssn = ["%s" % nc.chartostring(ss) for ss in dat.variables["ss_names"][:]]
        vn = ["%s" % nc.chartostring(nn) for nn in dat.variables["name_nod_var"][:]]
        idx_ss = ssn.index("bottomwall")
        idx_tw = vn.index("tau_wall")
        idx_p = vn.index("pressure")

        # Get wall x coordinates
        try:
            wall_elem_idx = dat.variables["elem_ss{0:d}".format(idx_ss + 1)][:] - 1
        except:
            continue
        wall_connect1 = dat.variables["connect1"][wall_elem_idx].flatten() - 1
        wall_coordx = dat.variables["coordx"][wall_connect1]
        wall_coordz = dat.variables["coordz"][wall_connect1]

        actual_idx = np.where(wall_coordz < 1e-14)
        wall_x = wall_coordx[actual_idx]
        wall_connect1 = wall_connect1[actual_idx]

        wall_x, unique_idx = np.unique(wall_x, return_index=True)
        wall_connect1 = wall_connect1[unique_idx]

        # Get wall variables
        tau_wall = dat.variables["vals_nod_var{0:d}".format(idx_tw + 1)][
            -1, wall_connect1
        ]
        pressure = dat.variables["vals_nod_var{0:d}".format(idx_p + 1)][
            -1, wall_connect1
        ]

        # To dataframe
        df = pd.DataFrame(
            data=np.vstack((wall_x, tau_wall, pressure)).T,
            columns=["x", "tau_wall", "pressure"],
        )
        lst.append(df)

    # Save
    dfw = pd.concat(lst, ignore_index=True)
    dfw = dfw.sort_values(by=["x"])
    dfw = dfw.drop_duplicates(subset="x")
    return dfw.reset_index(drop=True)


# ========================================================================
def get_fields_ss(enames, fields, ssname):
    """Get fields on the sideset of the Exodus file."""

    lst = []

    for ename in enames:
        dat = nc.Dataset(ename)

        ssn = ["%s" % nc.chartostring(ss) for ss in dat.variables["ss_names"][:]]
        vn = ["%s" % nc.chartostring(nn) for nn in dat.variables["name_nod_var"][:]]
        idx = ssn.index(ssname)

        try:
            elem_idx = dat.variables["elem_ss{0:d}".format(idx + 1)][:] - 1
        except KeyError:
            continue

        connect1 = dat.variables["connect1"][elem_idx].flatten() - 1
        coordx = dat.variables["coordx"][connect1]
        coordy = dat.variables["coordy"][connect1]
        coordz = dat.variables["coordz"][connect1]

        actual_idx = np.where(coordy >= 0.0)[0]
        x = coordx[actual_idx]
        z = coordz[actual_idx]
        connect1 = connect1[actual_idx]

        df = pd.DataFrame(data=np.vstack((x, z)).T, columns=["x", "z"])
        for key, val in fields.items():
            idx_field = vn.index(val)
            field_values = dat.variables["vals_nod_var{0:d}".format(idx_field + 1)][
                -1, connect1
            ]
            df[key] = field_values
        lst.append(df)

    # Save
    df = pd.concat(lst, ignore_index=True)
    df.x[np.fabs(df["x"]) < 1e-14] = 0.0  # make true zeros
    df.z[np.fabs(df["z"]) < 1e-14] = 0.0  # make true zeros
    df = df.sort_values(by=["x", "z"])
    df = df.drop_duplicates(subset=["x", "z"])
    df = df[df["x"] >= 0]  # remove everything before the plate
    return df.reset_index(drop=True)


# ========================================================================
def get_profiles(enames, fields, xloc):
    """Get field profile at xloc."""

    lst = []

    for ename in enames:
        dat = nc.Dataset(ename)

        ssn = ["%s" % nc.chartostring(ss) for ss in dat.variables["ss_names"][:]]
        vn = ["%s" % nc.chartostring(nn) for nn in dat.variables["name_nod_var"][:]]
        idx = ssn.index("back")

        try:
            elem_idx = dat.variables["elem_ss{0:d}".format(idx + 1)][:] - 1
        except KeyError:
            continue

        connect1 = dat.variables["connect1"][elem_idx].flatten() - 1
        coordx = dat.variables["coordx"][connect1]
        coordy = dat.variables["coordy"][connect1]
        coordz = dat.variables["coordz"][connect1]

        actual_idx = np.where((coordy >= 0.0) & (np.fabs(coordx - xloc) < 1e-3))[0]
        if len(actual_idx) > 0:
            x = coordx[actual_idx]
            z = coordz[actual_idx]
            connect1 = connect1[actual_idx]

            df = pd.DataFrame(data=np.vstack((x, z)).T, columns=["x", "z"])
            for key, val in fields.items():
                idx_field = vn.index(val)
                field_values = dat.variables["vals_nod_var{0:d}".format(idx_field + 1)][
                    -1, connect1
                ]
                df[key] = field_values
            lst.append(df)

    # Save
    df = pd.concat(lst, ignore_index=True)
    df.x[np.fabs(df["x"]) < 1e-14] = 0.0  # make true zeros
    df.z[np.fabs(df["z"]) < 1e-14] = 0.0  # make true zeros
    df = df.sort_values(by=["x", "z"])
    df = df.drop_duplicates(subset=["x", "z"])
    return df.reset_index(drop=True)


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == "__main__":

    # ========================================================================
    # Parse arguments
    parser = argparse.ArgumentParser(description="Postprocess Nalu data")
    args = parser.parse_args()

    # ========================================================================
    # Setup
    ppdirs = ["35x25", "69x49", "137x97", "273x193", "545x385"]
    fdirs = [os.path.abspath(fdir) for fdir in ppdirs]
    rdirs = [os.path.join(fdir, "results") for fdir in fdirs]
    ocname = os.path.join(os.path.abspath("."), "coeffs.dat")
    colnames = ["N2", "h", "cf", "cp", "cd", "cl"]
    dfc = pd.DataFrame(index=ppdirs, columns=colnames)

    # ========================================================================
    # Post-process
    for ppdir, fdir, rdir in zip(ppdirs, fdirs, rdirs):
        print("Post-processing directory:", ppdir)
        fname = os.path.join(fdir, "flatPlate.dat")
        yname = os.path.join(fdir, "flatPlate.yaml")
        enames = glob.glob(os.path.join(rdir, "flatPlate.e*"))
        owname = os.path.join(rdir, "wall_coeffs.dat")
        ouname = os.path.join(rdir, "yp_up.dat")
        opname = os.path.join(rdir, "profiles.dat")

        # Derived quantities
        L = 2.0
        W = 1.0
        area = L * W
        mach = 0.2
        u0, rho0, mu = parse_ic(yname)
        ainf2 = (u0 / mach) ** 2
        dynPres = rho0 * 0.5 * u0 * u0

        # ---------------------------------------------
        # Get wall values, coefficients, etc
        dfw = get_wall_values(enames)

        # Get integrated wall values
        df = pd.read_csv(fname, delim_whitespace=True)

        # Calculate coefficients
        df["cl"] = (df["Fpz"] + df["Fvz"]) / (dynPres * area)
        df["cd"] = (df["Fpx"] + df["Fvx"]) / (dynPres * area)
        dfw["cf"] = dfw["tau_wall"] / dynPres
        dfw["cp"] = dfw["pressure"] / dynPres
        xslice = 0.97008
        cf_slice = np.interp(xslice, dfw["x"], dfw["cf"])
        cp_slice = np.interp(xslice, dfw["x"], dfw["cp"])
        print(cf_slice, cp_slice)

        # Get some profiles
        profiles = get_profiles(
            enames,
            {
                "ux": "velocity_x",
                "mut": "turbulent_viscosity",
                "F1": "sst_f_one_blending",
                "tke": "turbulent_ke",
                "sdr": "specific_dissipation_rate",
                "tw": "tau_wall",
            },
            xslice,
        )
        profiles.ux /= u0
        profiles.mut /= mu
        profiles.tke /= ainf2
        profiles.sdr *= mu / (rho0 * ainf2)
        tw = profiles.tw.max()
        profiles["uxp"] = profiles.ux * u0 / np.sqrt(tw / rho0)
        profiles["yp"] = profiles.z * np.sqrt(tw / rho0) / (mu / rho0)

        # ---------------------------------------------
        # Also calculate Re_theta and u+
        dfu = get_fields_ss(enames, {"ux": "velocity_x"}, "back")

        # Reshapes
        res = [int(num) for num in re.findall(r"\d+", ppdir)]
        nx, nz = int(dfu.shape[0] / res[1]), res[1]
        x = dfu["x"].values.reshape((nx, nz))
        z = dfu["z"].values.reshape((nx, nz))
        ux = dfu["ux"].values.reshape((nx, nz))

        # Integrate to get theta and other quantities (and save)
        # The integral bound is where u is less than 99.5% of the freestream
        # velocity
        theta = []
        for i in range(ux.shape[0]):
            idx = ux[i, :] < 0.995 * u0
            theta.append(spi.simps(ux[i, idx] / u0 * (1 - ux[i, idx] / u0), z[i, idx]))
        dfw["theta"] = theta
        dfw["retheta"] = rho0 * u0 * dfw["theta"] / mu

        # Get y+ and u+ at Re-theta = 10000 (or next best thing)
        idx = np.argmin(np.fabs(dfw["retheta"].values - 10000))
        up = ux[idx, :] / np.sqrt(dfw["tau_wall"][idx] / rho0)
        yp = z[idx, :] * np.sqrt(dfw["tau_wall"][idx] / rho0) / (mu / rho0)

        odf = pd.DataFrame(
            data=np.vstack((z[idx, :], yp, up)).T, columns=["z", "yp", "up"]
        )

        # ---------------------------------------------
        # Save
        res = [int(num) for num in re.findall(r"\d+", ppdir)]
        N2 = (res[0] - 1) * (res[1] - 1)
        h = np.sqrt(1.0 / N2)
        dfc.loc[ppdir] = [
            N2,
            h,
            cf_slice,
            cp_slice,
            df["cd"].iloc[-1],
            df["cl"].iloc[-1],
        ]
        dfw.to_csv(owname, index=False)
        odf.to_csv(ouname, index=False)
        profiles.to_csv(opname, index=False)

    # Save the coefficients in a convenient table
    dfc.to_csv(ocname)
