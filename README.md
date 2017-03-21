
# Flat plate verification and validation

This presents verification and validation efforts for the SST
implementation in [Nalu](https://github.com/NaluCFD/Nalu) using
NASA's
[2D Zero Pressure Gradient Flat Plate](https://turbmodels.larc.nasa.gov/flatplate.html). The
setup, grids, and NASA code results are all taken from that website.

The initial conditions are chosen to ensure the same flow
conditions. The Mach number and Reynolds number (based on a reference
length scale of 1m) are 0.2 and 5 million. The length of the flat
plate is 2m. The density of air is 1.177 kg/m^3 and the viscosity is
1.846e-5 kg/ms. The inflow velocity is therefore 78.4197 m/s. No
pressure gradient is imposed between the inlet and outlet (though the
NASA setup does indicate a small pressure drop across the domain). To
ensure as close a setup as the NASA test cases, no wall function is
used to model the SST wall boundary conditions and the BC for the SST
model are set according to
the
[NASA specifications](https://turbmodels.larc.nasa.gov/flatplate_sst.html).

## Using this repository
A.  Generating the meshes

1. Get CGNS mesh from [the NASA website](https://turbmodels.larc.nasa.gov/flatplate_grids.html)
2. Use Pointwise to label the surfaces and set the BC
3. Use Shreyas' [near-distance-to-wall calculator](https://github.com/NaluCFD/NaluWindUtils) to get the NDTW

B. Running

```
mpiexec -np 1 ./naluX -i flatPlate.i
```

## SST 

### Verification

There is good agreement between Nalu's SST implementation
and
[NASA's SST implementation](https://turbmodels.larc.nasa.gov/flatplate_sst.html).

#### Convergence of skin friction coefficient at x = 0.97
<img src="./cf.png" alt="Cf" width="400">

#### Convergence of drag coefficient
<img src="./cd.png" alt="Cd" width="400">

#### Skin friction coefficient along the plate at t = 0.5 (545x385 mesh)
<img src="./wall_cf.png" alt="wall_Cf" width="400">

### Validation
Nalu results are compared to flat plate theoretical results and NASA
code
results. See
[here](https://turbmodels.larc.nasa.gov/flatplate_val.html)
and
[this paper](https://turbmodels.larc.nasa.gov/NAS_Technical_Report_NAS-2016-01.pdf) for
validation details. There is good agreement between theoretical
results and Nalu simulation results (545x385 mesh).

#### Skin friction: theory and simulation
<img src="./retheta_cf.png" alt="retheta_cf" width="400">

#### Velocity laws: theory and simulation (at Re_theta = 10000)
<img src="./yp_up.png" alt="yp_up" width="400">


## Thanks
Thanks to Shreyas Ananthan, Ganesh Vijayakumar, and Matt Barone for
their helpful insight and input.
