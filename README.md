
# Flat plate verification and validation

The setup, grids, results are taken from
the
[flat plate case](https://turbmodels.larc.nasa.gov/flatplate.html). The
NASA meshes have a given size and I adapt in initial conditions to
ensure the same flow conditions. The Mach number and Reynolds number
(based on a reference length scale of 1m) are 0.2 and 5 million. The
length of the flat plate is 2m. The density of air is 1.177 kg/m^3 and
the viscosity is 1.846e-5 kg/ms. The inflow velocity is therefore
78.4197 m/s. Some other boundary conditions for SST are
specified [here](https://turbmodels.larc.nasa.gov/flatplate_sst.html).

## Meshing workflow

1. Get CGNS mesh from [here](https://turbmodels.larc.nasa.gov/flatplate_grids.html)
2. Use Pointwise to label the surfaces and set the BC
3. Use Shreyas' [near-distance-to-wall calculator](https://github.com/NaluCFD/NaluWindUtils) to get the NDTW

## Running

```
mpiexec -np 1 ./naluX -i flatPlate.i
```

