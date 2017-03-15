#!/bin/bash

cwd=`pwd`
mpibin=/Users/mhenryde/spack/opt/spack/darwin-elcapitan-x86_64/gcc-6.3.0/openmpi-2.0.2-62xsvcb5fa4rhlndyzdxdozaf7illg3s/bin/mpiexec

declare -a fdirs=("35x25" "69x49" "137x97" "273x193" "545x385")

for fdir in "${fdirs[@]}"
do
cd "$fdir"
echo `pwd`
${mpibin} -np 1 ${HOME}/wind/NaluWindUtils/build/preprocessing/nalu_preprocess -i add_ndtw.yaml
cd "$cwd"
done
