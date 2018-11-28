#!/bin/bash

i=$1
echo "#PBS -N ${i}"                             > ${i}_submit.pbs
echo "#PBS -l pmem=8gb,nodes=1:ppn=4,walltime=250:00:00"      >> ${i}_submit.pbs
echo "#PBS -A yiheh_flux"                        >> ${i}_submit.pbs
echo "#PBS -q flux"                            >> ${i}_submit.pbs
echo "#PBS -M prith@umich.edu"                 >> ${i}_submit.pbs
echo "#PBS -e output/${i}.e"                          >> ${i}_submit.pbs
echo "#PBS -o output/${i}.o"                          >> ${i}_submit.pbs
echo "#PBS -m ea"                              >> ${i}_submit.pbs
echo "#PBS -V"                                 >> ${i}_submit.pbs
echo "ulimit -s unlimited"                     >> ${i}_submit.pbs
echo "module load julia/1.0.0"                 >> ${i}_submit.pbs
# echo "module load gcc/5.4.0 hdf5/1.8.16/gcc/5.4.0" >> ${i}_submit.pbs
echo "julia JuliaSEM/run.jl"                            >> ${i}_submit.pbs
qsub ${i}_submit.pbs
