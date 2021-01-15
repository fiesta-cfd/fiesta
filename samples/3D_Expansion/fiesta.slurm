#!/bin/bash

#SBATCH --time=5:00:00   # walltime
#SBATCH --qos=standard
#SBATCH --nodes=16   # number of nodes
#SBATCH --ntasks=64  # number of processor cores (i.e. tasks)
#SBATCH --ntasks-per-node=4   # number of tasks per node
#SBATCH --no-requeue   # do not requeue when preempted and on node failure
#SBATCH --signal=23@60  # send signal to job at [seconds] before end
#SBATCH --job-name=3dExp # job name
#SBATCH --output=slurm-%A.log

base_dir="/lustre/scratch3/turquoise/beromer/3D_Expansion"
infile="${base_dir}/fiesta.lua"

fiesta_cmd="/users/beromer/src/fiesta/build/fiesta -n 4"

module load gcc/8.3.0
module load cudatoolkit/10.2
module load openmpi/3.1.6
module load hdf5-parallel/1.8.16

mpirun -n ${SLURM_NPROCS} -wdir ${base_dir} ${fiesta_cmd} ${infile}
