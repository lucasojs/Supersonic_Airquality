#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH -t 3-00:00
#SBATCH --mem=40000
#SBATCH -J GC-Standard
#SBATCH -p edr
#SBATCH --hint=nomultithread

if [[ $# -ne 1 ]]; then
    echo "Year must be specified"
    exit 100
fi

yr=$1
if [[ $yr -lt 2000 ]]; then
    echo "Year must be at least 2000"
    exit 101
elif [[ $yr -gt 2013 ]]; then
    echo "Year must be at most 2013"
    exit 101
fi

./setup_dir.py $yr
rc=$?
if [[ $rc -ne 0 ]]; then
    echo "Failed to run setup!"
    exit $rc
fi

if [[ -z $SLURM_CPUS_PER_TASK ]]; then
   echo "Must be in a slurm environment"
   exit 100
fi

if [[ -z $MPI_ROOT ]]; then
   source gcc.env
fi

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

log=log.run
if [[ -e $log ]]; then rm $log; fi
touch $log

last_rc=10
while [[ $last_rc -ne 0 ]]; do
    ./geos_hex 2>&1 | tee $log
    last_rc=$?
    if [[ $last_rc -ne 0 ]]; then
        # Check if this was just a screw-up
        grep libifcoremt log.run > /dev/null 2>&1
        if [[ $? -ne 0 ]]; then
            echo "Unexplained failure!"
            exit 100
        fi
    fi
done
exit $?
