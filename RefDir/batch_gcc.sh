#!/bin/bash
job=0

if [[ $# -gt 0 ]]; then
    afterjob=$1
else
    afterjob=0
fi
if [[ $# -gt 1 ]]; then
    slurm_base_args="$2"
else
    slurm_base_args=""
fi
if [[ $# -gt 2 ]]; then
    start_yr=$3
else
    start_yr=2000
fi

echo "Starting sequence from year $start_yr. Additional slurm arguments: $slurm_base_args"
if [[ $afterjob -gt 0 ]]; then
    afterjob=$(( $afterjob + $start_yr - 2000 ))
    echo " --> Dependency: $afterjob"
fi

echo '#!/bin/bash' > kill_script.sh
sim_name="$( basename $( dirname $PWD ) | cut -d'_' -f2 )_$( basename $PWD | cut -d'_' -f2 )"
for yr in $( seq $start_yr 2013 ); do
    job_name=GCC_${sim_name}_${yr}
    slurm_args=$slurm_base_args
    if [[ $job == 0 ]]; then
        if [[ $afterjob -gt 0 ]]; then
            slurm_args="$slurm_args --dependency=afterok:$afterjob"
            # Assume jobs were contiguous
            afterjob=$(( $afterjob + 1 ))
        fi
    else
        slurm_args="$slurm_args --dependency=afterok:$job"
        if [[ $afterjob -gt 0 ]]; then
            slurm_args="${slurm_args}:$afterjob"
            # Assume jobs were contiguous
            afterjob=$(( $afterjob + 1 ))
        fi
    fi
    job=$( sbatch -J $job_name --parsable $slurm_args run_gcc.sh $yr )
    echo "scancel $job" >> kill_script.sh
done
chmod +x kill_script.sh
