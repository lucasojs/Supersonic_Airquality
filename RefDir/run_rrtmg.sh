#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 4-00:00:00
#SBATCH --mem=16000
#SBATCH -J GC-RRTMG
#SBATCH -p normal
#SBATCH --hint=nomultithread

if [[ $# -lt 1 ]]; then
   echo "Need at least one run year!"
   exit 100
fi

if [[ -z $SLURM_CPUS_PER_TASK ]]; then
   echo "Must be in a slurm environment"
   exit 1
fi

if [[ -z $MPI_ROOT ]]; then
   source gcc.env
fi

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

log=log.run
if [[ -e $log ]]; then rm $log; fi
touch $log

rc=0

for targ_yr in "$@"; do
   for mo in $( seq 1 12 ); do
      printf -v x_mo "%02d" $mo
      # Is the input ready?
      if [[ ! -e restarts/tracers/trac_rst.merra_L72.T93.${targ_yr}${x_mo}010000 ]]; then
         echo "Input not found for ${targ_yr}-${x_mo}"
         continue
      fi
      # Is there already output?
      out_file="output/trac_avg.merra_4x5_UCX.${targ_yr}${x_mo}010000"
      if [[ -f $out_file ]]; then
         # Check the size. If it's tiny, there's a problem
         n_bytes=$( wc -c $out_file | awk '{print $1}' )
         if [[ $n_bytes -le 8000 ]]; then
            echo "Output exists but is tiny. Deleting"
            mv $out_file BAD_$out_file
         else
            echo "Output already exists for ${targ_yr}-${x_mo}"
            continue
         fi
      fi
      printf " ----- STARTING RUN FOR %04d-%02d ---- " $targ_yr $mo >> $log
      echo "Date at start: $( date )" >> $log
      #./set_rrtmg_run.py $targ_yr $mo
      ./set_rundates.py $targ_yr $mo
      ./geos_2x25_hex 2>&1 | tee -a $log
      rc=$?
      echo "Date at end: $( date )" >> $log
      printf " ---- COMPLETED RUN FOR %04d-%02d ---- " $targ_yr $mo >> $log
      echo ""
   done
done
exit $rc
