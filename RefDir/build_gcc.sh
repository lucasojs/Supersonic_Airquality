#!/bin/bash

if [[ -e geos_hex ]]; then
   rm geos_hex
fi

log=$PWD/log.build
if [[ -e $log ]]; then rm $log; fi

if [[ -z $MPI_ROOT ]]; then
   source gcc.env
fi

touch $log
res=2x25

cd CodeDir
make realclean
#make -j4 UCX=yes GRID=F40 MET=ModelE CHEM=UCX TRACEBACK=yes 2>&1 | tee $log
make -j4 UCX=yes GRID=${res} RRTMG=yes MET=MERRA2 CHEM=UCX TRACEBACK=yes 2>&1 | tee $log
cd ..

if [[ ! -e CodeDir/bin/geos ]]; then
   echo "Build failed!"
   exit 1
fi
cp CodeDir/bin/geos geos_hex

exit 0
