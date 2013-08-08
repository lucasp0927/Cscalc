#!/bin/bash
#PBS -N cscalc
#PBS -o ./out
#PBS -e ./err
#PBS -q ibm 
#PBS -l nodes=1:ppn=1,pmem=20gb
cd $PBS_O_WORKDIR

### ENV
## Python 2.7.3
#export PATH=/opt/python-2.7.3/bin:$PATH
export PATH=/home/twliu/working/epd-7.3-1-rh3-x86_64/bin:$PATH
#export LD_LIBRARY_PATH=/opt/python-2.7.3/lib/python2.7:/opt/python-2.7.3/lib:$LD_LIBRARY_PATH
## Lapack, Blas
#export LD_LIBRARY_PATH=/opt/lapack-3.4.1-gnu:/opt/BLAS-gnu:/usr/lib64:$LD_LIBRARY_PATH
###
echo '======================================================='
echo Working directory is $PBS_O_WORKDIR
echo "Starting on `hostname` at `date`"

if [ -n "$PBS_NODEFILE" ]; then
   if [ -f $PBS_NODEFILE ]; then
      echo "Nodes used for this job:"
      cat ${PBS_NODEFILE}
      NPROCS=`wc -l < $PBS_NODEFILE`
   fi
fi

#put your jobs here!
#python iterate.py -m d1/d1.txt &> log
#python Pulse.py d1/d11.p &> log

echo "Job Ended at `date`"
echo '======================================================='
