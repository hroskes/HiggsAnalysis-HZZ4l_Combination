#!/bin/bash

#SBATCH --time=72:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=shared
#SBATCH --mail-type=FAIL,TIME_LIMIT_80
#SBATCH --mail-user=usarica1@jhu.edu

cd ${SLURM_SUBMIT_DIR}
echo "SLURM job running in: " `pwd`

eval `scram runtime -sh`

echo $CMSSW_VERSION

name=$2
firstpoint=$3
lastpoint=$4

cmdadd=""
if [[ "$firstpoint" != "" ]];then
cmdadd=$cmdadd" --firstPoint "$firstpoint
name=$name"_"$firstpoint
fi
if [[ "$lastpoint" != "" ]];then
cmdadd=$cmdadd" --lastPoint "$lastpoint
name=$name"_"$lastpoint
fi


cmd="-M MultiDimFit "$1" --algo=grid --points 200 -S 1 -t -1 -m 125 --saveNLL --saveSpecifiedNuis=all -v 3 -n "$name" "$cmdadd


combine $cmd
