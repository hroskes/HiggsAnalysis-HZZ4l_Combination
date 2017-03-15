#!/bin/bash
NMIN=0
NMAX=40
INCREMENT=5
COUNTER=$NMIN
CTP=0
fname=$1
wname="../../125/"$2

scr="scanExp.slurm.sh"
if [[ "$fname" == *"noSyst"* ]];then
	echo "Running version with no systematics"
	scr="scanExp_noSyst.slurm.sh"
fi

rm -rf $fname
mkdir -p $fname
cp $scr $fname"/"
pushd $fname

mkdir -p Logs
while [  $COUNTER -lt $NMAX ];
do
	let CTP=$COUNTER+1
	let minVar=$COUNTER*$INCREMENT
	let maxVar=$CTP*$INCREMENT
	sbatch --output="./Logs/lsflog_ScanExp_"$fname"_"$minVar"_"$maxVar".txt" --error="./Logs/lsferr_ScanExp_"$fname"_"$minVar"_"$maxVar".err" $scr $wname $fname $minVar $maxVar
	let COUNTER=COUNTER+1
done
popd
