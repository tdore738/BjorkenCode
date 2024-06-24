#!/bin/bash

eosList=$1
totalJobs=$(wc -l < '../../YukariHydro/'$2)
echo $(wc -l < '../../YukariHydro/'$2)


hydroDir=/projects/jnorhos/tdore/Hydro_Bjorken/finiteMuB/FreezeOutScan/DNMR/cs2_paper

while read d; do
   cd $hydroDir
   mkdir -p EOS_$d
   cd EOS_$d
   mkdir -p cluster_output
   resultDir=$hydroDir/EOS_$d/cluster_output
   cd ../../../Code/SlurmScripts/Yukari
   outFi=$resultDir/output.o%a
   errFi=$resultDir/err.e%a
   sbatch -vv --array=1-$totalJobs -p qgp -o $outFi -e $errFi --export=ALL,eos=$d,rhoList=$2 runYukari.sh
done <$eosList

