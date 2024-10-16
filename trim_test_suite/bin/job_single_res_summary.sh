#!/bin/bash
#SBATCH --job-name=summarise
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --partition=acomputeq
#SBATCH --time=0-1:00:00


### SLURM JOB TO SUMMARISE RESULTS AFTER ALL TEST JOBS HAVE COMPLETED ###

# Details of current test run
runlab=TEST_RUN_LABEL
LOGFILE=${runlab}_singleres.log
reslist=TESTED_RESIDUES
sellist=TESTED_SUBSETS
cleanup=CLEAN_COMBS_FILES

# Remove job list needed for dependency and temporary message in log file
rm jobdep_$runlab.log
sed -i "/Testing in progress/d" $LOGFILE

# Check and log results
echo "---RESULTS---" >> $LOGFILE
echo -n "Results checked: " >> $LOGFILE
date >> $LOGFILE
echo -n "      " >> $LOGFILE
for s in $sellist; do
	hd=$(printf %-9s $s)
	echo -n "$hd" >> $LOGFILE
done
echo >> $LOGFILE

errcount=0
for R in $reslist; do
	echo -n "$R   " >> $LOGFILE
	for FG in $sellist; do
		if [[ "$FG" == "SC"* && "$R" == "Gly" ]]; then
			echo -n "N/A      " >> $LOGFILE
		else
			fails=$(tail -n +3 single-residues/$R/$FG/${FG}_${runlab}.log | awk '{ tot +=$4 } END {print tot }')
			if [[ "$fails" == "0" ]]; then
				echo -n "all ok   " >> $LOGFILE
			else
				echo -n "problem! " >> $LOGFILE
				errcount=$(( errcount + 1 ))
			fi
			if [[ "$cleanup" == "true" ]]; then
				rm -rf single-residues/$R/$FG/combs
			fi
		fi
	done
	echo >> $LOGFILE
done
echo >> $LOGFILE

echo "---OUTCOME---" >> $LOGFILE
if [[ "$errcount" == "0" ]]; then
	echo "Success! Everything is exactly as expected" >> $LOGFILE
else
	echo "Uh oh, something isn't doing what it's meant to..." >> $LOGFILE
fi
