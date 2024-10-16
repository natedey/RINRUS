#!/bin/bash

### Script to run testing suite ###
### for each individual residue ###
### Created: D Wappett Sep 2024 ###


### Define help message ###
Help()
{
	# Display Help
	echo "Test how a single residue is trimmed by RINRUS based on the selected atoms"
	echo "Syntax: test_single_residues.sh [-r <residues>] [-s <atom subsets>] [-c <'keep'/'reuse'>]"
	echo "Specify multiple choices for options as comma-separated list, no spaces"
	echo "Options:"
	echo "r   Select residue(s) to test [default: all]"
	echo "s   Select atom subset(s) to test [default: all]"
	echo "    Subsets: SC,NT,CT,SC_NT,SC_CT,NT_CT,SC_NT_CT"
	echo "c   'keep' to keep combinations files or 'reuse' to use existing files"
	echo "h   Print this information"
}


### Set option defaults ###
reslist="Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val"
sellist="SC NT CT SC_NT SC_CT NT_CT SC_NT_CT"

### Read input arguments ###
while getopts ":r:s:c:h:" option; do
	case $option in
		r)  # choose residues
	 	    resinp=$(echo $OPTARG | sed "s/,/ /g")
		    reslist=""
		    for i in $resinp; do
			ilower=${i,,}
			reslist="$reslist ${ilower^}"
		    done
		    unset i ilower
		    ;;
		s)  # choose atom subset
		    selinp=$(echo $OPTARG | sed "s/,/ /g")
		    sellist=""
		    for i in $selinp; do
                        sellist="$sellist ${i^^}"
                    done
		    unset i
		    ;;
		c)  # options for combinations files
		    cflag=${OPTARG,,}
		    ;;
		h)  # display Help
		    Help
		    exit;;
		\?) # Invalid option
         	    echo "Error: Invalid option"
		    exit;;
	esac
done

### Set base directory for tests ###
SCRIPT=$(realpath "$0"); SCRIPTPATH=$(dirname "$SCRIPT")
BASEDIR=$(cd $SCRIPTPATH && cd ../ && pwd)
TESTDIR=$(cd $SCRIPTPATH && cd ../single-residues && pwd)
cd $BASEDIR

### Start testing procedure ###

# Set up overall summary log
runlab=$(date +"%Y-%m-%d")
runlab=$(echo "${runlab}_test")
if [[ -f "${runlab}_singleres.log" ]]; then
	let "a = $(ls -l ${runlab}*_singleres.log | wc -l) + 1"
	runlab=$(echo "${runlab}$a")
fi
LOGFILE=$BASEDIR/${runlab}_singleres.log

echo "---TESTING---" > $LOGFILE
echo -n "Tests started: " >> $LOGFILE
date >> $LOGFILE
echo "Residues tested: $reslist" >> $LOGFILE
echo "Atom subsets tested: $sellist" >> $LOGFILE
if [[ "$cflag" == "keep" ]]; then
	echo "Combinations files will NOT be deleted" >> $LOGFILE
elif [[ "$cflag" == "reuse" ]]; then
	echo "Using existing combinations files instead of generating new ones" >> $LOGFILE
else
	echo "Atom combinations files will be made as needed and subsequently deleted to save space" >> $LOGFILE
fi
echo >> $LOGFILE
echo "Testing in progress... Come back later..." >> $LOGFILE

# Set up slurm jobs for tests and submit to queue, set up individual res/atom subset logs
for R in $reslist; do
	cd $TESTDIR/$R
	echo -n "Tests started: " > ${R}_${runlab}.log
	date >> ${R}_${runlab}.log
	for FG in $sellist; do
		if [[ "$FG" == "SC"* && "$R" == "Gly" ]]; then
			continue
		else 
			cd $TESTDIR/$R/$FG
			date > ${FG}_${runlab}.log
                        # Generate combs files if not using existing ones
			if [[ "$cflag" != "reuse" ]]; then
				$SCRIPTPATH/single_res_combs.sh -r $R -s $FG
			fi
			echo "Set of combs    Total   Fine    Error  " >> ${FG}_${runlab}.log
                        NF=$(ls combs/*.dat | wc -l)
                        cp $SCRIPTPATH/job_single_res_test.sh trim_$runlab.slurm
                        sed -i "s/TEST_RUN_LABEL/$runlab/; s/RESIDUE/$R/g; s/FGSET/$FG/g; s/NUMFILES/$NF/g" trim_${runlab}.slurm
			if [[ "$FG" == "SC"* ]] && [[ "$R" == "Arg" || "$R" == "Trp" ]]; then
				sed -i "s/0-4:00:00/2-0:00:00/" trim_${runlab}.slurm
			elif [[ "$FG" == "SC"* && "$R" != "Ala" && "$R" != "Cys" && "$R" != "Ser" && "$R" != "Asp" ]]; then
				sed -i "s/0-4:00:00/1-0:00:00/" trim_${runlab}.slurm
			fi
                        jid=$(sbatch --parsable trim_$runlab.slurm)
			padname=$(printf %-15s $FG)
                        echo -n "$padname | JOBID $jid | submitted: " >> $TESTDIR/$R/${R}_${runlab}.log
                        date >> $TESTDIR/$R/${R}_${runlab}.log
                        echo -n ":$jid" >> $BASEDIR/jobdep_$runlab.log
                        cd $TESTDIR/$R
		fi
	done
	cd $TESTDIR
done

# Submit checking/summarising/cleanup job (dependency flag holds until everything else done)
cd $BASEDIR
cp $SCRIPTPATH/job_single_res_summary.sh summary_singleres_$runlab.slurm
deps=$(cat jobdep_$runlab.log)
if [[ "$cflag" == "keep" ]]; then
	clean="false"
else
	clean="true"
fi
sed -i "s/TEST_RUN_LABEL/$runlab/; s/TESTED_RESIDUES/\"$reslist\"/; s/TESTED_SUBSETS/\"$sellist\"/; s/CLEAN_COMBS_FILES/$clean/" summary_singleres_$runlab.slurm
sbatch --dependency=afterany$deps summary_singleres_$runlab.slurm
