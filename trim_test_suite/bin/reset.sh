#!/bin/bash

### RINRUS TRIMMING TEST SUITE - CREATED D. WAPPETT 2024 ###
### Clean up script - gets rid of files made by previous runs of the test suite ###

### Usage: scripts/reset.sh 

SCRIPT=$(realpath "$0"); SCRIPTPATH=$(dirname "$SCRIPT")
BASEDIR=$(cd $SCRIPTPATH && cd ../ && pwd)
cd $BASEDIR

find . -name "*.log" -delete
find . -name "slurm*.out" -delete
find . -name "*.slurm" -delete

cd single-residues/
for R in Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val; do
	cd $R
	for FG in SC NT CT SC_NT SC_CT NT_CT SC_NT_CT; do
		if [[ "$FG" == "SC"* && "$R" == "Gly" ]]; then
			continue
		else
			cd $FG
			rm -rf combs
			rm -rf workdir*
			cd ..
		fi
	done
	cd ..
done
cd ../

