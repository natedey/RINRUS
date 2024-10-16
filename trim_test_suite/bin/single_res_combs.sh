#!/bin/bash

### Set up required combs files ###
### Created: D Wappett Sep 2024 ###


### Define help message ###
Help()
{
	# Display Help
	echo "Create file [set].dat containing all subset atom combinations"
	echo "Syntax: 1-scripts/set_up_combs.sh [-r <residues>] [-s <atom subsets>]"
	echo "Specify multiple choices for options as comma-separated list, no spaces"
	echo "Options:"
	echo "r   Select residue(s) to test [default: all]"
	echo "s   Select atom subset(s) to test [default: all]"
	echo "    Subsets: SC,NT,CT,SC_NT,SC_CT,NT_CT,SC_NT_CT"
	echo "h   Print this information"
}


### Set option defaults ###
reslist="Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val"
sellist="SC NT CT SC_NT SC_CT NT_CT SC_NT_CT"

### Read input arguments ###
while getopts ":r:s:h:" option; do
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
BASEDIR=$(cd $SCRIPTPATH && cd ../single-residues && pwd)
cd $BASEDIR

for R in $reslist; do
	cd $BASEDIR/$R
	for FG in $sellist; do
		# Main chain sets
		if [[ "$FG" == "NT" || "$FG" == "CT" || "$FG" == "NT_CT" ]]; then
			cd $BASEDIR/$R/$FG
			if [ ! -d combs ]; then mkdir combs; else rm -f combs/*.dat; fi
			if [[ "$FG" == "NT"* && "$R" == "Pro" ]]; then
				cp $SCRIPTPATH/comb_templates/${FG}_Pro.dat combs/$FG.dat
				cp $SCRIPTPATH/comb_templates/${FG}_A_Pro.dat combs/${FG}_A.dat
			elif [[ "$R" == "Gly" ]]; then
				cp $SCRIPTPATH/comb_templates/${FG}.dat combs/$FG.dat
                                cp $SCRIPTPATH/comb_templates/${FG}_A_Gly.dat combs/${FG}_A.dat
			else
				cp $SCRIPTPATH/comb_templates/$FG.dat $SCRIPTPATH/comb_templates/${FG}_A.dat combs/
			fi
		# Skip if combo is Gly and subset with SC atoms
		elif [[ "$FG" == "SC"* && "$R" == "Gly" ]]; then
			continue
		# SC only set
		elif [[ "$FG" == "SC" ]]; then
			cd $BASEDIR/$R/$FG
                        if [ ! -d combs ]; then mkdir combs; else rm -f combs/*.dat; fi
			python3 $SCRIPTPATH/gen_SC_combs.py -res $R -file combs/SC.dat
			for A in A1 A2 A3; do
				cp combs/SC.dat combs/SC_$A.dat
				ats=$(grep "$A" $SCRIPTPATH/comb_templates/MC_labels.dat | cut -d" " -f2-)
				sed -i "s/$/ $ats/" combs/SC_$A.dat
			done
		# SC plus one side of MC
		elif [[ "$FG" == "SC_NT" || "$FG" == "SC_CT" ]]; then
			cd $BASEDIR/$R/$FG
                        if [ ! -d combs ]; then mkdir combs; else rm -f combs/*.dat; fi
			if [ -f $BASEDIR/$R/SC/combs/SC.dat ]; then
                                cp $BASEDIR/$R/SC/combs/SC.dat combs/SC.dat
                        else
                                python3 $SCRIPTPATH/gen_SC_combs.py -res $R -file combs/SC.dat
                        fi
			side=${FG: -2}
			if [[ "$FG" == "SC_NT" && "$R" == "Pro" ]]; then
				cp combs/SC.dat combs/SC_NT1.dat
				sed -i "s/$/ N/" combs/SC_NT1.dat
			else
				for j in 1 2 3; do
					cp combs/SC.dat combs/SC_${side}${j}.dat
                                	ats=$(grep "${side}${j}" $SCRIPTPATH/comb_templates/MC_labels.dat | cut -d" " -f2-)
                                	sed -i "s/$/ $ats/" combs/SC_${side}${j}.dat
                        		for A in A1 A2 A3; do
                                		cp combs/SC.dat combs/SC_${side}${j}_${A}.dat
                                		ats2=$(grep "$A" $SCRIPTPATH/comb_templates/MC_labels.dat | cut -d" " -f2-)
                                		sed -i "s/$/ $ats $ats2/" combs/SC_${side}${j}_${A}.dat
					done
                        	done
			fi
			rm combs/SC.dat
		# SC plus both MCs
		elif [[ "$FG" == "SC_NT_CT" ]]; then
			cd $BASEDIR/$R/$FG
                        if [ ! -d combs ]; then mkdir combs; else rm -f combs/*.dat; fi
			if [ -f $BASEDIR/$R/SC/combs/SC.dat ]; then 
				cp $BASEDIR/$R/SC/combs/SC.dat combs/SC.dat
			else
				python3 $SCRIPTPATH/gen_SC_combs.py -res $R -file combs/SC.dat
			fi
			if [[ "$R" == "Pro" ]]; then
				for C in CT1 CT2 CT3; do
                                        atsC=$(grep "$C" $SCRIPTPATH/comb_templates/MC_labels.dat | cut -d" " -f2-)
                                        cp combs/SC.dat combs/SC_NT1_${C}.dat
                                        sed -i "s/$/ N $atsC/" combs/SC_NT1_${C}.dat
                                        for A in A1 A2 A3; do
                                                cp combs/SC.dat combs/SC_NT1_${C}_${A}.dat
                                                atsA=$(grep "$A" $SCRIPTPATH/comb_templates/MC_labels.dat | cut -d" " -f2-)
                                                sed -i "s/$/ N $atsC $atsA/" combs/SC_NT1_${C}_${A}.dat
                                        done
                                done
			else
				for N in NT1 NT2 NT3; do
					atsN=$(grep "$N" $SCRIPTPATH/comb_templates/MC_labels.dat | cut -d" " -f2-)
					for C in CT1 CT2 CT3; do
						atsC=$(grep "$C" $SCRIPTPATH/comb_templates/MC_labels.dat | cut -d" " -f2-)
						cp combs/SC.dat combs/SC_${N}_${C}.dat
                        	        	sed -i "s/$/ $atsN $atsC/" combs/SC_${N}_${C}.dat
                        	        	for A in A1 A2 A3; do
                        	                	cp combs/SC.dat combs/SC_${N}_${C}_${A}.dat
                        	                	atsA=$(grep "$A" $SCRIPTPATH/comb_templates/MC_labels.dat | cut -d" " -f2-)
                        	                	sed -i "s/$/ $atsN $atsC $atsA/" combs/SC_${N}_${C}_${A}.dat
                        	        	done
					done
				done
			fi
			rm combs/SC.dat
		fi
	done
	cd $BASEDIR
done

