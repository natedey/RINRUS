#!/bin/bash
#SBATCH --job-name=RESIDUE_FGSET
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --partition=acomputeq
#SBATCH --time=0-4:00:00
#SBATCH --array=1-NUMFILES


### SLURM JOB FOR RUNNING TRIMMING TESTS ###

# date label being used for this round of tests
runlab=TEST_RUN_LABEL

MAINDIR=$(pwd)
RESDIR=$(cd ../ && pwd)
echo "${SLURM_JOB_ID}" >> joblist_${runlab}.log

# Each combs file for a residue/atom subset selection is run separately within this array job
fl=$(ls combs/*.dat | head -${SLURM_ARRAY_TASK_ID} | tail -1)
fn=${fl:6:-4}

# Reset working directory so nothing old raises errors
mkdir workdir_$runlab
mkdir workdir_$runlab/$fn
cp $RESDIR/A?A.pdb workdir_$runlab/$fn/

# Start error counter
errs=1

Ncomb=$(cat $fl | wc -l)
cd workdir_$runlab/$fn/

# Run trimming procedure on each combination of atoms, compare output to expected model (and save/label any unexpected models created)
for ((j=1;j<=$Ncomb;j++)); do
        ats=$(head -$j $MAINDIR/$fl | tail -1)
        echo "B  501  1  O2 C1 O1 C2 C4 O3 C3 C2 C4 O5 O6 H2 H41 H42 H43 H33 H1 H21 H22 H23 H5" > res_atoms.dat
        echo "B  47   1  $ats" >> res_atoms.dat
        python3 ~/git/RINRUS/bin/rinrus_trim2_pdb.py -pdb A?A.pdb -s "B:501" -model "max" > /dev/null
	if cmp -s res_?.pdb $MAINDIR/FGSET_expected_noH.pdb
        then
                echo "Works: $ats" >> ALLGOOD.log
        else
                echo -n "Error ${errs}: $ats " >> UHOH.log
		if [[ "$errs" -gt 1 ]]; then
			for ((k=1;k<=$errs;k++)); do
				if cmp -s res_?.pdb wrongmod_${k}.pdb; then
					echo "> makes wrongmod_${k}" >> UHOH.log
					unique="no"
					break
				fi
			done
			if [[ "$unique" != "no" ]]; then
				echo "> makes wrongmod_${errs}" >> UHOH.log
				mv res_?.pdb wrongmod_${errs}.pdb
				mv res_?_atom_info.dat wrongmod_${errs}_atom_info.dat
				errs=$(( errs + 1 ))
			fi
		else
			echo "> makes wrongmod_${errs}" >> UHOH.log
                	mv res_?.pdb wrongmod_${errs}.pdb
                	mv res_?_atom_info.dat wrongmod_${errs}_atom_info.dat
                	errs=$(( errs + 1 ))
		fi
        fi
done

# Log how many combinations make expected model and how many do not
Nworked=$(printf %-7s $(cat ALLGOOD.log | wc -l))
Nerr=$(printf %-7s $(cat UHOH.log | wc -l))
cd ../..
padname=$(printf %-15s $fn)
padtot=$(printf %-7s $Ncomb)
echo "$padname $padtot $Nworked $Nerr" >> FGSET_${runlab}.log

# Remove file of combinations that worked to save space
rm workdir_$runlab/$fn/ALLGOOD.log

# Log time job finished
fd=$(date)
echo "$padname | JOBID ${SLURM_JOB_ID} | finished: $fd" >> ../RESIDUE_${runlab}.log
