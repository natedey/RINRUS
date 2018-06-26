#printing metadata of PDB file

pdb_input=$1
output_file=$2

grep 'HETNAM' $pdb_input >> $output_file


echo $(grep 'HEADER' $pdb_input | sed 's/^.*: //') >> $output_file
organism="Organism: "
x="$(grep 'ORGANISM_SCIENTIFIC:' $pdb_input | sed 's/^.*: //')"
organism="$organism $x"
echo $organism>> $output_file

gene="Gene: "
y="$(grep 'GENE:' $pdb_input | sed 's/^.*: //')"
gene="$gene $y"
echo $gene>> $output_file


echo $(grep 'TITLE' $pdb_input)>> $output_file