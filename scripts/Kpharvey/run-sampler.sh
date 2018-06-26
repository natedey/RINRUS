#created by Krysten Harvey

# this wrapper bash script creates a randomly subsampled ‘pairs’ file that is about n times as small as the original, n is parameter in command line

#!/bin/bash

input_pairs=$1
index_output=$2
divisor=$3

#find index of first read and use it to create head and tail files for sorting

index=$(grep -n '#chromsize' $input_pairs | tail -1| cut -f1 -d:)

 
let index=$index+1

tail -n +$index $input_pairs > $index_output

let index=$index-1

sed -i -e -n 1,"$index"p $input_pairs

#find total number of reads, randomly sort, extract specified proportion of reads to output file

result=$(wc -l $index_output|awk '{print $1}')

let result=$result/$divisor

sort -R $index_output -o $index_output

sed -i -e -n 1,"$result"p index_output

#recombine files
cat $index_output >> $input_pairs

cp $input_pairs $index_output.pairs.gz

