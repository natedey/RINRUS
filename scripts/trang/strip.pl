#!/usr/bin/perl

use File::Copy;
use File::Basename;

use strict;
use warnings;

my $args = @ARGV;
if ($args < 2)
{
    print "Usage: filename_to_process node1,node1,...,node_n\n" .
           "Example: pdb1v0y_h.sif A:340:_:SER\n" .
           "Example: pdb1v0y_h.sif A:340:_:SER,A:334:_:LYS\n"; 
    exit -1;
}

my($filename, $path, $extension) = fileparse($ARGV[0], qr/\.[^.]*/);
my @select_nodes = split(',',$ARGV[1]);

copy($ARGV[0],$path."tempfile") or die "Copy failed: $!\n";

print "" .
      "--------------------------------------------------------------------------------\n" .
      "Processing file ($filename$extension) for nodes ($ARGV[1]):\n" .
      "--------------------------------------------------------------------------------\n";


my @nodes_set;
my @edges_list;
my @selection;
my @csv_output;
my @selected_nodes;
my @nodes;


open FILE, "<$path"."tempfile" or die "Error: Could not open temporary file for processing.\n";
while (my $line = <FILE>) 
{
    foreach (@select_nodes)
    {
        # if the line read from file doesn't contain
        # combi:all_all, proceed
        if (index($line,'combi:all_all') == -1)
        {        
            # if the line does contain the node from selected nodes
            if (index($line, $_) != -1) 
            {
                chomp($line); # removes \n
                push (@nodes_set, $line);
            } 
        }                    
    }
}
close (FILE);

unlink($path."tempfile");

foreach (@select_nodes)
{
    @nodes = split(/:/,$_);
    push(@selected_nodes,$nodes[1]);      
}

print "\n" .
      "--------------------------------------------------------------------------------\n" .
      "List of first neighbors for nodes ($ARGV[1]):\n" .
      "--------------------------------------------------------------------------------\n";
      
foreach (@nodes_set)
{
    print "$_\n";
}

# regular expression here instead of a split?
my $source;
my $edge;
my $target;
my $counter = 0;

#--------------------------------------------
# Gets the different types of edges and loads
# into an array

# Ignores combi:all_all
#---------------------------------------------

foreach (@nodes_set)
{
    $source = $edge = $target = "";
    ($source, $edge, $target) = split(/ /,$_);
    
    if ($edge ne "combi:all_all")
    {    
        if (!($edge ~~ @edges_list))
        {
            push(@edges_list, $edge);
            $counter++;
        }
    }
}

print "\n" .
      "--------------------------------------------------------------------------------\n" .
      "List of distinct edge types for nodes ($ARGV[1]):\n" .
      "--------------------------------------------------------------------------------\n";

foreach (@edges_list)
{
    print "$_\n";
}

print "\n" .
      "--------------------------------------------------------------------------------\n" .
      "Generating the output file for selected nodes ($ARGV[1]):\n" .
      "--------------------------------------------------------------------------------\n";

sub combine;

my @sorted_selection;
my $row = 0;
my $column = 0;
my $max_rows = 0;
my $next_combo_row = 0;

for (my $j=1; $j<=$counter; $j++)
{
    foreach my $combo (combine [@edges_list],$j)
    {
       
        @selection = ();
        $row = $next_combo_row;
        
        # will probably need to reset max_rows here, to keep it orderly

        foreach my $node_selection (@nodes_set)
        {
            $source = $edge = $target = "";
            ($source, $edge, $target) = split(/ /,$node_selection);     
        
            for (my $k = 0; $k<$j; $k++)
            {
                # pushes the list of edge types from combo to the 
                # csv output, each row (k) adds a combo line
                $csv_output[$k+$next_combo_row][$column] = "@$combo[$k]";
                if ($edge eq @$combo[$k])
                {
                    # check if the source node is in our selection list
                    if (!($source ~~ @selection))
                    {
                        push(@selection,$source);
                    }
                
                    # check if the target node is in our selection list
                    if (!($target ~~ @selection))
                    {
                        push (@selection,$target);
                    }
                }             
            }
        }
        
        # advance the column to print the data beside
        # the list of edge combos
        $column++;
              
        # sort the selection list
        @sorted_selection = sort @selection;

        # for each node in the selection list, add them
        # to the appropriate row and column
        foreach (@sorted_selection)
        {         
            $csv_output[$row][$column] = "$_";
            $row++;
            $max_rows = $row if ($row > $max_rows);
        }  
        
        # to ensure we don't go past the max number
        # of excel columns (using less than max, but can go as high
        # as 4 from max)
        if ($column >= 16300)
        {
            $column=0;
            $next_combo_row = $max_rows+4;
        }
        else
        {        
            $column+=2; #for a blank column between                     
        }
    }

    # move down n number of rows
    $next_combo_row = $max_rows + 4;
    
    # return to the first column
    $column = 0;
   
}

# build the filename string
my $list_name = join('-',@selected_nodes);

# Write the file
open OUT, ">$path".$filename."-".$list_name.".csv" or die "Error: Could not generate output file ($path"."$filename"."-".$list_name.".csv)";

# Loop through the multidimensional array containing the csv output
# and join each set of columns (@$output returns the full compliment
# of columns for a given row)
foreach my $output (@csv_output)
{
    $_ ||= ' ' for @$output; #remove undef from array
    print OUT join(',',@$output)."\n";
   
}
close (OUT);

print "\n" .
      "Success ! \n";

#http://stackoverflow.com/questions/10299961/in-perl-how-can-i-generate-all-possible-combinations-of-a-list
sub combine {
    my ($list, $n) = @_;

    die "Too many combos for string size" if $n > @$list;

    return map [$_], @$list if $n <= 1;

    my @comb;

    for (my $i = 0; $i+$n <= @$list; ++$i) {
        my $val  = $list->[$i];
        my @rest = @$list[$i+1..$#$list];
        push @comb, [$val, @$_] for combine \@rest, $n-1;
    }
    
    return @comb;
}

