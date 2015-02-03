#! /usr/bin/perl -w
#############################################################################################
# output history within the specified blocks
# By Sky, Nov. 2014 
#############################################################################################

use strict;
$|=1;  # don't buffer output

if (scalar(@ARGV) != 4) 
{
        print "ABORT:: parameters required!\a\n";
        print "Usage: $0 <history.info> <output_file> <timestep_begin>(included) <timestep_end>(included)\n";
	die;
}

my ($history_file, $output_file, $bg, $ed) = @ARGV;

if (! -e $history_file)	{ die("ABORT::history file <$history_file> doesn't exist!\n");}
else			{ open(IN, $history_file)   || die ("Cannot open t0.xyz file for reading\n");}

if (-e $output_file)	{ die("ABORT::output file <$output_file> exists! <protection mode>\n");}
else			{ open(OUT, ">$output_file") || die ("Cannot open file for writing\n");}

my $line= -1;
my $timestep= 1;

print "Dump begins!\n";
while (my $buff=<IN>){
	$line ++;

	print OUT "$buff" if($timestep >= $bg);

	if($line%3==2){
		print "\r$timestep";
		die("\n********** Dump completed!! **********\n") if($ts_in==$ed);
		
		$timestep ++;
	}
}
