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

my @data;
my $line_in_block= 0;
my $nlines;
my ($dump, $timestep, $totaltime)= (0)x3;
print "Dump begins!\n";
while (my $buff=<IN>){
	$line_in_block ++;

	if(1==$line_in_block){
		$data[1]= $buff;

		chomp($buff);          # remove '\n'
		$buff =~s/^(\s*)//;    # remove space at the beginning
		$nlines= $buff;
	}
		
	($dump, $timestep, $totaltime) = split(/\s+/, $buff) if(2==$line_in_block);
	
	$data[$line_in_block]= $buff if($line_in_block != 1 && $timestep>=$bg && $timestep<=$ed);

	if($line_in_block==$nlines+2){
		if($timestep>=$bg && $timestep<=$ed){
			for(my $i=1; $i<=$nlines+2; $i ++){
				print OUT "$data[$i]";
			}
		}

		print "\r$timestep";
	 
		last if($timestep==$ed);
		
		$line_in_block= 0;
	}
}

die("\n********** Dump completed!! **********\n")
