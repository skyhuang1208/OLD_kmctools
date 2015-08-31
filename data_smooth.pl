#! /usr/bin/perl -w
#############################################################################################
# output history within the specified blocks
# By Sky, Nov. 2014 
#############################################################################################

use strict;
$|=1;  # don't buffer output

if (scalar(@ARGV) != 2) 
{
        print "ABORT:: parameters required!\a\n";
        print "Usage: $0 <input_file> #Total numbers of the range#\n";
        print "numbers of the range must be an odd number\n";
	die;
}

my ($history_file, $N_smooth) = @ARGV;
die("N_smooth must be an odd number") if(0==$N_smooth%2);
my $N_range= ($N_smooth-1)/2;

if (! -e $history_file)	{ die("ABORT::history file <$history_file> doesn't exist!\n");}
else			{ open(IN, $history_file)   || die ("Cannot open t0.xyz file for reading\n");}

my $line= 0;
my $N_data;
my @store;
while (my $buff=<IN>){
	$line ++;

	my (@data) = split(/\s+/, $buff);

    $N_data= scalar(@data) if(1==$line);
    die("Number data incorrect at line $line\n") if($N_data != scalar(@data));
    
    for(my $i=0; $i<$N_data; $i ++){
        $store[$line-1][$i]= $data[$i];
    }
}

for(my $i=0; $i<$line; $i ++){
    my @data= (0)x$N_data;
    my $lower= $i-$N_range; $lower= 0 if($i-$N_range <0);
    my $upper= $i+$N_range; $upper= $line-1 if($i+$N_range >$line-1);

    my $count= 0;
    for(my $j= $lower; $j<= $upper; $j ++){
        $count ++;
        for(my $k=0; $k<$N_data; $k ++){
            $data[$k] += $store[$j][$k];
        }
    }
    for(my $k=0; $k<$N_data; $k ++){
        printf("%f ", $data[$k] /= $count);
    }
    print "\n";
}
