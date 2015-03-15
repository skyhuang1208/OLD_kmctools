#! /usr/bin/perl -w
use strict;
$|=1;  # don't buffer output

my $line=0;
foreach my $name (@ARGV){
	$line ++;

	system("rm $name") if($line%2==0);
}

