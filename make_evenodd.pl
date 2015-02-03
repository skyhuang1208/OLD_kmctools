#! /usr/bin/perl -w
use strict;
$|=1;  # don't buffer output

my $line= 0;
my $even=-1;
my $odd=-1;
my (@evendata, @odddata);
while (my $buff=<STDIN>){
	$line ++;

	next if($line<=2);

	my ($state, $x, $y, $z) = split(/\s+/, $buff);

	if(($x+$y+$z) =~ m/^\d+$/){
		$even ++;
		$evendata[$even]= $buff;
	}
	else{
		$odd ++;
		$odddata[$odd]= $buff;
	}
}

printf ("Even: %d\n", $even+1);
printf ("Odd:  %d\n", $odd+1);

open(OUT, ">even.xyz");
open(OU2, ">odd.xyz");

for(my $i=0; $i<=$even; $i ++){
	print OUT "$evendata[$i]";
}
for(my $i=0; $i<=$odd; $i ++){
	print OU2 "$odddata[$i]";
}
