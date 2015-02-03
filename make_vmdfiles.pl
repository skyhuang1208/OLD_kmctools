#! /usr/bin/perl -w
#############################################################################################
# make .xyz files that can be read by VMD
# By Sky, sept 2014 
#############################################################################################

use strict;
$|=1;  # don't buffer output

if (scalar(@ARGV) != 2) 
{
        print "ABORT:: parameters required!\a\n";
        print "Usage: $0 [<t0.xyz>||<t0.ltcp>] <history.info>\n";
	die;
}

my ($history_file, $history_file2) = @ARGV;

if (! -e $history_file)	{ die("ABORT::history file <$history_file> doesn't exist!\n");}
else			{ open(IN, $history_file)   || die ("Cannot open t0.xyz file for reading\n");}
if (! -e $history_file2){ die("ABORT::history file <$history_file2> doesn't exist!\n");}
else			{ open(IN2, $history_file2) || die ("Cannot open history.info file for reading\n");}

my ($ismovie, $isltcp)= (0)x2;
my ($min_time, $max_time, $bk_step, $ts_out);
my ($N_vac, $N_at1, $N_at2)=(0)x3;

print "Output a movie or a snapshot? (1:movie)\n";
$ismovie= <STDIN>;
if(1==$ismovie){
	print "insert min_time\n"; 
	$min_time= <STDIN>;
	print "insert max_time\n"; 
	$max_time= <STDIN>;
	print "insert bk_steps\n";
	$bk_step = <STDIN>;
	print "start calculating...\n";
}
else{
	print "Output for calculations(.ltcp) or VMD(.xyz)? (1:ltcp)\n";
	$isltcp= <STDIN>;
	print "1 is chosen, output .ltcp file (WARNING: input file of t0 should be t0.ltcp)\n" if(1==$isltcp);
	print "insert the output timestep\n";
	$ts_out= <STDIN>;
	chop($ts_out);
	print "start calculating...\n";
}

my $N_atoms;
my (@states, @x, @y, @z);

my $count= -1;
my $count_out= 0;

if($isltcp != 1){
	open(VAC, "> vac.xyz") || die ("Cannot open vac.xyz for reading\n");
	#open(A01, "> a01.xyz") || die ("Cannot open a01.xyz for reading\n");
	open(A02, "> a02.xyz") || die ("Cannot open a02.xyz for reading\n");
}
else{
	open(LTC, "> $ts_out.ltcp") || die ("Cannot open $ts_out.ltcp for reading\n");
}

sub writedata($){
	my ($time_)=@_;
	
	print VAC "$N_vac\ntime: $time_\n";
#	print A01 "$N_at1\ntime: $time_\n";
	print A02 "$N_at2\ntime: $time_\n";
	my ($check_vac, $check_at1, $check_at2)= (0)x3;
	for(my $i=0; $i<$N_atoms; $i ++){
		if(0==$states[$i]) {print VAC "0 $x[$i] $y[$i] $z[$i]\n"; $check_vac ++;}
		if(1==$states[$i]) {$check_at1 ++;}
#		if(1==$states[$i]) {print A01 "1 $x[$i] $y[$i] $z[$i]\n"; $check_at1 ++;}
		if(-1==$states[$i]) {print A02 "-1 $x[$i] $y[$i] $z[$i]\n"; $check_at2 ++;}
		
		die("wrong states\n") if($states[$i]!=0 && $states[$i] !=1 && $states[$i] != -1);
	}

	die("vac number inconsistent\n") if($N_vac != $check_vac);
	die("at1 number inconsistent\n") if($N_at1 != $check_at1);
	die("at2 number inconsistent\n") if($N_at2 != $check_at2);
}

sub writeltcp($){
	my ($time_)=@_;
	
	print LTC "$N_atoms\nt0.ltcp\n";
	for(my $i=0; $i<$N_atoms; $i ++){
		print LTC "$states[$i] $x[$i] $y[$i] $z[$i]\n";
		
		die("wrong states\n") if($states[$i]!=0 && $states[$i] !=1 && $states[$i] != -1);
	}
}

sub output($$){
	my ($timestep_, $time_)=@_;

	if(1==$ismovie){
		if(($timestep_ >= $min_time) && ($timestep_ <= $max_time)){
			$count ++;

			if(0==($count%$bk_step)){
				$count_out ++;
				die("no more than 1000 steps can be dumped to history\n") if ($count_out>1000);

				writedata($time_);
			}
		}

		die("\nJob completed for $count snapshots in the movie.xyz file\n") if($timestep_ > $max_time);
	}
	else{
		if($timestep_==$ts_out){
			if(1==$isltcp){
				writeltcp($time_);
			}
			else{
				writedata($time_);
			}
			die("\nJob completed: only at timestep $timestep_ were written\n");
		}
	}
}

my $line=0;
while (my $buff=<IN>){
	$line ++;
	
	chomp($buff);          # remove '\n'
        $buff =~s/^(\s*)//;    # remove space at the beginning

	if   (1==$line){
		$N_atoms= $buff;
	}
	elsif($line>=3){
		my ($type, $xi, $yi, $zi) = split(/\s+/, $buff);
		$states[$line-3]= $type; 
		$x[$line-3]= $xi;
		$y[$line-3]= $yi;
		$z[$line-3]= $zi;

		$N_vac ++ if(0==$states[$line-3]);
		$N_at1 ++ if(1==$states[$line-3]);
		$N_at2 ++ if(-1==$states[$line-3]);
	}
	elsif($line>=($N_atoms+3) && ! ($buff=~ /\s+/)){
		die("error: t0.xyz has more lines than N_atoms in the first line\n");
	}
}
die("vac+at1+at2 != N_atoms\n") if(($N_vac+$N_at1+$N_at2) != $N_atoms);
print "reading t0.xyz's done\n";
print "start reading info and output data\n";

my $block=1;    # block number
my $line_in_block=0;

output(0, 0);

my ($timestep, $time)= (0)x2;
while (my $buff=<IN2>){
	$line_in_block ++;
        
	chomp($buff);          # remove '\n'
        $buff =~s/^(\s*)//;    # remove space at the beginning

	if($buff=~ /dt:/){
		die("line_in_block is $line_in_block and isnt 3, change code if it shouldnt be 3\n") if($line_in_block != 3);
		my ($dump, $dt) = split(/\s+/, $buff);
		
		$timestep ++;
		$time += $dt;

		output($timestep, $time);
		print "\r$timestep";

		$block ++;
		$line_in_block= 0;
	}
	else{
		my ($type, $id) = split(/\s+/, $buff);
		$states[$id]= $type;
	}
} # while()


close(IN);
close(IN2);
print "\noutput completed!\n";
