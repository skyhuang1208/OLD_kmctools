#! /usr/bin/perl -w
#############################################################################################
# make .xyz files that can be read by VMD
# By Sky, sept 2014 
#############################################################################################

use strict;
$|=1;  # don't buffer output

if (scalar(@ARGV) != 3 && scalar(@ARGV) != 1) 
{
        print "ABORT:: parameters required!\a\n";
        print "Usage: $0 [<t0.xyz>||<t0.ltcp>] <history.def> <history.sol>\n";
        print "if only input t0.xyz then output t0.xyz directly\n";
	die;
}

my ($history_file, $history_file2, $history_file3) = @ARGV;
if(1==scalar(@ARGV)){
	$history_file2= $history_file;
	$history_file3= $history_file;
}

if (! -e $history_file)	{ die("ABORT::history file <$history_file> doesn't exist!\n");}
else			{ open(IN, $history_file)   || die ("Cannot open t0.xyz file for reading\n");}
if (! -e $history_file2){ die("ABORT::history file <$history_file2> doesn't exist!\n");}
else			{ open(IN2, $history_file2) || die ("Cannot open history.def file for reading\n");}
if (! -e $history_file3){ die("ABORT::history file <$history_file3> doesn't exist!\n");}
else			{ open(IN3, $history_file3) || die ("Cannot open history.sol file for reading\n");}

my ($ismovie, $isltcp)= (0)x2;
my ($min_time, $max_time, $bk_step, $ts_out);
my ($N_vac, $N_at1, $N_at2)=(0)x3;
my $exe_cal;

if(1==scalar(@ARGV)){
	print "Are you sure you want to output $history_file alone? (yes/no)\n";
	chomp( my $inn= <STDIN>);
	die("Do nothing~\n") if($inn ne "yes");
}
else{
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
}

my ($timestep, $time);
my $N_atoms;
my (@states, @x, @y, @z);

my $count= 0;
my $bk_count= -1;

if($isltcp != 1){
	open(VAC, "> vac.xyz") || die ("Cannot open vac.xyz for reading\n");
	open(ITL, "> itl.xyz") || die ("Cannot open itl.xyz for reading\n");
	#open(A01, "> a01.xyz") || die ("Cannot open a01.xyz for reading\n");
	open(A02, "> a02.xyz") || die ("Cannot open a02.xyz for reading\n");
}
else{
	open(LTC, "> $ts_out.ltcp") || die ("Cannot open $ts_out.ltcp for reading\n");
}

sub writedata($$){
	my ($timestep_, $time_)=@_;
	
	print VAC "$N_vac\ntime: $timestep_ $time_\n";
#	print A01 "$N_at1\ntime: $time_\n";
	print A02 "$N_at2\ntime: $timestep_ $time_\n";

	my $count_at2= 0;
	my $lasti;
	for(my $i=0; $i<$N_atoms; $i ++){
		print VAC "0 $x[$i] $y[$i] $z[$i]\n"  if( 0==$states[$i]);
		if(-1==$states[$i]){ print A02 "-1 $x[$i] $y[$i] $z[$i]\n"; $count_at2 ++; $lasti= $i;}
		print ITL "2 $x[$i] $y[$i] $z[$i]\n"  if( 3==$states[$i] || 2==$states[$i] || -2==$states[$i]);
	}
	for(my $i=$count_at2+1; $i<=$N_at2; $i ++){ print A02 "-1 $x[$lasti] $y[$lasti] $z[$lasti]\n"}
}

sub writeltcp($$){
	my ($timestep_, $time_)=@_;
	
	print LTC "$N_atoms\nT: $timestep_ $time_\n";
	for(my $i=0; $i<$N_atoms; $i ++){
		if(3==$states[$i] || 2==$states[$i] || -2==$states[$i]){
			print LTC "$states[$i] $x[$i] $y[$i] $z[$i] 0 0 0 0 0\n";
		}
		elsif(0==$states[$i]){
			print LTC "$states[$i] $x[$i] $y[$i] $z[$i] 0 0 0\n";
		}
		else{
			print LTC "$states[$i] $x[$i] $y[$i] $z[$i]\n";
		}
	}
}

sub exe_check($$){
	my ($timestep_, $time_)=@_;

	if(1==$ismovie){
		if(($timestep_ >= $min_time) && ($timestep_ <= $max_time)){
			$bk_count ++;

			if(0==($bk_count%$bk_step)){
				$count ++;
				die("no more than 1000 steps can be dumped to history\n") if ($count>1000);
				
				return 1;
			}
		}

		die("\nJob completed for $count snapshots in the movie.xyz file\n") if($timestep_ > $max_time);
	}
	else{
		return 1 if($timestep_==$ts_out);
		die("\nJob completed: only at timestep $timestep_ were written\n") if($timestep_ > $ts_out);
	}

	return 0;
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
print "WARNING!! More than 1 vcc in the system!! \n" if ($N_vac != 1);

if(1==scalar(@ARGV)){
	writedata(0, 0);
	die("Output only 1 snapshot: $history_file\nJob completed\n");
}
$exe_cal= exe_check(0, 0);
if(1==$exe_cal){
	if(1==$isltcp){	writeltcp(0, 0);}
	else{		writedata(0, 0);}
}

print "start reading history and output data\n";

sub read_def(){ # reading defects
	my $nlines= 0;
	my $line_in_block=0;
	
	do{
		$line_in_block ++;
       
		my $buff= <IN2>;
		chomp($buff);          # remove '\n'
		$buff =~s/^(\s*)//;    # remove space at the beginning

		if   ($line_in_block==1){ 
			$nlines= $buff;
			die("unexpected eof in reading def\n") if(eof);
		}
		elsif($line_in_block==2){
			my ($dump, $step_def, $time_def) = split(/\s+/, $buff);
			die("timestep inconsistent, (in sol) (in def): $timestep $step_def\n")	if($step_def != $timestep);
			die("time     inconsistent, (in sol) (in def): $time $time_def\n")	if($time_def != $time);
		}
		elsif(1==$exe_cal){
			my ($type, $ltcp, @dump) = split(/\s+/, $buff);
			die("defect is updating an occupied ltcp, (timestep) (ltcp) (type): $timestep $ltcp $states[$ltcp]\n") if($states[$ltcp] != 1);
			$states[$ltcp]= $type;
		}

	}while($line_in_block!=$nlines+2);
} # read_def

my $block=0;    # block number
my $nlines=0;
my $line_in_block=0;
while (my $buff=<IN3>){
	$line_in_block ++;
        
	chomp($buff);          # remove '\n'
        $buff =~s/^(\s*)//;    # remove space at the beginning

	if   ($line_in_block==1){ 
		$block ++;
		$nlines= $buff;
	}
	elsif($line_in_block==2){
		(my $dump, $timestep, $time) = split(/\s+/, $buff);
		print "\r$timestep";
		
		$exe_cal= exe_check($timestep, $time);
	 	if (1==$exe_cal){
			for(my $i=0; $i<$N_atoms; $i ++){ $states[$i]= 1;}
		}
	}
	elsif(1==$exe_cal){
		my ($ltcp) = split(/\s+/, $buff);
		$states[$ltcp]= -1;
	}

	if($line_in_block==$nlines+2){
		read_def();
		
		if(1==$exe_cal){
		
			if(1==$isltcp){
				writeltcp($timestep, $time);
			}
			else{
				writedata($timestep, $time);
			}
		}
		
		$line_in_block= 0;
	}
} # while()


close(IN);
close(IN2);
close(IN3);
print "\noutput completed!\n:) happy !!!\n";
