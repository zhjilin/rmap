#!/usr/perl -w
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use POSIX qw(strftime);
#update 2017 Apr29	add function to parse the stem loop IC from di nucleotid matrix


#script to calculate the information content for stem and loop separately.
#three modes: number [n], frequence [f] and dinucleotide dependency matrix [d]. 

my ($Mode,$Help);
GetOptions(
	'm:s'=>\$Mode,
	'Help!'=>\$Help
);
die "usage perl $0 --m [n|f|d] dir \n$!" if(defined $Help || @ARGV< 1);
die "mode not specified \n$!" if(not defined $Mode);
my $pfm_dir=shift;
my @pfms=glob("$pfm_dir/*.pfm");
my $now_string = strftime "%Y%b%d", localtime;
open OUT,">Bits_$Mode\_$now_string.table";
for my $pfm (@pfms){
	my @array;
	@array=&pfmLoader($pfm) if($Mode=~/n/);
	@array=&pssmLoader($pfm) if($Mode=~/f/);
	@array=&diLoader($pfm) if($Mode=~/d/);
	my $pfmbase=basename $pfm;
	foreach (@array){
		print OUT "$pfmbase\t$_\n";
	}
}

sub diLoader{
	my $f=shift;
	my $p;
	open $p,"<",$f;
	my (@freq,@IC,@Bits);
	my $count=0;
	my $motif_length;
	while(<$p>){
		chomp;
		my @l1=split/\t/;
		$motif_length=$#l1+1 if($count==0);
		push @freq,[@l1];
		$count++;
	}
	close $p;
	for(my $j=0; $j < $motif_length; $j++){
		for(my $i=0; $i< $count; $i++){
			$IC[$j]+=($freq[$i][$j]+0.0001)*log($freq[$i][$j]+0.0001)/log(2);
		}
		$Bits[$j]=4+$IC[$j];
	}

	return @Bits;
}



sub pssmLoader{
	my $file=shift;
	my $p;
	open $p,'<',$file;
	my (@freq,@IC,@Bits);
	while(<$p>){
		chomp;
		my @l1=split/\t/;
		my $line2=<$p>;
		my $line3=<$p>;
		my $line4=<$p>;
		my @l2=split/\t/,$line2;
		my @l3=split/\t/,$line3;
		my @l4=split/\t/,$line4;
		for (my $i=0; $i<=$#l1; $i++){
			$IC[$i]+=($l1[$i]+0.0001)*log($l1[$i]+0.0001)/log(2);
			$IC[$i]+=($l2[$i]+0.0001)*log($l2[$i]+0.0001)/log(2);
			$IC[$i]+=($l3[$i]+0.0001)*log($l3[$i]+0.0001)/log(2);
			$IC[$i]+=($l4[$i]+0.0001)*log($l4[$i]+0.0001)/log(2);
			$Bits[$i]=2+ $IC[$i];
		}
	}
	return @Bits;
}

sub pfmLoader{
	my $file=shift;
	my $p;
	open $p,'<',$file;
	my (@freq,@IC,@Bits);
	while(<$p>){
		chomp;
		my @l1=split/\t/;
		my $line2=<$p>;
		my $line3=<$p>;
		my $line4=<$p>;
		my @l2=split/\t/,$line2;
		my @l3=split/\t/,$line3;
		my @l4=split/\t/,$line4;
		for (my $i=0; $i<=$#l1; $i++){

			$freq[$i][0]=($l1[$i]+0.25)/($l1[$i]+$l4[$i]+$l3[$i]+$l2[$i]+1);
			$freq[$i][1]=($l2[$i]+0.25)/($l1[$i]+$l4[$i]+$l3[$i]+$l2[$i]+1);
			$freq[$i][2]=($l3[$i]+0.25)/($l1[$i]+$l4[$i]+$l3[$i]+$l2[$i]+1);
			$freq[$i][3]=($l4[$i]+0.25)/($l1[$i]+$l4[$i]+$l3[$i]+$l2[$i]+1);
			foreach my $c (0..3){
				$IC[$i]+=$freq[$i][$c]*log($freq[$i][$c])/log(2);
			}
			$Bits[$i]=2+ $IC[$i];
		}
	}
	return @Bits;
}
