
#!/usr/bin/perl -w

#script to calculate the frequency of dominant base at each position in the pfm
use strict;
use Data::Dumper;


my $pfm_dir=shift;

my @pfm_files=glob("$pfm_dir/*.pfm");
my %domseq;

for my $file (@pfm_files){
	my $tseq=&pfmString($file);	
	for (my $j=1; $j<=3; $j++){
	for(my $i=0; $i<= (length $tseq) -$j; $i ++){
		my $sub_seq=substr($tseq,$i,$j);
		$domseq{$sub_seq}+=1;
	}
	}
}
for my $k( sort { (length $a) <=> (length $b)} keys %domseq){
	my $label="mono";
	$label="di" if(length $k ==2);
	$label="tri" if(length $k ==3);
	print $k,"\t", $domseq{$k},"\t $label\n";
	
	}

sub pfmString{
my $f=shift;
my $string;
open PFM,$f;
while(<PFM>){
	
	chomp;
	my $l2=<PFM>;
	my $l3=<PFM>;
	my $l4=<PFM>;
	chomp $l2; chomp $l3; chomp $l4;
	my @t1=split "\t", $_;
	my @t2=split "\t", $l2;
	my @t3=split "\t", $l3;
	my @t4=split "\t", $l4;

	for my $i (0..$#t1){
		my ($max_ini,$base_ini)=($t1[$i],"A");
		my ($max,$base);
		($t2[$i] > $max_ini)?(($max,$base)=($t2[$i],"C")):(($max,$base)=($max_ini,$base_ini));
		($t3[$i] > $max)?(($max,$base)=($t3[$i],"G")):(($max,$base)=($max,$base));
		($t4[$i] > $max)?(($max,$base)=($t4[$i],"T")):(($max,$base)=($max,$base));
		$string.=$base;
	}
}
return $string;
}
#print Dumper(%data);
