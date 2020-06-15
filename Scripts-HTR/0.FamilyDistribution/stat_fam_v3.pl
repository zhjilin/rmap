#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use POSIX qw(strftime);

my $direct;
GetOptions(
	'd' =>\$direct
);
my $table=shift;
my $nowstring = strftime "%Y%b%d", localtime;
my $table_base=basename $table;
$table_base=~s/tsv/fam/;

my %htr_class=('RRM'=>0,'PUF'=>0,'CSD'=>0,'KH'=>0,'HEXIM'=>0,'Znf'=>0,'CCCH'=>0,'La'=>0,'THUMP'=>0,'YTH'=>0,'SAM'=>0,'Other'=>0,"RanBP"=>0,"S1"=>0);
my %known=('RRM'=>0,'PUF'=>0,'CSD'=>0,'KH'=>0,'HEXIM'=>0,'Znf'=>0,'CCCH'=>0,'La'=>0,'THUMP'=>0,'YTH'=>0,'SAM'=>0,'Other'=>0,"RanBP"=>0,"S1"=>0);
my %unid_class=('RRM'=>0,'PUF'=>0,'CSD'=>0,'KH'=>0,'HEXIM'=>0,'Znf'=>0,'CCCH'=>0,'La'=>0,'THUMP'=>0,'YTH'=>0,'SAM'=>0,'Other'=>0,"RanBP"=>0,"S1"=>0);
open my $f,"<",$table;
while(<$f>){
	my @t=split/\t/;
	if($t[6] ne "D" && ($t[2] !~/HTR/ && $t[3] !~/SELEX/ && $t[4] !~/RNACompete/ && $t[5] !~/RBNS/) ){
	$unid_class{'RRM'}+=1 if($t[1]=~/RRM/);
	$unid_class{'PUF'}+=1 if($t[1]=~/PUF/);
	$unid_class{'CSD'}+=1 if($t[1]=~/CSD/);
	$unid_class{'KH'}+=1 if($t[1]=~/KH/);
	$unid_class{'HEXIM'}+=1 if($t[1]=~/HEXIM/);
	$unid_class{'Znf'}+=1 if($t[1]=~/ZF/ && $t[1] !~/CCCH/);
	$unid_class{'CCCH'}+=1 if($t[1]=~/CCCH/);
	$unid_class{'La'}+=1 if($t[1]=~/La/);
	$unid_class{'THUMP'}+=1 if($t[1]=~/THUMP/);
	$unid_class{'SAM'}+=1 if($t[1]=~/SAM/);
	$unid_class{'YTH'}+=1 if($t[1]=~/YTH/);
	$unid_class{'S1'}+=1 if($t[1]=~/S1/);
	$unid_class{'RanBP'}+=1 if($t[1]=~/RanBP/);
	$unid_class{'Other'}+=1 if(/unknown|other/i);
	next;
	}
	else{
	$known{'RRM'}+=1 if($t[1]=~/RRM/);
	$known{'CSD'}+=1 if($t[1]=~/CSD/);
	$known{'PUF'}+=1 if($t[1]=~/PUF/);
	$known{'KH'}+=1 if($t[1]=~/KH/);
	$known{'HEXIM'}+=1 if($t[1]=~/HEXIM/);
	$known{'Znf'}+=1 if($t[1]=~/ZF/ && $t[1]!~/CCCH/);
	$known{'CCCH'}+=1 if($t[1]=~/CCCH/);
	$known{'La'}+=1 if($t[1]=~/La/);
	$known{'THUMP'}+=1 if($t[1]=~/THUMP/);
	$known{'SAM'}+=1 if($t[1]=~/SAM/);
	$known{'YTH'}+=1 if($t[1]=~/YTH/);
	$known{'S1'}+=1 if($t[1]=~/S1/);
	$known{'RanBP'}+=1 if($t[1]=~/RanBP/);
	$known{'Other'}+=1 if($t[1]=~/unknown|other/i);
	
	$htr_class{'RRM'}+=1 if($t[1]=~/RRM/ && $t[2]=~/HTR/);
	$htr_class{'CSD'}+=1 if($t[1]=~/CSD/ && $t[2]=~/HTR/);
	$htr_class{'PUF'}+=1 if($t[1]=~/PUF/ && $t[2]=~/HTR/);
	$htr_class{'KH'}+=1 if($t[1]=~/KH/&& $t[2]=~/HTR/);
	$htr_class{'HEXIM'}+=1 if($t[1]=~/HEXIM/&& $t[2]=~/HTR/);
	$htr_class{'Znf'}+=1 if($t[1]=~/ZF/ && $t[1]!~/CCCH/&& $t[2]=~/HTR/);
	$htr_class{'CCCH'}+=1 if($t[1]=~/CCCH/&& $t[2]=~/HTR/);
	$htr_class{'La'}+=1 if($t[1]=~/La/&& $t[2]=~/HTR/);
	$htr_class{'THUMP'}+=1 if($t[1]=~/THUMP/&& $t[2]=~/HTR/);
	$htr_class{'SAM'}+=1 if($t[1]=~/SAM/&& $t[2]=~/HTR/);
	$htr_class{'YTH'}+=1 if($t[1]=~/YTH/&& $t[2]=~/HTR/);
	$htr_class{'S1'}+=1 if($t[1]=~/S1/&& $t[2]=~/HTR/);
	$htr_class{'RanBP'}+=1 if($t[1]=~/RanBP/&& $t[2]=~/HTR/);
	$htr_class{'Other'}+=1 if($t[1]=~/unknown|other/i && $t[2]=~/HTR/);
	}
}
close $f;

&HashOut($table_base."htr".$nowstring, \%htr_class);
&HashOut($table_base."known".$nowstring, \%known);
&HashOut($table_base."unid".$nowstring, \%unid_class);
#print Dumper(%class);

sub HashOut{
	my $of=shift;
	my $h=shift;
open my$ot, ">./$of";
for my $k(keys %$h){
	print $ot "$k\t$h->{$k}\n";
}
close $ot;
}
