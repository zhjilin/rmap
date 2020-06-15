use strict;
use Data::Dumper;
use Getopt::Long;



my ($half,$help);

GetOptions(
	'f:i' => \$half,
	'h!' => \$help
);

my $IC=shift;
my $slm=shift;
#$half ||=1; #default it will mark and print out all parts 

my (%slm,%info);

&Tabtoarray($IC,\%info);
&Tabtoarray($slm,\%slm);
#print Dumper(%info);

for my $k(keys %slm){
	my (@loop, @stem);
	if( $half eq 1 ){
		@stem=splice( @{$info{$k}},0,$slm{$k}[0]);
		if($#{$info{$k}} +1 > $slm{$k}[1]){
			@loop=splice( @{$info{$k}},0,$slm{$k}[1]);
		}
		else{
		@loop=@{$info{$k}};
		}
	}
	else{
	@loop=splice( @{$info{$k}}, $slm{$k}[0],$slm{$k}[1]);
	@stem=@{$info{$k}};
	}
	#push (@stem,splice( @{$info{$k}}, 0, $slm{$k}[0]));
	#push (@stem,splice( @{$info{$k}}, $slm{$k}[0]+$slm{$k}[1],$slm{$k}[0]));	
	
	&Tagoutput(\@loop,$k,"LOOP");
	&Tagoutput(\@stem,$k,"STEM");
}


sub Tagoutput{
	my $x=shift;
	my $name=shift;
	my $label=shift;
	my $i=1;
	foreach my $t ( @$x){
		print "$name.$i\t$t\t$label\n";
		$i++;
	}
}

sub Tabtoarray{

my $f=shift;
my $h=shift;
open SLM, $f;

while(<SLM>){
	
	chomp;
	my @t=split /\t/;
	my $id=shift @t;
	push @{$h->{$id}},@t;
}


close SLM;

}
