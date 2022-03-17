#!/usr/bin/perl
use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;

#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen

### This script is to normalize the feature read counts based on the ERCC normal factor
### requirements:
#####1.	normal factor: the α, β  ---   α logYERCC counts – β. 
#####2. the column number of feature read counts


my %opts;
GetOptions(\%opts,"i:s","o:s","g:s","c:s","t:s","e:s","m:s","h:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{a} ||!defined $opts{b}||!defined $opts{c} ||!defined $opts{o} || defined $opts{h}) {
	die "************************************************************************
	Usage: $0.pl
	
	
	Function: Normalize feature read counts based on ERCC normal factor.
	
	Request Parameters:
	 -a The normal factor α Slope
	 -b The normal factor β Intercept
	 -c Read counts column number in the input file
	 -i the input file before normalization 
	 -o the output file after ERCC normalization
	
	Optional Restriction Parameters:
	-h Help
	
************************************************************************\n";
}


my $input=$opts{i};
my $output=$opts{o};
my $slope=$opts{a};
my $intercept=$opts{b};

my $column=(defined $opts{c})?$opts{c}:1;

open IN,"$input" or die $!;
open OUT,">$output" or die $!;
while (<IN>){
	chomp;
	my @array=(split/\t/,$_);
	my $cov=$array[$column];
	my $trans=(log($cov+1)/log(10))*$slope-$intercept;
	my $normcov=int(10**($tranf))-1;
	$array[$column]=$normcov;
	my $finalstr=join "\t",@array;
	print OUT "$finalstr\n";
}

close OUT;
close IN;
