#!/usr/bin/perl -w
use strict;



my @file=glob("./*log");
my $line="";
foreach my $i (50,100,200,500,1000){
	my @tmp=sort @file;
	foreach my $j (@tmp){
		if($j=~/$i\_/){
			open IN,$j;
				while(<IN>){
					$line=$_;
					next;
					#print $line
				}
			close IN;
			print STDERR $j,"\t",$line;
		}
	
	}
}
