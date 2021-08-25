#!/usr/bin/perl -w
use strict;
use List::MoreUtils qw(uniq);
use Algorithm::Evolve::Util ':str';
use Math::Random;
use Math::Random qw(random_beta
                    random_chi_square
                    random_exponential
                    random_f
                    random_gamma
                    random_multivariate_normal
                    random_multinomial
                    random_noncentral_chi_square
                    random_noncentral_f
                    random_normal
                    random_permutation
                    random_permuted_index
                    random_uniform
                    random_poisson
                    random_uniform_integer
                    random_negative_binomial
                    random_binomial
                    random_seed_from_phrase
                    random_get_seed
                    random_set_seed_from_phrase
                    random_set_seed );

my $seqw="ATGGACAAGAAGGTGGTGCTGATCACAGGCTGCTCCTCGGGAATCGGTCTCAGCCTGGCTGTCCGGTTAGCATCTGACCCCGACAAAACATTCAAAGGTAACAGCCCACATCTACACACAATGATGGACTTTGTTGTATTTTCCGAGAGAGCTTGCACTGGGTAAGGGATCTGTGGTTAATCTGTGGTTTTAACTCCACGGAGACACTGTATAATAGTATGTAATGCTATACACTCATCAAGCACTTTATTAGAAACAGCACACTAATGCCAGGTAGGGCTTTCCTTTGCTCTCACAGCAGCTTCAGGTCATTGTGCATGCACAGATTCCACATGATGTTGGAAATATTCCTTTGAGATTCTGGTCCATGTTGACATGATCGATTGTGTTTAGTCTATGCCACAATGAGGAACCTGGCCAAGAAGGAGCGTCTTTTAGAGTGCGTGAAAGGCCTGCACAAGGACACCTTGGACATTCTCCAAATGGACGTGACTGACCGACAGTCCATTCTTGATGCAAGGGACAAGGTTGTGGAGAAGCGTGTGGACATTCTGGGTATGCCTAAGTGTCTGTGTGTTGTTAGGAAGATTGCAATTTCTGTACGTTAATGAGCTTTATCCTCACAACTTACAGTGTGTAATGCTGGTGTGGGTTGGATGGGGCCGCTGGAGCTGCAGTCCTTGGACTCCATGAGGCACATTCTGGAGGTCAACCTCTTAGGTACCATCCAGACCATTCAGGCCTTCCTACCAGAGATGAAGGCTGAGAGCCAGGGCCGCATTCTGGTCACTGGCAGCACCGGAGGGCTTCACGGTGAGACAAGCAAACGGTGGAGGAGTCTCTGAATGTCAACAGCATTCCACATGTGCAGTTTCTGTAATGATTTGATAGAAATCATGTAAACAAACCTTGGTAAGACTTGCTGCCATGCATGACATTGCCAAAATAACCCATCTCAGTGCTTTCTGTGTTCTCACAGGTCTCCCTTTTAATGAGGTGTACTGTGCCAGTAAATTTGCAATTGAGGGAGCATGTGAGAGTCTGGCTGTCCTCCTGCAACACTTCAATATCCAGTGAGTTTTGAAGCTTGCTGCATCTTGCTCCCTGAGTTTAGTGCTATTGGCTGTGTATGGTTTATTAATTTTTATATTGTTGATGGGGAGCGACTATTGATTTGATGCTTCATCCTGCAGTGTGAGTCTCATTGAGTGTGGTCCAGTCAACACTGACTTCCTGCTCAACCTGCAGAAGGCGGAGCTTGGGGATACAACTCACCAACAGGTTGATACCCGCACACTCAGCCTCTATGAAAAATACCTGCATCACTGTGGCTCAGTTTTCCAAAATGCAGCACAGGACACTGAGGACATTGTAAAGGTATGTTGTACATTGACAGCAGATCATAACTGCAGTAGTAAAATGTATTGATGAGTTTATCACAACATAAATAATATATAAAATTCAGTGTGACAGTATAAGAGCATCTCTAAATGTCACTGTTGTGTTTATTTTTGTCTCCAAGGTATTTCTAGATGCCATCCAGTCACCCAGCCCTGCATTCAGATACTTCACCAGTGGTGTCGTTCCACCTCTCACCCAACTGAAGATCACAGAGCCAGATGGCTCGCAGTGCATCCGTGCTATGAGTAAAATAATCTTCTCAGCTGAGGAACAATAA";
my $seqz="ATGGACAAGAAGGTGGTGCTGATCACAGGCTGCTCCTCGGGAATCGGTCTCAGCCTGGCTGTCCGGTTAGCATCTGACCCCGACAAAACATTCAAAGATAACAGCCCACATCTACACACAATGATGGACTTTGTTGTATTTTCCGAGAGAGCTTGCACTGGGTAAGGGATCTGTGGTTAATCTGTGGTTTTAACTCCACGGAGACACTGTATAATAGTATGTAATGCTATACACTCATCAAGCACTTTATTAGAAACAGCACACTAATGCCAGGTAGGGCTTTCCTTTGCTCTCACAGCAGCTTCAGGTCATTGTGCATGCACAGATTCCACATGATGTTGGAAATATTCCTTTGAGATTCTGGTCCATGTTGACATGATCGATTGTGTTTAGTCTATGCCACAATGAGGAACCTGGCCAAGAAGGAGCGTCTTTTAGAGTGCGTGAAAGGCCTGCACAAGGACACCTTGGACATTCTCCAAATGGACGTGACTGACCGACAGTCCATTCTTGATGCAAGGGACAAGGTTGTGGAGAAGCGTGTGGACATTCTGGGTATGCCTAAGTGTCTGTGTGTTGTTAGGAAGATTGCAATTTCTGTACGTTAATGAGCTTTATCCTCACAACTTACAGTGTGTAATGCTGGTGTGGGTTGGATGGGGCCGCTGGAGCTGCAGTCCTTGGACTCCATGAGGCACATTCTGGAGGTCAACCTCTTAGGTACCATCCAGACCATTCAGGCCTTCCTACCAGAGATGAAGGCTGAGAGCCAGGGCCGCATTCTGGTCACTGGCAGCACCGGAGGGCTTCACGGTGAGACAAGCAAACGGTGGAGGAGTCTCTGAATGTCAACAGCATTCCACATGTGCAGTTTCTGTAATGATTTGATAGAAATCATGTAAACAAACCTTGGTAAGACTTGCTGCCATGCATGACATTGCCAAAATAACCCATCTCAGTGCTTTCTGTGTTCTCACAGGTCTCCCTTTTAATGAGGTGTACTGTGCCAGTAAATTTGCAATTGAGGGAGCATGTGAGAGTCTGGCTGTCCTCCTGCAACACTTCAATATCCAGTGAGTTTTGAAGCTTGCTGCATCTTGCTCCCTGAGTTTAGTGCTATTGGCTGTGTATGGTTTATTAATTTTTATATTGTTGATGGGGAGCGACTATTGATTTGATGCTTCATCCTGCAGTGTGAGTCTCATTGAGTGTGGTCCAGTCAACACTGACTTCCTGCTCAACCTGCAGAAGGCGGAGCTTGGGGATACAACTCACCAACAGGTTGATACCCGCACACTCAGCCTCTATGAAAAATACCTGCATCACTGTGGCTCAGTTTTCCAAAATGCAGCACAGGACACTGAGGACATTGTAAAGGTATGTTGTACATTGACAGCAGATCATAACTGCAGTAGTAAAATGTATTGATGAGTTTATCACAACATAAATAATATATAAAATTCAGTGTGACAGTATAAGAGCATCTCTAAATGTCACTGTTGTGTTTATTTTTGTCTCCAAGGTATTTCTAGATGCCATCCAGTCACCCAGCCCTGCATTCAGATACTTCACCAGTGGTGTCGTTCCACCTCTCACCCAACTGAAGATCACAGAGCCAGATGGCTCGCAGTGCATCCGTGCTATGAGTAAAATAATCTTCTCAGCTGAGGAACAATAA";
my @bw=split(//,$seqw);
my @bz=split(//,$seqz);
my $cds="ATGGACAAGAAGGTGGTGCTGATCACAGGCTGCTCCTCGGGAATCGGTCTCAGCCTGGCTGTCCGGTTAGCATCTGACCCCGACAAAACATTCAAAGTCTATGCCACAATGAGGAACCTGGCCAAGAAGGAGCGTCTTTTAGAGTGCGTGAAAGGCCTGCACAAGGACACCTTGGACATTCTCCAAATGGACGTGACTGACCGACAGTCCATTCTTGATGCAAGGGACAAGGTTGTGGAGAAGCGTGTGGACATTCTGGTGTGTAATGCTGGTGTGGGTTGGATGGGGCCGCTGGAGCTGCAGTCCTTGGACTCCATGAGGCACATTCTGGAGGTCAACCTCTTAGGTACCATCCAGACCATTCAGGCCTTCCTACCAGAGATGAAGGCTGAGAGCCAGGGCCGCATTCTGGTCACTGGCAGCACCGGAGGGCTTCACGGTCTCCCTTTTAATGAGGTGTACTGTGCCAGTAAATTTGCAATTGAGGGAGCATGTGAGAGTCTGGCTGTCCTCCTGCAACACTTCAATATCCATGTGAGTCTCATTGAGTGTGGTCCAGTCAACACTGACTTCCTGCTCAACCTGCAGAAGGCGGAGCTTGGGGATACAACTCACCAACAGGTTGATACCCGCACACTCAGCCTCTATGAAAAATACCTGCATCACTGTGGCTCAGTTTTCCAAAATGCAGCACAGGACACTGAGGACATTGTAAAGGTATTTCTAGATGCCATCCAGTCACCCAGCCCTGCATTCAGATACTTCACCAGTGGTGTCGTTCCACCTCTCACCCAACTGAAGATCACAGAGCCAGATGGCTCGCAGTGCATCCGTGCTATGAGTAAAATAATCTTCTCAGCTGAGGAACAATAA";
my @cdsw=split(//,$cds);
my %code=(
"TTT"=>"F","TTC"=>"F",
"TTA"=>"L","TTG"=>"L","CTC"=>"L","CTA"=>"L","CTG"=>"L","CTT"=>"L",
"TCT"=>"S","TCC"=>"S","TCA"=>"S","TCG"=>"S","AGT"=>"S","AGC"=>"S",
"TAT"=>"Y","TAC"=>"Y",
"TAA"=>"TMC","TAG"=>"TMC","TGA"=>"TMC",
"TGT"=>"C","TGC"=>"C",
"TGG"=>"W",
"CCT"=>"P","CCC"=>"P","CCA"=>"P","CCG"=>"P",
"CAT"=>"H","CAC"=>"H",
"CAA"=>"Q","CAG"=>"Q",
"ATT"=>"I","ATC"=>"I","ATA"=>"I",
"ATG"=>"M",
"ACT"=>"T","ACC"=>"T","ACA"=>"T","ACG"=>"T",
"AAT"=>"N","AAC"=>"N",
"AAA"=>"K","AAG"=>"K",
"AGA"=>"R","AGG"=>"R","CGT"=>"R","CGC"=>"R","CGA"=>"R","CGG"=>"R",
"GTT"=>"V","GTC"=>"V","GTA"=>"V","GTG"=>"V",
"GCT"=>"A","GCC"=>"A","GCA"=>"A","GCG"=>"A",
"GAT"=>"D","GAC"=>"D",
"GAA"=>"E","GAG"=>"E",
"GGT"=>"G","GGC"=>"G","GGG"=>"G","GGA"=>"G",
);

my %base=(
"0"=>"A",
"1"=>"T",
"2"=>"C",
"3"=>"G"
);

my $G=20000;
my $mu=$ARGV[2];#1E-5;
my $N=$ARGV[0];
my $R=$ARGV[3];#1E-3;
my $i="";
my $j="";
my $x="";
my $y="";
my $tmp="";
my $tmp1="";
my $tmp2="";
my %seq=();
my @tmpseq=();
my @tmpseqw=();
my @tmpseqz=();
my @tmpseqa=();
my %ha=();
my @switch=();
my $dif=0;
my $c=0;
my $piw="";
my $piz="";
my $pia="";
my $Fst="";
my $dxy="";
my $mutation=0;
my $recom=0;
my $r=$ARGV[1];
my @uniq=();
###初始化
#print STDERR "Initialization begin...\t";
for $i (0 .. (4*$N-1)){
	$seq{0}{$i}{seq}=$seqw;
	$seq{0}{$i}{chr}="w";
}
#print STDERR "done\n";
###
$tmp1=$mu;$tmp1=~s/1E(-\d)/$1/;
$tmp2=$R;$tmp2=~s/1E(-\d)/$1/;
$tmp="Wsimulate".$tmp1.$tmp2."_".$N."_".$r.".out";
open OUT, ">$tmp";
print OUT "Pos\tRef\tAlt\tAltNo.\tAllNo.\tFre\n";
print STDERR "generation\tpiw\n";
for $i (1 .. $G){	
	#print OUT $i,"\t";
	#	print STDERR "Generation ",$i,"\n";
	@tmpseqw=();
	$mutation=0;
	for $j (0 .. (4*$N-1)){
		if($seq{$i-1}{$j}{chr} eq "w"){
			push(@tmpseqw,$seq{$i-1}{$j}{seq})
		}
	}
	#	print STDERR $#tmpseqw,"\t",$#tmpseqz,"\n";
	for $j (0 .. (4*$N-1)){
                $seq{$i-1}{$j}{seq}=();
	}
	for $j (0 .. 4*$N-1){
		$tmp=int(rand($#tmpseqw+1));
		@tmpseq=split(//,$tmpseqw[$tmp]);
		for $x (0 .. $#tmpseq){
			$y=random_poisson(1,$mu);
			if($y==0){
				$seq{$i}{$j}{seq}.=$tmpseq[$x];
			}
			if($y>=1){
				$mutation++;
				#print STDERR "MUTATION\n";
				$tmp = int(rand(4));
				$seq{$i}{$j}{seq}.=$base{$tmp};
			}
		}
	}
	##recombination
	$recom=0;
	for $j (0 .. 2*$N-1){
		$tmp1=$seq{$i}{$j}{seq};
		$tmp2=$seq{$i}{$j+2*$N}{seq};
		$dif=($#bw + 1 - str_agreement( $tmp1, $tmp2 ))/($#bw + 1);
		for $x (0 .. $#bw){
			$y=random_poisson(1,$R*(1-0.03*$dif));
			if($y>=1){
				$tmp=$x;
				$recom++;
				# print STDERR $tmp,"\n";
				$seq{$i}{$j}{seq}=substr($tmp1,0,$tmp).substr($tmp2,$tmp);
				$seq{$i}{$j+2*$N}{seq}=substr($tmp2,0,$tmp).substr($tmp1,$tmp);
				last;
			}
		}
	}

	#################3
	if($mutation==0 and $recom==0){
		for $j (0 .. 4*$N-1){
			$seq{$i}{$j}{chr}="w";
		}
	}
	elsif($mutation>0 or $recom>0){
		for $j (0 .. (4*$N-1)){
			@switch=();
			$tmp=substr($seq{$i}{$j}{seq},0,97).substr($seq{$i}{$j}{seq},393,162).substr($seq{$i}{$j}{seq},633,180).substr($seq{$i}{$j}{seq},979,94).substr($seq{$i}{$j}{seq},1193,184).substr($seq{$i}{$j}{seq},1522,156);
			#print STDERR $tmp,"\n";
			@tmpseq=split(//,$tmp);
			for $x (1 .. length($tmp)/3){
				$y=($x-1)*3;
				$tmp1=$cdsw[$y].$cdsw[$y+1].$cdsw[$y+2];
				$tmp2=$tmpseq[$y].$tmpseq[$y+1].$tmpseq[$y+2];
				if($code{$tmp2} ne $code{$tmp1}){
					if($code{$cdsw[$y].$cdsw[$y+1].$cdsw[$y+2]} ne $code{$tmpseq[$y].$cdsw[$y+1].$cdsw[$y+2]}){
						push(@switch,$tmpseq[$y]);
					}
					elsif($code{$cdsw[$y].$cdsw[$y+1].$cdsw[$y+2]} ne $code{$cdsw[$y].$tmpseq[$y+1].$cdsw[$y+2]}){
						push(@switch,$tmpseq[$y+1]);
					}
					elsif($code{$cdsw[$y].$cdsw[$y+1].$cdsw[$y+2]} ne $code{$cdsw[$y].$cdsw[$y+1].$tmpseq[$y+2]}){
						push(@switch,$tmpseq[$y+2]);
					}
				}
			}
			@tmpseq=split(//,$seq{$i}{$j}{seq});
			foreach $x (97,98,391,292,555,556,631,632,813,814,977,978,1073,1074,1191,1192,1377,1378,1520,1521){
				if($tmpseq[$x] ne $bw[$x]){
					push(@switch,$tmpseq[$x]);
				}
			}
			if(@switch){
				$seq{$i}{$j}{chr}="z";
			}
			else{
				$seq{$i}{$j}{chr}="w";
			}
		}
	}
	######recombination
	
	
	##
	if($i%10==0){
	###pi
	#print STDERR "pi caculation begin\n";
		@tmpseqw=();
		for $j (0 .. (4*$N-1)){
			if($seq{$i}{$j}{chr} eq "w"){
				if(length($seq{$i}{$j}{seq})!=1678){die}
				push(@tmpseqw,$seq{$i}{$j}{seq});
			}
		}
		#print STDERR @tmpseqz,@tmpseqw,"\n";
		$dif=0;$c=0;$piw=0;
		foreach $j (0 .. $#tmpseqw-1){
			foreach $x ($j+1 ..  $#tmpseqw){
				$c++;
				$dif+=($#bw + 1 - str_agreement( $tmpseqw[$j], $tmpseqw[$x] ))/($#bw + 1)
			}
		}
	
		$piw=$dif/$c;

		@switch=();
		print OUT "Generation\t$i\n";
		foreach $j (@tmpseqw){
			@tmpseq=split(//,$j);
			foreach $x (0..$#tmpseq){
				if($tmpseq[$x] ne $bw[$x]){
					push(@switch,$x);
					$ha{$x}=$tmpseq[$x];
				}
			}
		}
		@uniq= sort {$a<=>$b} uniq @switch;
		foreach $j (@uniq){
			$tmp=0;
			foreach $x (@tmpseqw){
				@tmpseq=split(//,$x);
				if($tmpseq[$j] ne $bw[$j]){
					$tmp++;
				}
			}
			print OUT $j,"\t",$bw[$j],"\t",$ha{$j},"\t",$tmp,"\t",$#tmpseqw+1,"\t",$tmp/($#tmpseqw+1),"\n";
		}
		
		print STDERR $i,"\t",$piw,"\n";
	}
}
close OUT;

