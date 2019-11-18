#! /usr/bin/perl
# Copyright (c) 2016-2017 Geng Liu , Dongli Li

use strict;
use FindBin '$Bin';

die "---------------------------------------------------------------------------------------
                                     PSSMHCpan
The usage of PSSMHCpan:

Description

  predict class I HLA-peptide binding affinity;V1.0

Contact

  liugeng liugeng\@genomics.cn
  lidongli lidongli\@genomics.cn

usage

     perl $0 <Peptides> <length> <HLA_Allele> <PSSM>
     Peptide sequences in fasta format,requair
     The length of peptide sequences,requair
     HLA allele type,requair. for example: HLA-A0201
     File contain PSSM. Format: PSSM_Path Allele  Length

sample

     perl ./PSSMHCpan-1.0.pl test/test.fa 9 HLA-A0201 database/PSSM/pssm_file.list
     perl ./PSSMHCpan-1.0.pl test/test.fa 9 HLA-A0102 database/PSSM/pssm_file.list

# Threshold for Strong binding peptides (IC50)  50.000
# Threshold for Weak binding peptides (IC50)    500.000

---------------------------------------------------------------------------------------\n" if @ARGV < 3;

my($peptide, $length, $allele, $pssm)=@ARGV;
$pssm ||= "$Bin/database/PSSM/pssm_file.list";

my $weight = "$Bin/database/WeighitMatrix.txt";

my %weight_hash;
open(FH, $weight) || die $!;
while(<FH>)
{
	chomp;
	my @tmp = split /\s+/;
	push @{$weight_hash{$tmp[0]}}, "$tmp[1]\t$tmp[2]\n";
}
close FH;

my %pssm_hash;
open PSSM, "$pssm" or die $!;
while(<PSSM>)
{
      chomp $_;
      my @tmp = split /\s+/;
      $pssm_hash{"$tmp[1]:$tmp[2]"} = $tmp[0];
}
close PSSM;

print "Method\tAllele\tID\tPeptide\tIC50\tLevel\n";

if(exists $pssm_hash{"$allele:$length"})
{
	 my @result=&ScorePeptides($peptide,$pssm_hash{"$allele:$length"},$allele,$length);
         print @result;
}
else
{
	if(exists $weight_hash{$allele})
	{
		my %sum;
		my $count = 0;
		foreach (@{$weight_hash{$allele}})
		{
			if(/(HLA-\S+)\t(\S+)/) 
			{
				last if($2 < 650);
				if(exists $pssm_hash{"$1:$length"})
				{
					my @result=&ScorePeptides($peptide,$pssm_hash{"$1:$length"},$1,$length);
					foreach my $r (@result)
					{
						chomp($r);
						my @case = split /\t/, $r;
						$sum{"$case[2]\t$case[3]"}[0] += $2 * $case[4];
						$sum{"$case[2]\t$case[3]"}[1] += $2;
					}
					$count ++;
				}
			}
			last if($count > 1);
		}
		if($count == 0)
		{
			print STDERR "PSSMHCpan could not predict $allele with the length of $length.\n";
		}
		else
		{
			for my $tag (keys %sum)
			{
				my ($id, $seq) = split /\t/, $tag;
				my $r = $sum{$tag}[0] / $sum{$tag}[1];
				my $level = "Positive";
				$level = "Negative" if($r > 500);
				print "PSSMHCpan\t$allele\t$id\t$seq\t$r\t$level\n";
			}
		}
	}
	else
	{
		print STDERR "PSSMHCpan could not predict $allele with the length of $length.\n";
	}
}

sub ScorePeptides{
my @lib = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');
my ($protein, $pssm, $allele, $length) = @_;
my %pssm;
local $/ = "\n>";
open(FH, $pssm) || die $!;
while(<FH>)
{
        chomp;
        s/^>//;
        my @tmp = split /\n/;
        die "Wrong PSSM, number of amino acid unequal to 20.\n" if($#tmp != 20);
        for my $i (1..$#tmp)
        {
                @{$pssm{$tmp[0]}{$lib[$i-1]}} = split /\s+/, $tmp[$i];
        }
}
close FH;

my @output=();
my $max = 0.8;
my $min = 0.8 * (1 - log(50000) / log(500));
local $/ = "\n>";
open(FH, $protein) || die $!;
while(<FH>)
{
        chomp;
        s/^>//;
        my ($tag, $content) = split(/\n/, $_, 2);
        my ($id) = $tag =~ /^(\S+)/;
        $content =~ s/\n//g;

        my $pssm_id = $allele." ".$length;
        if(exists $pssm{$pssm_id})
        {
                my $score = 0;
                for my $i (0..length($content)-1)
                {
                        my $char = substr($content, $i, 1);
                        $score += $pssm{$pssm_id}{$char}[$i];
                }
                $score = $score / $length;
                $score = $max if($score > $max);
                $score = $min if($score < $min);
                my $ic50 = 50000 **(($max - $score) / ($max - $min));
                my $level = "Positive";
                $level = "Negative" if($ic50 > 500);
                push(@output,"PSSMHCpan\t$allele\t$id\t$content\t$ic50\t$level\n");

        }
}
return @output;
close FH;
local $/ = "\n";
}
