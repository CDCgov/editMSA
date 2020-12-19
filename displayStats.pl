#!/usr/bin/env perl
# Samuel Shepard - 2019.03.11
# Version 1.0
# Codon corrects sequences from user-supplied statistics.

use strict;
use warnings;
use Storable;
use Getopt::Long;

my $verbose;
GetOptions(
		'verbose|V' => \$verbose
	);

if ( scalar(@ARGV) != 1 ) {
	my $message 	= "Usage:\n\tperl $0 <STATS.sto> [options]\n";
	$message 	.= "\t\t--verbose|V\tMake display verbose.\n";
	die($message."\n");
}


# Obtain alignment-specific codon statistics if available
my $stats = $ARGV[0];
my %codonStats = %{retrieve($stats)} or die("Cannot open statistics file '$stats'.\n");

if ( defined($verbose) ) {
	my @keys = keys(%codonStats);
	if ( $keys[0] =~ /^\d+$/ ) {
		foreach my $p ( sort { $a <=> $b } @keys ) {
			foreach my $c ( sort { $codonStats{$p}{$b} <=> $codonStats{$p}{$a} } keys(%{$codonStats{$p}}) ) {
				print STDOUT $p,"\t",$c,"\t",$codonStats{$p}{$c},"\n";
			}
		}
	} else {
		foreach my $g ( sort(@keys) ) {
			print STDOUT $g,"\n";
			foreach my $p ( sort { $a <=> $b } keys(%{$codonStats{$g}}) ) {
				foreach my $c ( sort { $codonStats{$g}{$p}{$b} <=> $codonStats{$g}{$p}{$a} } keys(%{$codonStats{$g}{$p}}) ) {
					print STDOUT $p,"\t",$c,"\t",$codonStats{$g}{$p}{$c},"\n";
				}
			}
			print STDOUT "\n";
		}
	}
} else {
	my $maxC = 0;
	my @keys = keys(%codonStats);
	if ( $keys[0] =~ /^\d+$/ ) {
		@keys = sort { $a <=> $b } @keys;
		print STDOUT "Number of Codons\tMaximum Codon Depth\n";
		foreach my $p ( @keys ) {
			foreach my $c ( sort { $codonStats{$p}{$b} <=> $codonStats{$p}{$a} } keys(%{$codonStats{$p}}) ) {
				if ( $codonStats{$p}{$c} > $maxC ) {
					$maxC = $codonStats{$p}{$c};
				}	
			}
		}
		print STDOUT $#keys,"\t",$maxC,"\n";
	} else {
		print STDOUT "Group\tNumber of Codons\tMaximum Codon Depth\n";
		foreach my $g ( sort(@keys) ) {
			my @positions = sort { $a <=> $b } keys(%{$codonStats{$g}});
			$maxC = 0;

			foreach my $p ( @positions ) {
				foreach my $c ( sort { $codonStats{$g}{$p}{$b} <=> $codonStats{$g}{$p}{$a} } keys(%{$codonStats{$g}{$p}}) ) {
					if ( $codonStats{$g}{$p}{$c} > $maxC ) {
						$maxC = $codonStats{$g}{$p}{$c};
					}	
				}
			}
			print STDOUT $g,"\t",$#positions,"\t",$maxC,"\n";
		}
	}
}
