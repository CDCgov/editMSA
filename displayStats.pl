#!/usr/bin/env perl
# Samuel Shepard - 2019.03.11
# Version 1.0
# Codon corrects sequences from user-supplied statistics.

use Storable;
use File::Basename;
use Getopt::Long;
use constant { NIL => '__NIL__' };
GetOptions(
		'verbose|S=s' => \$verbose
	);

if ( scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\tperl $0 <STATS.sto> [options]\n";
	$message .= "\t\t--verbose|V\tMake display verbose.\n";
	die($message."\n");
}


# Obtain alignment-specific codon statistics if available
my $stats = $ARGV[0];
my %codonStats = %{retrieve($stats)} or die("Cannot open statistics file '$stats'.\n");
my @keys = keys(%codonStats);
if ( $keys[0] =~ /^\d+%/ ) {
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
