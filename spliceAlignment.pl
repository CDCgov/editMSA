#!/usr/bin/env perl
# Sam Shepard - 2016

%gc = (
	'TCA' => 'S', # Serine
	'TCC' => 'S', # Serine
	'TCG' => 'S', # Serine
	'TCT' => 'S', # Serine
	'TTC' => 'F', # Phenylalanine
	'TTT' => 'F', # Phenylalanine
	'TTA' => 'L', # Leucine
	'TTG' => 'L', # Leucine
	'TAC' => 'Y', # Tyrosine
	'TAT' => 'Y', # Tyrosine
	'TAA' => '*', # Stop
	'TAG' => '*', # Stop
	'TGC' => 'C', # Cysteine
	'TGT' => 'C', # Cysteine
	'TGA' => '*', # Stop
	'TGG' => 'W', # Tryptophan
	'CTA' => 'L', # Leucine
	'CTC' => 'L', # Leucine
	'CTG' => 'L', # Leucine
	'CTT' => 'L', # Leucine
	'CCA' => 'P', # Proline
	'CAT' => 'H', # Histidine
	'CAA' => 'Q', # Glutamine
	'CAG' => 'Q', # Glutamine
	'CGA' => 'R', # Arginine
	'CGC' => 'R', # Arginine
	'CGG' => 'R', # Arginine
	'CGT' => 'R', # Arginine
	'ATA' => 'I', # Isoleucine
	'ATC' => 'I', # Isoleucine
	'ATT' => 'I', # Isoleucine
	'ATG' => 'M', # Methionine
	'ACA' => 'T', # Threonine
	'ACC' => 'T', # Threonine
	'ACG' => 'T', # Threonine
	'ACT' => 'T', # Threonine
	'AAC' => 'N', # Asparagine
	'AAT' => 'N', # Asparagine
	'AAA' => 'K', # Lysine
	'AAG' => 'K', # Lysine
	'AGC' => 'S', # Serine
	'AGT' => 'S', # Serine
	'AGA' => 'R', # Arginine
	'AGG' => 'R', # Arginine
	'CCC' => 'P', # Proline
	'CCG' => 'P', # Proline
	'CCT' => 'P', # Proline
	'CAC' => 'H', # Histidine
	'GTA' => 'V', # Valine
	'GTC' => 'V', # Valine
	'GTG' => 'V', # Valine
	'GTT' => 'V', # Valine
	'GCA' => 'A', # Alanine
	'GCC' => 'A', # Alanine
	'GCG' => 'A', # Alanine
	'GCT' => 'A', # Alanine
	'GAC' => 'D', # Aspartic Acid
	'GAT' => 'D', # Aspartic Acid
	'GAA' => 'E', # Glutamic Acid
	'GAG' => 'E', # Glutamic Acid
	'GGA' => 'G', # Glycine
	'GGC' => 'G', # Glycine
	'GGG' => 'G', # Glycine
	'GGT' => 'G'  # Glycine
);

use File::Basename;
use Getopt::Long;
GetOptions(	'protein-seq|A=s'=> \$aaFile, 'fill-missing|F' => \$fillMissing
		);

if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <diced> <miniNT> [options]\n";
	$message .= "\t\t-A|--protein-seq <FILE>\t\tAmino acid alignment.\n";
	$message .= "\t\t-F|--fill-missing\t\tFill missing data with gaps rather than excluding.\n";
	die($message."\n");
}

open(DICE,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
open(NTS,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");

$/ = ">";
while($record = <NTS>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = join('',@lines);
	$length = length($sequence);
	
	if ( $length == 0 ) { next; }
	($cluster,$count) = split('\|',$header);
	$ntByCluster{$cluster} = $sequence;
	$fillLength = $length;
}
close(NTS);

if ( defined($aaFile) ) {
	open(AA,'<',$aaFile) or die("Cannot open $aaFile for reading.\n");
	$/ = '>';
	while($record = <NTS>) {
		chomp($record);
		@lines = split(/\r\n|\n|\r/, $record);
		$header = shift(@lines);
		$sequence = join('',@lines);
		$length = length($sequence);
		
		if ( $length == 0 ) { next; }
		($clusterList,$count) = split('\|',$header);
		@clusters = split('-',$clusterList);
		foreach $cluster ( @clusters ) {
			$aaByCluster{$cluster} = $sequence;
		}
		$fillLengthAA = $length;
	}
	close(AA);
}

$/ = "\n";
while($line = <DICE> ) {
	chomp($line);
	($header,$left,$cluster,$right) = ('','','','');
	($header,$left,$cluster,$right) = split("\t",$line);
	if ( defined($ntByCluster{$cluster}) ) {
		$mid = $ntByCluster{$cluster};
		if ( defined($aaByCluster{$cluster}) ) {
			$mid = realignByAA($fillLength,$fillLengthAA,$mid,$aaByCluster{$cluster});	
		}
		print '>',$header,"\n",$left,$mid,$right,"\n";
	} elsif( $fillMissing ) {
		print '>',$header,"\n",$left,('-'x$fillLength),$right,"\n";
	}
}
close(DICE);


sub realignByAA($$$$) {
	my ($ntL,$aaL,$ntS,$aaS) = @_;
	# rewrite code
	return $ntS;
}

sub translate($) {
	my $nt = uc($_[0]);
	my $aa = '';
	my $i = 0;
	my $codon = '';
	for($i=0;$i<length($nt);$i+=3) {
		$codon = substr($nt,$i,3);
		if ( !defined($gc{$codon}) ) {
			$aa .= '?';
		} else {
			$aa .= $gc{$codon};
		}
	}
	return $aa;
}
