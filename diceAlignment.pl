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
GetOptions(	'begin|B=i'=> \$begin, 'end|E=i' => \$end,
		'shift-left|L' => \$shiftLeft, 'fix-coords|F' => \$fixCoords,
		'shift-right|R' => \$shiftRight, 'translate|T' => \$translate,
		'prefix|P=s' => \$prefix, 'sort-by-count|C' => \$sortByCount
		);

if ( scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\tperl $0 <in.fasta> <...>\n";
	$message .= "\t\t-P|--prefix <STR>\tPrefix used for output.\n";
	$message .= "\t\t-B|--begin <INT>\tBegin coordinate.\n";
	$message .= "\t\t-E|--end <INT>\t\tEnd coordinate.\n";
	$message .= "\t\t-R|--shift-right\tShift alignment right.\n";
	$message .= "\t\t-T|--translate\tTranslate data as well.\n";
	$message .= "\t\t-C|--sort-by-count\tSort data by their cluster counts.\n";
	die($message."\n");
}

if ( !defined($begin) ) {
	$begin = 0;
} else {
	$begin--;
}

if ( $fixCoords ) {
	if ( ($begin % 3) != 0 ) {
		$tmp = $begin - ($begin%3);
		if ( $tmp >= 0 ) {
			print STDERR "Moving starting coordinate up to ",($tmp+1),".\n";
			$begin = $tmp;
		}
	}
}

if ( !defined($end) ) {
	$end = 'L';
} else {
	$end--;
}

if ( !defined($prefix) ) {
	$filename = basename($ARGV[0]);
	print STDERR $filename,"\n";
	@pieces = split('\.',$filename);
	if ( scalar(@pieces) > 1 ) {
		pop(@pieces);
	}
	$prefix = join('.',@pieces);
}

open(DICED,'>',$prefix.'.diced.txt') or die("Cannot write to $prefix.diced.txt\n");
open(MINI,'>',$prefix.'.miniNT.fasta') or die("Cannot write to $prefix.mini.fasta\n");
%clusterCount = %cluster = ();
$/ = ">";
$cluster=0;
while($record = <>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = join('',@lines);
	$length = length($sequence);
	
	if ( $length == 0 ) { next; }
	if ( $end eq 'L' ) {
		$stop = $length;
	} else {
		$stop = $end;
	}
	$subLen = $stop - $begin + 1;
	if ( $fixCoords ) {
		$fixCoords = 0;
		if ( ($subLen % 3) != 0 ) {
			$tmp = $stop - ($subLen%3);
			if ( $stop > $begin ) {
				print STDERR "Moving stopping coordinate up to ",($tmp+1),".\n";
				$stop = $end = $tmp;
				$subLen = $stop - $begin + 1;
			}
		}
	}
	$pattern = $mid = lc(substr($sequence,$begin,$subLen));
	$pattern =~ tr/-.~//d;
	
	$left = $right = "";
	if ( $begin > 0 ) {
		$left = uc(substr($sequence,0,$begin));
	}

	if ( $stop < $length ) {
		$right = uc(substr($sequence,$stop+1,($length-$stop-1)));
	}

	if ( !defined($clusters{$pattern}) ) {
		$cluster++;
		$clusters{$pattern} = $cluster;
		$clusterCount{$cluster} = 1;
		$originals{$pattern} = $mid;
	} else {
		$clusterCount{$clusters{$pattern}}++;
	}
	print DICED $header,"\t",$left,"\t",$clusters{$pattern},"\t",$right,"\n";
}

@fillLens = @fillPats = ();
@patterns = keys(%originals);
if ( $shiftRight ) {
	foreach $i (0..$#patterns ) {
		$pattern = $patterns[$i];
		$fillLens[$i] = length($pattern);
		$fill = $subLen - $fillLens[$i];
		$fillPats[$i] = ( '-' x $fill ) . $pattern;
	}
} elsif ($shiftLeft ) {
	foreach $i (0..$#patterns ) {
		$pattern = $patterns[$i];
		$fillLens[$i] = length($pattern);
		$fill = $subLen - $fillLens[$i];
		$fillPats[$i] = $pattern .( '-' x $fill );
	}
} else {
	foreach $i (0..$#patterns ) {
		$fillLens[$i] = length($patterns[$i]);
		$fillPats[$i] = $originals{$patterns[$i]};
	}
}

if ( $sortByCount ) {
	@order = sort { $clusterCount{$clusters{$patterns[$b]}} <=> $clusterCount{$clusters{$patterns[$a]}} || $fillPats[$a] cmp $fillPats[$b] } 0..$#patterns;
} else {
	@order = sort { $fillLens[$a] <=> $fillLens[$b] || $fillPats[$a] cmp $fillPats[$b] } 0..$#patterns;
}

for $i (0..$#order) {
	$index = $order[$i];
	$pattern = $patterns[$index];
	$cluster = $clusters{$pattern};
	print MINI '>',$cluster,'|',$clusterCount{$cluster},"\n";
	print MINI $fillPats[$index],"\n";
}
close(MINI);
close(DICED);


if ( $translate ) {
	open(AA,'>',$prefix.'.miniAA.fasta') or die("Cannot write to $prefix.miniAA.fasta\n");
	$aaLength = 0;
	%aaFill = %ntClusterList = ();
	for $pattern ( @patterns ) {
		$aa = translate($pattern);
		if ( length($aa) > $aaLength ) {
			$aaLength = length($aa)
		}
		
		if ( !defined($aaClusterCount{$aa}) ) {
			$aaClusterCount{$aa} = $clusterCount{$clusters{$pattern}};
		} else {
			$aaClusterCount{$aa} += $clusterCount{$clusters{$pattern}};
		}
		push(@{$ntClusterList{$aa}},$clusters{$pattern});
	}

	for $aa ( keys(%aaClusterCount) ) {
		$gap = $aaLength - length($aa);
		if ( $shiftLeft ) {
			$fill = $aa . ('-'x$gap);
		} elsif ( $shiftRight ) {
			$fill = ('-'x$gap).$aa;
		} else {
			if ( $gap == 0 ) {
				$fill = $aa;
			} elsif ( $gap % 2 == 0 ) {
				$gap /= 2;
				$fill = ('-'x$gap) . $aa . ('-'x$gap);
			} else {
				$gap = int($gap/2);
				$fill = ('-'x$gap) . $aa . ('-'x($gap+1));
			}
		}
		$aaFill{$fill} = $aa;
	}

	if ( $sortByCount ) {
		@peptides = sort { $aaClusterCount{$aaFill{$b}} <=> $aaClusterCount{$aaFill{$a}} } keys(%aaFill);
	} else {
		@peptides = sort { $a cmp $b } keys(%aaFill);
	}
	for $fill ( @peptides ) {
		$pat = $aaFill{$fill};
		print AA '>',join('-',@{$ntClusterList{$pat}}),'|',$aaClusterCount{$pat},"\n",$fill,"\n";
	}
	close(AA);
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
