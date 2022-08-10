#!/usr/bin/env perl
# Sam Shepard - 2014

%gc = (
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G'     # Glycine
);

use Getopt::Long;
GetOptions( 'begin|B=i' => \$begin, 'end|E=i' => \$end, 'shift-right|R' => \$right, 'translate|T' => \$translate );

if ( scalar(@ARGV) != 1 ) {
    $message = "Usage:\n\tperl $0 <in.fasta> <...>\n";
    $message .= "\t\t-B|--begin <INT>\tBegin coordinate.\n";
    $message .= "\t\t-E|--end <INT>\t\tEnd coordinate.\n";
    $message .= "\t\t-R|--shift-right\tShift alignment right.\n";
    die( $message . "\n" );
}

if ( !defined($begin) ) {
    $begin = 0;
} else {
    $begin--;
}

if ( !defined($end) ) {
    $end = 'L';
} else {
    $end--;
}

@fillSeqs = @fillLens = ();
$/        = ">";
$i        = 0;
while ( $record = <> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/, $record );
    $header   = shift(@lines);
    $sequence = uc( join( '', @lines ) );
    $length   = length($sequence);

    if ( $length == 0 ) {
        next;
    }

    #	print '>',$header,"\n";

    if ( $end eq 'L' ) {
        $stop = $length;
    } else {
        $stop = $end;
    }

    $subLen = $stop - $begin + 1;
    $subSeq = lc( substr( $sequence, $begin, $subLen ) );
    $subSeq =~ tr/-.~//d;
    $fill = $subLen - length($subSeq);

    if ($translateSort) {
        $fillSeqs[$i] = translate($subSeq);
    } else {
        $fillSeqs[$i] = $subSeq;
    }
    $fillLens[$i] = length($subSeq);

    if ($right) {
        $subSeq = ( '-' x $fill ) . $subSeq;
    } else {
        $subSeq = $subSeq . ( '-' x $fill );
    }
    substr( $sequence, $begin, $subLen ) = $subSeq;

    #	print $sequence,"\n";

    @records[$i] = '>' . $header . "\n" . $sequence, "\n";
    $i++;
}
$N     = $i;
@order = sort { $fillLens[$a] <=> $fillLens[$b] || $fillSeqs[$a] cmp $fillSeqs[$b] } 0 .. $#fillSeqs;

for ( $i = 0; $i < $N; $i++ ) {
    print $records[$order[$i]], "\n";
}

sub translate($) {
    my $nt = uc( $_[0] );
    my $aa = '';
    my $i  = 0;
    for ( $i = 0; $i < length($nt); $i += 3 ) {
        $aa .= $gc{ substr( $nt, $i, 3 ) };
    }
    return $aa;
}
