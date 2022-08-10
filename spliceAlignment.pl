#!/usr/bin/env perl
# Sam Shepard - 2016

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

use File::Basename;
use Getopt::Long;
GetOptions(
            'amino-acid-msa|A=s' => \$aaFile,
            'fill-missing|F'     => \$fillMissing,
            'prefix|P=s'         => \$prefix,
            'splice-now|S'       => \$spliceNow
);

if ( scalar(@ARGV) != 2 ) {
    $message = "Usage:\n\tperl $0 <diced> <miniNT> [options]\n";
    $message .= "\t\t-F|--fill-missing\t\tFill missing data with gaps rather than excluding.\n";
    $message .= "\t\t-A|--amino-acid-msa <FILE>\tAmino acid alignment.\n";
    $message .= "\t\t-P|--prefix <STR>\t\tPrefix used for AAtoNT file output.\n";
    $message .= "\t\t-S|--splice-now\t\t\tOutput the splice without further editing of the AAtoNT file.\n";
    die( $message . "\n" );
}

open( DICE, '<', $ARGV[0] ) or die("Cannot open $ARGV[0] for reading.\n");
open( NTS,  '<', $ARGV[1] ) or die("Cannot open $ARGV[1] for reading.\n");

if ( !defined($prefix) ) {
    $prefix = dirname( $ARGV[1] );
    $prefix .= '/' . basename( $ARGV[1], '.miniNT.fasta' );
}

$/ = ">";
while ( $record = <NTS> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/, $record );
    $header   = shift(@lines);
    $sequence = join( '', @lines );
    $length   = length($sequence);

    if ( $length == 0 ) { next; }
    ( $cluster, $count ) = split( '\|', $header );
    $ntByCluster{$cluster}    = $sequence;
    $countByCluster{$cluster} = $count;
    $fillLength               = $length;
}
close(NTS);

if ( defined($aaFile) ) {
    open( AA, '<', $aaFile ) or die("Cannot open $aaFile for reading.\n");
    $/           = '>';
    %aaByCluster = ();
    while ( $record = <AA> ) {
        chomp($record);
        @lines    = split( /\r\n|\n|\r/, $record );
        $header   = shift(@lines);
        $sequence = join( '', @lines );
        $length   = length($sequence);

        if ( $length == 0 ) { next; }
        ( $clusterList, $count ) = split( '\|', $header );
        @clusters = split( '-', $clusterList );
        foreach $cluster (@clusters) {
            $aaByCluster{$cluster} = $sequence;
        }
        $fillLengthAA = $length;
    }
    close(AA);
    open( AATONT, '>', $prefix . '.AAtoNT.fasta' ) or die("Cannot open $prefix.'.AAtoNT.fasta' for writing.\n");
    foreach $cluster ( keys(%ntByCluster) ) {
        $ntByCluster{$cluster} = realignByAA( $ntByCluster{$cluster}, $aaByCluster{$cluster} );
        print AATONT '>', $cluster, '|', $countByCluster{$cluster}, "\n", $ntByCluster{$cluster}, "\n";
    }
    close(AATONT);
    if ( !defined($spliceNow) ) {
        print STDERR "Created '$prefix.AAtoNT.fasta' for editing.\n";
        exit 0;
    }
}

$/ = "\n";
while ( $line = <DICE> ) {
    chomp($line);
    ( $header, $left, $cluster, $right ) = ( '', '', '', '' );
    ( $header, $left, $cluster, $right ) = split( "\t", $line );
    if ( defined( $ntByCluster{$cluster} ) ) {
        $mid = $ntByCluster{$cluster};
        print '>', $header, "\n", $left, $mid, $right, "\n";
    } elsif ($fillMissing) {
        print '>', $header, "\n", $left, ( '-' x $fillLength ), $right, "\n";
    }
}
close(DICE);

sub realignByAA($$) {
    my ( $ntS, $aaS )      = ( uc( $_[0] ), uc( $_[1] ) );
    my ( $r, $s, $codons ) = ( '', '', '' );
    my ( $gap, $i )        = ( 0, 0 );
    my @codons   = ();
    my @residues = ();

    $ntS =~ tr/-//d;
    my $ntL = length($ntS);
    my $gap = 3 - ( $ntL % 3 );
    if ( $gap != 3 ) {
        $ntS .= ( '-' x $gap );
        $ntL += $gap;
    }
    for ( $i = 0; $i < $ntL; $i += 3 ) {
        $codon = substr( $ntS, $i, 3 );
        if ( !defined( $gc{$codon} ) ) {
            $r = '?';
        } else {
            $r = $gc{$codon};
        }
        push( @codons,   $codon );
        push( @residues, $r );
    }

    $i = 0;
    for $r ( split( '', $aaS ) ) {
        if ( $r eq '-' ) {
            $s .= '---';
        } elsif ( $r eq $residues[$i] ) {
            $s .= $codons[$i];
            $i++;
        } else {
            $i++;
            while ( $i <= $#residues ) {
                if ( $r eq $residues[$i] ) {
                    $s .= $codons[$i];
                    $i++;
                    last;
                } else {
                    $i++;
                }
            }
        }
    }
    return lc($s);
}
