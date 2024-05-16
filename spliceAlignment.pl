#!/usr/bin/env perl
# Sam Shepard - 2016

#<<< augmented translation table
my %gc = (
	'TAA'=>'*','TAG'=>'*','TAR'=>'*','TGA'=>'*','TRA'=>'*','GCA'=>'A','GCB'=>'A','GCC'=>'A','GCD'=>'A','GCG'=>'A','GCH'=>'A',
	'GCK'=>'A','GCM'=>'A','GCN'=>'A','GCR'=>'A','GCS'=>'A','GCT'=>'A','GCV'=>'A','GCW'=>'A','GCY'=>'A','TGC'=>'C','TGT'=>'C',
	'TGY'=>'C','GAC'=>'D','GAT'=>'D','GAY'=>'D','GAA'=>'E','GAG'=>'E','GAR'=>'E','TTC'=>'F','TTT'=>'F','TTY'=>'F','GGA'=>'G',
	'GGB'=>'G','GGC'=>'G','GGD'=>'G','GGG'=>'G','GGH'=>'G','GGK'=>'G','GGM'=>'G','GGN'=>'G','GGR'=>'G','GGS'=>'G','GGT'=>'G',
	'GGV'=>'G','GGW'=>'G','GGY'=>'G','CAC'=>'H','CAT'=>'H','CAY'=>'H','ATA'=>'I','ATC'=>'I','ATH'=>'I','ATM'=>'I','ATT'=>'I',
	'ATW'=>'I','ATY'=>'I','AAA'=>'K','AAG'=>'K','AAR'=>'K','CTA'=>'L','CTB'=>'L','CTC'=>'L','CTD'=>'L','CTG'=>'L','CTH'=>'L',
	'CTK'=>'L','CTM'=>'L','CTN'=>'L','CTR'=>'L','CTS'=>'L','CTT'=>'L','CTV'=>'L','CTW'=>'L','CTY'=>'L','TTA'=>'L','TTG'=>'L',
	'TTR'=>'L','YTA'=>'L','YTG'=>'L','YTR'=>'L','ATG'=>'M','AAC'=>'N','AAT'=>'N','AAY'=>'N','CCA'=>'P','CCB'=>'P','CCC'=>'P',
	'CCD'=>'P','CCG'=>'P','CCH'=>'P','CCK'=>'P','CCM'=>'P','CCN'=>'P','CCR'=>'P','CCS'=>'P','CCT'=>'P','CCV'=>'P','CCW'=>'P',
	'CCY'=>'P','CAA'=>'Q','CAG'=>'Q','CAR'=>'Q','AGA'=>'R','AGG'=>'R','AGR'=>'R','CGA'=>'R','CGB'=>'R','CGC'=>'R','CGD'=>'R',
	'CGG'=>'R','CGH'=>'R','CGK'=>'R','CGM'=>'R','CGN'=>'R','CGR'=>'R','CGS'=>'R','CGT'=>'R','CGV'=>'R','CGW'=>'R','CGY'=>'R',
	'MGA'=>'R','MGG'=>'R','MGR'=>'R','AGC'=>'S','AGT'=>'S','AGY'=>'S','TCA'=>'S','TCB'=>'S','TCC'=>'S','TCD'=>'S','TCG'=>'S',
	'TCH'=>'S','TCK'=>'S','TCM'=>'S','TCN'=>'S','TCR'=>'S','TCS'=>'S','TCT'=>'S','TCV'=>'S','TCW'=>'S','TCY'=>'S','ACA'=>'T',
	'ACB'=>'T','ACC'=>'T','ACD'=>'T','ACG'=>'T','ACH'=>'T','ACK'=>'T','ACM'=>'T','ACN'=>'T','ACR'=>'T','ACS'=>'T','ACT'=>'T',
	'ACV'=>'T','ACW'=>'T','ACY'=>'T','GTA'=>'V','GTB'=>'V','GTC'=>'V','GTD'=>'V','GTG'=>'V','GTH'=>'V','GTK'=>'V','GTM'=>'V',
	'GTN'=>'V','GTR'=>'V','GTS'=>'V','GTT'=>'V','GTV'=>'V','GTW'=>'V','GTY'=>'V','TGG'=>'W','TAC'=>'Y','TAT'=>'Y','TAY'=>'Y'
);
#>>>

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
