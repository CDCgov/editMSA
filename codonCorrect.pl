#!/usr/bin/env perl
# Samuel Shepard - 4.18.2016
# Version 1.0
# Codon corrects sequences.

use File::Basename;
use Getopt::Long;
GetOptions( 'insertion-table|I=s' => \$insertionTable,
            'output-table|O=s'    => \$outputTable );

if ( -t STDIN && scalar(@ARGV) != 1 ) {
    $message = "Usage:\n\tperl $0 <nts.fasta> [options]\n";
    $message .= "\t\t--insertion-table|-I <STR>\tInsertion table for insertion corrections.\n";
    $message .= "\t\t--output-table|-O <STR>\t\tOutput file for the insertion table.\n";
    die( $message . "\n" );
}

if ($insertionTable) {
    $/       = "\n";
    %inserts = ();
    open( INS, '<', $insertionTable ) or die("Cannot open $insertionTable for reading.\n");
    @lines = <INS>;
    chomp(@lines);
    foreach $line (@lines) {
        ( $id, $pos, $insert ) = split( "\t", $line );
        $inserts{$id}{$pos} = lc($insert);
    }
    close(INS);
}

if ($outputTable) {
    open( TABL, '>', $outputTable ) or die("Cannot open $outputTable for writing.\n");
} else {
    *TABL = *STDERR;
}

$PROG      = basename( $0, '.pl' );
$firstFile = $ARGV[0];

# Process records.
%codonStats = %sequences = ();
@patterns   = ();
$/          = ">";
while ( $record = <> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/, $record );
    $id       = shift(@lines);
    $sequence = lc( join( '', @lines ) );
    $length   = length($sequence);

    if ( $length == 0 ) {
        next;
    } elsif ( $length % 3 != 0 ) {
        die("$PROG:\tNot a codon alignment! See:\n\t$firstFile\n");
    } else {
        $sequences{$id} = $sequence;
        $counts{$sequence}++;
    }
}

@patterns = keys(%counts);
foreach $id ( keys(%sequences) ) {
    $sequence = $sequences{$id};
    $seqLimit = length($sequence) - 2;    # TO-DO: what about the second to last insertion opportunity?
    if ( defined( $inserts{$id} ) ) {

        # position is 1 based
        foreach $pos ( keys( %{ $inserts{$id} } ) ) {
            if ( $pos > $seqLimit ) {

                # should be in frame if at length
                if ( $pos == length($sequence) ) { print TABL $id, "\t", $pos, "\t", uc($insert), "\n"; }
                next;
            }
            $insert = $inserts{$id}{$pos};
            $iDivis = length($insert) % 3;
            $iFrame = $pos % 3;
            if ( ( ( $iDivis == 0 ) && ( $iFrame != 0 ) ) ) {

                # zero based codon number
                $codonNumber = int( ( $pos - 1 ) / 3 );
                $codonStart  = $codonNumber * 3;
                if ( !defined( $codonStats{$codonNumber} ) ) {
                    foreach $pattern (@patterns) {
                        $codon = substr( $pattern, $codonStart, 3 );
                        $codonStats{$codonNumber}{$codon} += $counts{$pattern};
                    }
                }

                # A2 insertion
                if ( $iFrame == 2 ) {
                    $A2L2 = substr( $insert,   -2 ) . substr( $sequence, $codonStart + 2, 1 );
                    $A2R1 = substr( $sequence, $codonStart, 2 ) . substr( $insert, 0, 1 );

                    if ( $codonStats{$codonNumber}{$A2R1} >= $codonStats{$codonNumber}{$A2L2} ) {

                        #A2::R1 shift
                        $newInsert = substr( $insert, 1 ) . substr( $sequence, $pos, 1 );
                        $newCodon  = $A2R1;
                        $newPos    = $codonStart + 3;                                       # +2 +1 => +3 given 1-based
                    } else {

                        #A2::L2 shift
                        $newInsert = substr( $sequence, $codonStart, 2 ) . substr( $insert, 0, -2 );
                        $newCodon  = $A2L2;
                        $newPos    = $codonStart;                                           # -1 + 1 => + 0 given 1-based
                    }

                    # A1 insertion
                } else {
                    $A1L1 = substr( $insert,   -1 ) . substr( $sequence, $codonStart + 1, 2 );
                    $A1R2 = substr( $sequence, $codonStart, 1 ) . substr( $insert, 0, 2 );
                    if ( $codonStats{$codonNumber}{$A1L1} >= $codonStats{$codonNumber}{$A1R2} ) {

                        #A1::L1 shift
                        $newInsert = substr( $sequence, $codonStart, 1 ) . substr( $insert, 0, -1 );
                        $newCodon  = $A1L1;
                        $newPos    = $codonStart;                                           # -1 + 1 => + 0 given 1-based
                    } else {

                        #A1::R2 shift
                        $newInsert = substr( $insert, 2 ) . substr( $sequence, $pos, 2 );
                        $newCodon  = $A1R2;
                        $newPos    = $codonStart + 3;                                       # +2 +1 => +3 given 1-based
                    }
                }

                print TABL $id, "\t", $newPos, "\t", uc($newInsert), "\n";
                substr( $sequence, $codonStart, 3 ) = $newCodon;
            } else {
                print TABL $id, "\t", $pos, "\t", uc($insert), "\n";
            }
        }
    }

    # perform after insertion corrections
    while ( $sequence =~ /([A-Za-z]{3})((---)+)([A-Za-z]{3})/g ) {
        ( $left, $gaps, $gapT, $right ) = ( $1, $2, $3, $4 );
        $frame = $-[1] % 3;
        if ( $frame == 1 ) {

            # RIGHT SHIFT
            $replacement = substr( $left, 0, -1 ) . $gaps . substr( $left, -1 );
            substr( $sequence, $-[1], length($replacement) ) = $replacement;
        } elsif ( $frame == 2 ) {

            # LEFT SHIFT
            $replacement = $left . substr( $right, 0, 1 ) . $gaps;
            substr( $sequence, $-[1], length($replacement) ) = $replacement;
        }
    }

    print '>', $id, "\n", $sequence, "\n";
}

if ($outputTable) {
    close(TABL);
}
