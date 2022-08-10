#!/usr/bin/env perl
# Samuel Shepard - 4.18.2016
# Version 1.0
# Fill out

use Getopt::Long;
GetOptions( 'fill-triplets|F' => \$fillTriplets );

if ( scalar(@ARGV) != 2 ) {
    $message = "Usage:\n\tperl $0 <fasta> <ins.txt> [-F|--fill-triplets]\n";
    die( $message . "\n" );
}

%inserts        = ();
$/              = "\n";
$insertionTable = $ARGV[1];
open( INS, '<', $insertionTable ) or die("Cannot open $insertionTable for reading.\n");
@lines = <INS>;
chomp(@lines);
foreach $line (@lines) {
    ( $id, $pos, $insert ) = split( "\t", $line );
    $inserts{$id}{$pos} = $insert;
}
close(INS);

# Process records.
$/ = ">";
open( FASTA, '<', $ARGV[0] ) or die("Cannot open $ARGV[0] for reading.\n");
while ( $record = <FASTA> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/, $record );
    $id       = shift(@lines);
    $sequence = lc( join( '', @lines ) );
    $length   = length($sequence);

    if ( $length == 0 ) {
        next;
    }

    $offset = 0;
    foreach $pos ( sort { $a <=> $b } keys( %{ $inserts{$id} } ) ) {
        $insert = $inserts{$id}{$pos};

        if ( $fillTriplets && length($insert) % 3 != 0 ) {
            next;
        }

        substr( $sequence, $pos + $offset, 0 ) = $insert;
        $offset += length($insert);
    }
    print '>', $id, "\n", $sequence, "\n";
}
close(FASTA);
