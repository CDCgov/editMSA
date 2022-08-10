#!/usr/bin/env perl
# Sam Shepard - treeOrderedFasta.pl - 5.2011
# Gets a fasta file and orders by tree order.
#
# Version 1.0.3
# 1.0.3 + Improved capatibility with MEGA4's header mutilations.
# 1.0.2 + Added a quoted ID option for MEGA5.
# 1.0.2 * Possible warning mix-up.
# 1.0.1 + Gives warnings for fasta headers not found.
# 1.0.1 * Only outputs fasta headers with sequences.
# 1.0.1 * Fixed a problem with apostrophes in headers.

use Getopt::Long;
GetOptions( 'quoted-ids|Q' => \$quoted );

if ( scalar(@ARGV) != 3 ) {
    die("Usage:\n\t$0 <tree.nwk> <sequences.fasta> <output.fasta> [-Q|--quoted-ids]\n");
}

open( IN,  '<', $ARGV[0] ) or die("Cannot open $ARGV[0].\n");
open( SEQ, '<', $ARGV[1] ) or die("Cannot open $ARGV[1].\n");
open( OUT, '>', $ARGV[2] ) or die("Cannot open $ARGV[2].\n");
$debug = 0;

# Process Newick file.
@ordered = ();
%tree    = ();

if ($quoted) {
    $/    = "\n";
    $line = <IN>;
    chomp($line);

    #if ( $line =~ s/([^'\d])(A\/.+?)([,:)])/\1'\2'\3/g ) {
    #	print STDERR "TREEORDER: Warning, sequence adjustment was required!\n";
    #}

    while ( $line =~ m/'(.+?)':/g ) {
        $ordered[$i] = $1;
        $tree{$1} = 'T';
        $i++;
    }
} else {
    $/     = ':';
    $block = <IN>;
    chomp($block);
    $block =~ tr/'//d;
    $block =~ /\(*(.+)$/;
    $tree{$1}   = 'T';
    $ordered[0] = $1;
    $i          = 1;

    while ( $block = <IN> ) {
        chomp($block);
        $block =~ tr/'//d;
        @parts = split( ',', $block );

        if ( scalar(@parts) > 1 ) {
            $parts[1] =~ /\(*(.+)$/;
            $ordered[$i] = $1;
            $tree{$1} = 'T';
            $i++;
        }
    }
}
close(IN);

if ($debug) {
    open( DEBUG, '>', 'debugTree.log' );
    foreach $item ( sort(@ordered) ) {
        print DEBUG $item, "\n";
    }
    close(DEBUG);
}

# Process Fasta file.
$/         = ">";
%sequences = ();
while ( $record = <SEQ> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/, $record );
    $header   = shift(@lines);
    $sequence = join( '', @lines );

    if ( length($sequence) == 0 ) {
        next;
    }

    if ($quoted) {
        $header =~ tr/'//d;
        $header =~ tr/ /_/;
    }

    if ( $tree{$header} eq 'T' ) {
        $sequences{$header} = $sequence;
    } else {

        #if( $header =~ /(.*?H\d{1,2}N\d(_NP){0,1}).*$/ ) {
        #	$key = $1;
        #} else {
        #	$key = $header;
        #}
        if ( $header =~ /(.*?)_\{.+\}/ ) {
            $header = $1;
        }

        $key = substr( $header, 0, 40 );
        $key =~ tr/:)(//d;
        if ( substr( $key, -1 ) eq '_' ) {
            chop($key);
        }
        $key =~ s/__/_/g;

        $sequences{$key} = $sequence;
    }
}
close(SEQ);

if ($debug) {
    open( DEBUG, '>', 'debugFasta.log' );
    foreach $item ( sort( keys(%sequences) ) ) {
        print DEBUG $item, "\n";
    }
    close(DEBUG);
}

# Process output
$notFound = 0;
foreach $key (@ordered) {
    if ( !defined( $sequences{$key} ) ) {
        $notFound++;
        print("WARNING! Cannot find header in sequence file: $key\n");
    } else {
        print OUT '>', $key, "\n";
        print OUT $sequences{$key}, "\n";
    }
}
close(OUT);

if ( $notFound > 1 ) {
    print "There were $notFound sequences not found in the fasta sequence file.\n";
}
