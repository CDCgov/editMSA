#!/usr/bin/env perl
# removeGapColumns.pl - Version 1.0
# Removes gap columns (100% gaps) from FASTA alignments.
# Warning, this program overwrites the source file.
#
# Samuel S. Shepard (vfn4@cdc.gov)
# 2012, Centers for Disease Control & Prevention

use Getopt::Long;
GetOptions(
            'display-k-tons|D=i'    => \$ktons,
            'remove-k-tons|R'       => \$removeKtons,
            'tab-formatted|T'       => \$tabDelimited,
            'show-column-support|S' => \$showSupport,
            'a2m-insertions-only|A' => \$a2mInsertionsOnly
);

if ( scalar(@ARGV) != 1 ) {
    $message = "Usage:\n\tperl $0 <file.fasta> [options]\n";
    $message .= "\t\t-D|-display-k-tons\tDisplays sequences where for columns with 0 < bases <= K.\n";
    $message .= "\t\t-T|-tab-formatted\tPrint out tab-delimited data.\n";
    $message .= "\t\t-S|-show-column-support\tShow column support.\n";
    die($message);

}

open( IN, '<', $ARGV[0] ) or die("Cannot open $ARGV[0].\n");

# PROCESS fasta data
$/ = ">";
$i = 0;
my $first   = 1;
my @allowed = ();
my $del     = '-';
while ( $record = <IN> ) {
    chomp($record);
    @lines        = split( /\r\n|\n|\r/, $record );
    $header[$i]   = shift(@lines);
    $sequence[$i] = join( '', @lines );
    $length       = length( $sequence[$i] );

    if ( $length == 0 ) {
        next;
    }

    if ( defined($a2mInsertionsOnly) ) {
        if ($first) {
            my @C = split( '', $sequence[$i] );
            foreach my $j ( 0 .. $#C ) {
                if ( $C[$j] =~ /[a-z.]/ ) {
                    push( @allowed, $j );
                }
            }
            $first = 0;
            $del   = '.';
        }

        foreach my $j (@allowed) {
            if ( substr( $sequence[$i], $j, 1 ) eq '.' ) {
                $gaps{$j}++;
            }
        }
    } else {

        # unify gap characters
        $sequence[$i] =~ tr/:~./-/;

        $R = index( $sequence[$i], '-', 0 );
        while ( $R != -1 ) {
            $gaps{$R}++;
            $R = index( $sequence[$i], '-', $R + 1 );
        }
    }
    $i++;
}
close(IN);
$totalLength = length( $sequence[0] );
$numSeqs     = $i;

if ( $ktons > 0 ) {

    # Find and display sequences associated with k-ton support for each column.
    %ktons  = ();
    @sorted = sort { $a <=> $b } ( keys(%gaps) );
    foreach $g (@sorted) {
        $support = $numSeqs - $gaps{$g};
        if ( $support <= $ktons ) {
            print STDERR "$g:$support\n";
            $i = 0;
            while ( $support > 0 && $i < $numSeqs ) {
                $base = substr( $sequence[$i], $g, 1 );
                if ( $del ne $base ) {
                    print STDERR "\t", uc($base), " ", $header[$i], "\n";
                    if ($tabDelimited) {
                        print $g, "\t", $support, "\t", uc($base), "\t", $header[$i], "\n";
                    }
                    $support--;
                }
                $i++;
            }
        } else {
            delete( $gaps{$g} );
        }
    }

    if ( defined($removeKtons) ) {

        # Update file using an overwrite
        open( OUT, '>', $ARGV[0] ) or die("Cannot open $ARGV[0] for writing.\n");
        @sorted = sort { $b <=> $a } ( keys(%gaps) );
        for ( $i = 0; $i < $numSeqs; $i++ ) {
            foreach $R (@sorted) {
                substr( $sequence[$i], $R, 1, '' );
            }
            print OUT '>', $header[$i], "\n", $sequence[$i], "\n";
        }
        close(OUT);
    }
} else {

    # Find 100% gap columns
    foreach $g ( keys(%gaps) ) {
        if ( $gaps{$g} != $numSeqs ) {
            delete( $gaps{$g} );
        }
    }

    # Update file using an overwrite
    open( OUT, '>', $ARGV[0] ) or die("Cannot open $ARGV[0] for writing.\n");
    @sorted = sort { $b <=> $a } ( keys(%gaps) );
    for ( $i = 0; $i < $numSeqs; $i++ ) {
        foreach $R (@sorted) {
            substr( $sequence[$i], $R, 1, '' );
        }
        print OUT '>', $header[$i], "\n", $sequence[$i], "\n";
    }
    close(OUT);
}

$percentage = 0;
if ($showSupport) {
    for ( $i = 0; $i < $totalLength; $i++ ) {
        if ( defined( $gaps{$i} ) ) {
            $percentage = 1 - $gaps{$i} / $numSeqs;
            $percentage *= 100;
        } else {
            $percentage = 100;
        }
        print sprintf( "%d\t%5.1f\n", ( $i + 1 ), $percentage );
    }
}
