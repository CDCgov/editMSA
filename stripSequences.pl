#!/usr/bin/env perl
# Filename:         stripSequences.pl
# Description:      Removes extraneous sequences from FASTA.
#
# Date dedicated:   2023-02-10
# Author:           Samuel S. Shepard, Centers for Disease Control and Prevention
#
# Citation:         Unpublished
#
# =============================================================================
#
#                            PUBLIC DOMAIN NOTICE
#
#  This source code file or script constitutes a work of the United States
#  Government and is not subject to domestic copyright protection under 17 USC §
#  105. This file is in the public domain within the United States, and
#  copyright and related rights in the work worldwide are waived through the CC0
#  1.0 Universal public domain dedication:
#  https://creativecommons.org/publicdomain/zero/1.0/
#
#  The material embodied in this software is provided to you "as-is" and without
#  warranty of any kind, express, implied or otherwise, including without
#  limitation, any warranty of fitness for a particular purpose. In no event
#  shall the Centers for Disease Control and Prevention (CDC) or the United
#  States (U.S.) government be liable to you or anyone else for any direct,
#  special, incidental, indirect or consequential damages of any kind, or any
#  damages whatsoever, including without limitation, loss of profit, loss of
#  use, savings or revenue, or the claims of third parties, whether or not CDC
#  or the U.S. government has been advised of the possibility of such loss,
#  however caused and on any theory of liability, arising out of or in
#  connection with the possession, use or performance of this software.
#
#  Please provide appropriate attribution in any work or product based on this
#  material.

use English qw(-no_match_vars);
use Getopt::Long;
use Carp qw(croak);
use warnings;
use strict;

my ( $fixHeader, $stripLower, $stripBadBases ) = ( 0, 0, 0 );
GetOptions( 'fix-header|F' => \$fixHeader, 'strip-lower|L' => \$stripLower, 'remove-bad-bases|N' => \$stripBadBases );

if (    ( ( $stripLower || $stripBadBases ) && scalar @ARGV != 1 )
     || ( ( !$stripLower && !$stripBadBases ) && scalar @ARGV != 2 )
     || ( $stripLower && $stripBadBases ) ) {

    die(   "\nUsage:\n\t$PROGRAM_NAME <file.fas> {-N|-L|<quoted_characters_to_delete>}\n"
         . "\t\t-N|--remove-bad-bases\tRemoves invalid nucleotide characters from the sequence.\n"
         . "\t\t-L|--strip-lower\tRemoves lowercase letters from the sequence.\n"
         . "\t\t-F|--fix-header\t\tRemoves and replaces troublesome characters from the FASTA header.\n"
         . "\n" );
}

open( my $IN, '<', $ARGV[0] ) or die("$PROGRAM_NAME ERROR: Cannot open $ARGV[0].\n");

# PREPARE the strip deletion
my $strip;
if ( scalar @ARGV > 1 ) {
    $strip = quotemeta( $ARGV[1] );                      # save input
    $strip = '$sequence =~ tr/' . $strip . '//d; 1;';    # create safe eval
}

local $RS = ">";
while ( my $fasta_record = <$IN> ) {
    chomp($fasta_record);
    my @lines    = split( /\r\n|\n|\r/smx, $fasta_record );
    my $header   = shift(@lines);
    my $sequence = join( q{}, @lines );
    if ( length($sequence) == 0 ) { next; }

    if ($fixHeader) {
        $header =~ s/^\s*(.*?)\s*$/$1/smx;
        $header =~ s/[\s:]/_/gsmx;
        $header =~ tr/',//d;
    }

    if ($stripLower) {
        $sequence =~ tr/[a-z]//d;
    } elsif ($stripBadBases) {
        $sequence =~ s/[^gcatrykmbvdhunGCATRYKMBVDHUN~.-]//gsmx;
    } else {
        eval($strip) or croak("Error in eval: $strip\n");

    }

    if ( length($sequence) == 0 ) {
        next;
    } else {
        print STDOUT '>', $header, "\n", $sequence, "\n";
    }
}
close $IN or croak("Could not close file: $OS_ERROR\n");
