#!/usr/bin/env perl
# Filename:         displayStats
# Description:      Displays the codon weight matrix (frequency stats).
#
# Date dedicated:   2022-08-10
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

use strict;
use warnings;
use Storable;
use Getopt::Long;
use English qw( -no_match_vars);

my $verbose;
GetOptions( 'verbose|V' => \$verbose );

if ( scalar(@ARGV) != 1 ) {
    die( "Usage:\n\tperl $PROGRAM_NAME <STATS.sto> [options]\n" . "\t\t--verbose|V\tMake display verbose.\n\n" );
}

# Obtain alignment-specific codon statistics if available
my $stats      = $ARGV[0];
my %codonStats = %{ retrieve($stats) } or die("Cannot open statistics file '$stats'.\n");

if ( defined($verbose) ) {
    my @keys = keys(%codonStats);
    if ( $keys[0] =~ /^\d+$/smx ) {
        foreach my $p ( sort { $a <=> $b } @keys ) {
            foreach my $c ( sort { $codonStats{$p}{$b} <=> $codonStats{$p}{$a} } keys( %{ $codonStats{$p} } ) ) {
                print STDOUT $p, "\t", $c, "\t", $codonStats{$p}{$c}, "\n";
            }
        }
    } else {
        foreach my $g ( sort(@keys) ) {
            print STDOUT $g, "\n";
            foreach my $p ( sort { $a <=> $b } keys( %{ $codonStats{$g} } ) ) {
                foreach
                  my $c ( sort { $codonStats{$g}{$p}{$b} <=> $codonStats{$g}{$p}{$a} } keys( %{ $codonStats{$g}{$p} } ) ) {
                    print STDOUT $p, "\t", $c, "\t", $codonStats{$g}{$p}{$c}, "\n";
                }
            }
            print STDOUT "\n";
        }
    }
} else {
    my $maxC = 0;
    my @keys = keys(%codonStats);
    if ( $keys[0] =~ /^\d+$/smx ) {
        @keys = sort { $a <=> $b } @keys;
        print STDOUT "Number of Codons\tMaximum Codon Depth\n";
        foreach my $p (@keys) {
            foreach my $c ( sort { $codonStats{$p}{$b} <=> $codonStats{$p}{$a} } keys( %{ $codonStats{$p} } ) ) {
                if ( $codonStats{$p}{$c} > $maxC ) {
                    $maxC = $codonStats{$p}{$c};
                }
            }
        }
        print STDOUT $#keys, "\t", $maxC, "\n";
    } else {
        print STDOUT "Group\tNumber of Codons\tMaximum Codon Depth\n";
        foreach my $g ( sort(@keys) ) {
            my @positions = sort { $a <=> $b } keys( %{ $codonStats{$g} } );
            $maxC = 0;

            foreach my $p (@positions) {
                foreach
                  my $c ( sort { $codonStats{$g}{$p}{$b} <=> $codonStats{$g}{$p}{$a} } keys( %{ $codonStats{$g}{$p} } ) ) {
                    if ( $codonStats{$g}{$p}{$c} > $maxC ) {
                        $maxC = $codonStats{$g}{$p}{$c};
                    }
                }
            }
            print STDOUT $g, "\t", $#positions, "\t", $maxC, "\n";
        }
    }
}
