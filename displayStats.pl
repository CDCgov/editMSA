#!/usr/bin/env perl
# Displays the codon weight matrix (frequency stats).
#
# Samuel S. Shepard (vfn4@cdc.gov)
# 2022, Centers for Disease Control & Prevention

use strict;
use warnings;
use Storable;
use Getopt::Long;
use File::Basename;
use English qw( -no_match_vars);

my $verbose;
my $tomlFormat;
GetOptions( 'verbose|V' => \$verbose, 'toml-format|T' => \$tomlFormat );

if ( scalar(@ARGV) != 1 ) {
    die( "Usage:\n\tperl $PROGRAM_NAME <STATS.sto> [options]\n" . "\t\t--verbose|V\tMake display verbose.\n\n" );
}

# Obtain alignment-specific codon statistics if available
my $stats      = $ARGV[0];
my %codonStats = %{ retrieve($stats) } or die("Cannot open statistics file '$stats'.\n");

if ( defined $tomlFormat ) {
    my @keys = keys(%codonStats);
    print STDOUT 'module = "', basename( $stats, '.sto' ) . '"', "\n";
    if ( $keys[0] !~ /^\d+$/smx ) {
        print STDOUT "\n[[product]]\n";

        foreach my $g ( sort(@keys) ) {
            my ( $ref_id, $protein ) = split( '\|', $g );
            print STDOUT 'ref_id = "',              $ref_id,  "\"\n";
            print STDOUT 'protein = "',             $protein, "\"\n";
            print STDOUT 'codon_count_table = """', "\n";
            foreach my $p ( sort { $a <=> $b } keys( %{ $codonStats{$g} } ) ) {
                foreach
                  my $c ( sort { $codonStats{$g}{$p}{$b} <=> $codonStats{$g}{$p}{$a} } keys( %{ $codonStats{$g}{$p} } ) ) {
                    print STDOUT ( $p + 1 ), q{ }, $c, q{ }, $codonStats{$g}{$p}{$c}, "\n";
                }
            }
            print STDOUT '"""', "\n";
        }
    }

} elsif ( defined $verbose ) {
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
