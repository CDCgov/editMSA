#!/usr/bin/env perl
# Filename:         reviseTaxa.pl
# Description:      Allows taxa annotation manipulations for fasta files.
#
# Date dedicated:   2023-04-06
# Author:           Samuel S. Shepard, Centers for Disease Control and Prevention
#
# Citation:         Shepard SS, Davis CT, Bahl J, Rivailler P, York IA, Donis
#                   RO. LABEL: fast and accurate lineage assignment with
#                   assessment of H5N1 and H9N2 influenza A hemagglutinins. PLoS
#                   One. 2014;9(1):e86921. Published 2014 Jan 23.
#                   doi:10.1371/journal.pone.0086921
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

use warnings;
use strict;

use English qw(-no_match_vars);
use Getopt::Long;
use Carp qw(croak);
use File::Basename;

my ( $find, $replace, $joinAnnot, $addAnnot, $orderMode, $prefix, $suffix, $lastFieldDelim );
my $inPlace          = 0;
my $confirm          = 0;
my $deletePrev       = 0;
my $deleteSingle     = 0;
my $ignoreFastaAnnot = 0;
my $fuzzyMatch       = 0;
my $matchFiles       = 0;
my $prevInfix        = 0;
my $appendPipeAnnot  = 0;

Getopt::Long::Configure('no_ignore_case');
GetOptions(
            'find|F:s'             => \$find,
            'in-place|I'           => \$inPlace,
            'replace|R:s'          => \$replace,
            'join-to-end|J=s'      => \$joinAnnot,
            'confirm-prediction|C' => \$confirm,
            'delete-prev|D'        => \$deletePrev,
            'delete-single|S'      => \$deleteSingle,
            'add-annot|A=s'        => \$addAnnot,
            'ignore-fasta-annot|G' => \$ignoreFastaAnnot,
            'fuzzy-match|Z'        => \$fuzzyMatch,
            'order-mode|O:s'       => \$orderMode,
            'match-files|M'        => \$matchFiles,
            'previous-infix|N'     => \$prevInfix,
            'prefix|P=s'           => \$prefix,
            'append-pipe-annot|p'  => \$appendPipeAnnot,
            'suffix|X=s'           => \$suffix,
            'last-field|L=s'       => \$lastFieldDelim
);

if ( ( scalar(@ARGV) != 1 && !$matchFiles ) || ( $matchFiles && scalar(@ARGV) < 2 ) ) {
    die(  "Usage:\n\tperl $PROGRAM_NAME <input.fasta> [OPTIONS] [-M <file1 file2 ...>]\n"
        . "\t\t--confirm-prediction|-C\t\tConfirm predicted annotations.\n"
        . "\t\t--delete-prev|-D\t\tDelete previous annotations where there are two.\n"
        . "\t\t--delete-single|-S\t\tDelete previous annotations where there is one, filters AFTER delete-prev.\n"
        . "\t\t--find|-F <TEXT>\t\tSelect sequences including this annotation.\n"
        . "\t\t--replace|-R <TEXT>\t\tReplace the annotation with TEXT.\n"
        . "\t\t--add-annot|-A <FILE>\t\tAdd annotations based on tab delimited file (ID\tANNOT).\n"
        . "\t\t--ignore-fasta-annot|-G <FILE>\tIgnores previous annotation on FASTA headers vs. annotation file.\n"
        . "\t\t--fuzzy-match|-Z\t\tSearches for IDs in FASTA (-A option), matching if header contains the ID.\n"
        . "\t\t--order-mode|-O <out.file>\tAnnotation with ordinals for truncated names.\n"
        . "\t\t--match-files|-M <file1 ...>\tOutput filenames containing annotation names in the input file.\n"
        . "\t\t--previous-infix|N\t\tSuffix and prefix only apply to a 'previous' or first of a doublet annotation.\n"
        . "\t\t--prefix|-P <TEXT>\t\tPrefix for matching filenames with annotations OR adds prefix to header.\n"
        . "\t\t--suffix|-X <TEXT>\t\tSuffix for matching filenames with annotations OR adds suffix to header.\n"
        . "\t\t--in-place|-I\t\t\tRevise files in-place (could be dangerous).\n"
        . "\t\t--join-to-end|-J <TEXT>\t\tJoin annotation to end of the header (similar to add but without a file).\n"
        . "\t\t--last-field|-L <delim>\t\tClips the last field of the header and turns it into an annotation. Uses the specified delimiter to determine fields.\n"
        . "\t\t--append-pipe-annot|-p\t\tAppends annotation as a header field with pipe delim.\n" );
}

my $ORD;
if ( defined $orderMode ) {
    open( $ORD, '>', $orderMode ) or die("$PROGRAM_NAME ERROR: Cannot open $orderMode.\n");
}

my $lastField = defined $lastFieldDelim ? 1 : 0;
my %match     = ();
my %count     = ();
my %annotMap  = ();
my @annotIDs  = ();

if ($addAnnot) {
    open( my $ANNOT, '<', $addAnnot ) or die("$PROGRAM_NAME ERROR: Cannot open $addAnnot.\n");

    local $RS = "\n";
    while ( my $line = <$ANNOT> ) {
        chomp($line);
        my ( $key, $value ) = split( /\t/smx, $line );
        $annotMap{ uc($key) } = $value;
    }
    close $ANNOT or croak("Cannot close file: $OS_ERROR");
    @annotIDs = keys(%annotMap);
}

local $RS = ">";
open( my $IN, '<', $ARGV[0] ) or die("$PROGRAM_NAME ERROR: Cannot open $ARGV[0] for reading.\n");
my @records = <$IN>;
close $IN or croak("Cannot close file: $OS_ERROR");

my $OUT;
if ($inPlace) {
    open( $OUT, '>', $ARGV[0] ) or die("$PROGRAM_NAME ERROR: Cannot open $ARGV[0] for writing.\n");
}

# FNC - search for a header in the header hash database
sub headerInDB {
    my ( $ids, $keys, $header ) = @_;

    my $id = q{};

    if ( exists( $ids->{$header} ) ) {
        return ( 1, $header );
    }

    if ($fuzzyMatch) {
        foreach my $id ( @{$keys} ) {
            if ( $header =~ /\Q$id\E/smx ) {
                return ( 1, $id );
            }
        }
    }
    return ( 0, q{} );
}

# PROCESS stored fasta information
foreach my $fasta_record (@records) {
    chomp($fasta_record);
    my @lines    = split( /\r\n|\n|\r/smx, $fasta_record );
    my $id       = shift(@lines);
    my $sequence = join( q{}, @lines );
    my $annot    = q{};
    my @fields   = ();

    if ( length $sequence == 0 ) {
        next;
    }

    # Add annotations from a file
    if ($addAnnot) {
        my $tmp = $id;
        if ($ignoreFastaAnnot) {
            $tmp =~ s/_?\{.+?}//smx;
        }

        my ( $r, $newID ) = headerInDB( \%annotMap, \@annotIDs, uc($tmp) );
        if ($r) {
            if ($appendPipeAnnot) {
                $id = $id . '|' . $annotMap{ uc($newID) };
            } else {
                $id = $id . '{' . $annotMap{ uc($newID) } . '}';
            }
        }
    }

    # join to the end
    if ( defined $joinAnnot ) {
        $id .= '{' . $joinAnnot . '}';
    }

    # Get rid of a previous annotation.
    if ($deletePrev) {
        $id =~ s/{.+?}.*?\{(.+?)}/{$1}/smx;
    }

    if ($deleteSingle) {
        $id =~ s/_?\{.+?}//smx;
    }

    # Confirm prediction if they exist.
    if ($confirm) {
        $id =~ s/{PRED:(.+?)}$/{$1}/smx;
    }

    # Find and replace.
    if ( defined $find && defined $replace ) {
        $id =~ s/\Q{$find}\E/{$replace}/smx;

        # Always replace.
    } elsif ( defined $replace ) {
        if ( $id !~ s/{[^{}]*}\s*$/{$replace}/smx ) {
            $id .= "{$replace}";
        }
    }

    if ($lastField) {
        @fields = split( /\Q$lastFieldDelim\E/smx, $id );
        $annot  = pop(@fields);
        $id     = join( $lastFieldDelim, @fields ) . '{' . $annot . '}';
    }

    if ($orderMode) {
        if ( $id =~ /\{(.+?)\}/smx ) {
            $annot = $1;
        }

        if ( !defined( $count{$annot} ) ) {
            $count{$annot} = 1;
        } else {
            $count{$annot}++;
        }

        print $ORD $id, "\t", $count{$annot}, '_', $annot, "\n";
        $id = $count{$annot} . '_' . $annot;
    }

    # Check for prefix or suffix
    if ( ( defined $prefix || defined $suffix ) && !$matchFiles ) {
        if ($prevInfix) {
            $id =~ s/{(.+?)}\{/{$prefix$1$suffix}{/smx;
        } else {
            $id =~ s/{(.+?)}/{$prefix$1$suffix}/smx;
        }
    }

    if ($matchFiles) {
        if ( $id =~ /{(.+?)}/smx ) {
            $match{$1} = 1;
        }
    } elsif ($inPlace) {
        print $OUT '>', $id, "\n", $sequence, "\n";
    } else {
        print '>', $id, "\n", $sequence, "\n";
    }
}

if ($inPlace)             { close $OUT or croak("Cannot close file: $OS_ERROR"); }
if ( defined $orderMode ) { close $ORD or croak("Cannot close file: $OS_ERROR"); }

if ($matchFiles) {
    foreach my $i ( 1 .. $#ARGV ) {
        my $filename = basename( $ARGV[$i] );

        foreach my $annot ( keys(%match) ) {
            if ( $filename =~ /^$prefix$annot$suffix/ismx ) {
                print $ARGV[$i], q{ };
                last;
            }
        }
    }
}