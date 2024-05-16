#!/usr/bin/env perl
# Filename:         storeStats.pl
# Description:      Computes & stores the codon weight matrix from an MPSA.
# Author:           Samuel S. Shepard, Centers for Disease Control and Prevention

use strict;
use warnings;
use English qw( -no_match_vars);
use Storable qw(store_fd store);
use Getopt::Long;

my ( $outfile, $delim, $fieldSet );
GetOptions(
            'output|O=s' => \$outfile,
            'delim|D=s'  => \$delim,
            'field|F=s'  => \$fieldSet
);

if ( -t STDIN && scalar(@ARGV) != 1 ) {
    die(   "Usage:\n\tperl $PROGRAM_NAME <nts.fasta> [options]\n"
         . "\t\t--output|-O <file.sto>\tOutput file for storable object. Default: STDOUT\n"
         . "\t\t--delim|-D <CHAR>\tDelimiter for header fields. Default delim is '|'.\n"
         . "\t\t--field|-F <STR>\tComma-delimited set of fields to use for group. Default: no group.\n"
         . "\n" );
}

my $numberSelected = 0;
my @fields         = ();
if ( defined $fieldSet ) {
    @fields         = split( ',', $fieldSet );
    $numberSelected = scalar(@fields);
    foreach my $x (@fields) {
        if ( $x == 0 ) {
            die("$PROGRAM_NAME ERROR: field must be specified.\n");
        } elsif ( $x < 0 ) {
            die("$PROGRAM_NAME ERROR: field must be a positive number.\n");
        }
    }
    foreach my $x ( 0 .. ( $numberSelected - 1 ) ) {
        $fields[$x]--;
    }
}

if ( !defined $delim ) {
    $delim = '|';
} elsif ( $delim eq q{} ) {
    die("$PROGRAM_NAME ERROR: No delimiter argument detected.\n");
} elsif ( length($delim) > 1 ) {
    die("$PROGRAM_NAME ERROR: single character delimiter expected instead of '$delim'.\n");
}

my %counts = ();
local $RS = ">";
while ( my $fasta_record = <> ) {
    chomp($fasta_record);
    my @lines    = split( /\r\n|\n|\r/smx, $fasta_record );
    my $id       = shift(@lines);
    my $sequence = lc( join( q{}, @lines ) );
    my $length   = length($sequence);

    if ( $length == 0 ) {
        next;
    } elsif ( $length % 3 != 0 ) {
        die("$PROGRAM_NAME:\tNot a codon alignment, length $length not divisible by 3!\n");
    } else {
        if ( $numberSelected > 0 ) {
            my @values      = split( /\Q$delim\E/smx, $id );
            my $numberFound = scalar(@values);
            if ( $numberFound < $numberSelected ) {
                die("$PROGRAM_NAME ERROR: non-existant field specified. Wanted $numberSelected but found $numberFound\n");
            }
            $id = join( $delim, ( @values[@fields] ) );
        }

        my $group = $numberSelected > 0 ? $id : '__NIL__';
        $counts{$group}{$sequence}++;
    }
}

my @groups    = keys(%counts);
my %codonData = ();
if ( scalar @groups == 1 && $groups[0] eq '__NIL__' ) {
    my $group    = '__NIL__';
    my @patterns = keys( %{ $counts{$group} } );
    foreach my $pattern (@patterns) {
        my $length = length($pattern);
        my $count  = $counts{$group}{$pattern};

        ## no critic (ControlStructures::ProhibitCStyleForLoops)
        for ( my $pos = 0; $pos < $length; $pos += 3 ) {

            # zero based codon number
            my $codonNumber = int( $pos / 3 );
            my $codonStart  = $codonNumber * 3;
            my $codon       = substr( $pattern, $codonStart, 3 );
            if ( $codon =~ /[.nN-]/smx ) { next; }
            $codonData{$codonNumber}{$codon} += $count;
        }

    }
} else {
    foreach my $group ( keys(%counts) ) {
        my @patterns = keys( %{ $counts{$group} } );
        foreach my $pattern (@patterns) {
            my $length = length($pattern);
            my $count  = $counts{$group}{$pattern};

            ## no critic (ControlStructures::ProhibitCStyleForLoops)
            for ( my $pos = 0; $pos < $length; $pos += 3 ) {

                # zero based codon number
                my $codonNumber = int( $pos / 3 );
                my $codonStart  = $codonNumber * 3;
                my $codon       = substr( $pattern, $codonStart, 3 );
                if ( $codon =~ /[.nN-]/smx ) { next; }
                $codonData{$group}{$codonNumber}{$codon} += $count;
            }

        }
    }

}

if ( defined($outfile) ) {
    store( \%codonData, $outfile ) or die("Cannot write to '$outfile'.\n");
} else {
    store_fd( \%codonData, *STDOUT ) or die("Can't write to STDOUT.\n");
}
