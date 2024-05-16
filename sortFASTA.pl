#!/usr/bin/env perl
# sortFASTA - Version 1.1
# Sort the fasta file by fields.

# Samuel S. Shepard (vfn4@cdc.gov)
# 2012, Centers for Disease Control & Prevention

use Getopt::Long;
GetOptions(
            'ignore-case|C' => \$ignoreCase,
            'field|F=s'     => \$fieldSet,
            'delim|D=s'     => \$delim,
            'reverse|R'     => \$descending,
            'numerical|N'   => \$numerical
);

$numberSelected = 0;
if ( defined($fieldSet) ) {
    @fields         = split( ',', $fieldSet );
    $numberSelected = scalar(@fields);
    foreach $x (@fields) {
        if ( $x == 0 ) {
            die("$0 ERROR: field must be specified.\n");
        } elsif ( $x < 0 ) {
            die("$0 ERROR: field must be a positive number.\n");
        }
    }
    for ( $x = 0; $x < $numberSelected; $x++ ) { $fields[$x]--; }
}

if ( !defined($delim) ) {
    $delim = '|';
} elsif ( $delim eq '' ) {
    die("$0 ERROR: No delimiter argument detected.\n");
} elsif ( length($delim) > 1 ) {
    die("$0 ERROR: single character delimiter expected instead of '$delim'.\n");
}

if ( -t STDIN && !scalar(@ARGV) ) {
    $message = "Usage:\n\tperl $0 [FASTA ...] [OPTIONS]\n";
    $message .= "\t\t-C|--case-sensitive\tIgnore case.\n";
    $message .= "\t<STDIN> read if no FASTA given.\n";

    die($message);
}

# PROCESS fasta data
$/          = ">";
%keyByKeyID = ();
while ( $record = <> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/, $record );
    $header   = $h = shift(@lines);
    $sequence = join( '', @lines );

    if ( defined($fieldSet) ) {
        @values = split( /\Q$delim\E/, $h );
        $h      = $values[$fields[0]];
        foreach ( $i = 1; $i < $numberSelected; $i++ ) {
            if ( $fields[$i] > $#values ) {
                $fields[$i]++;
                die("$0 ERROR: non-existant field ($fields[$i]) specified.\n");
            }
            $k = $fields[$i];
            $h = $h . $delim . $values[$k];
        }
    }

    $length = length($sequence);
    if ( $length == 0 ) {
        next;
    } else {
        $seqByKeyID{$h}{$header} = $sequence;
    }
}
####################

@keys = ();
if ($numerical) {
    if ($descending) {
        @keys = sort { $b <=> $a } keys(%seqByKeyID);
    } else {
        @keys = sort { $a <=> $b } keys(%seqByKeyID);
    }
} else {
    if ($descending) {
        if ($ignoreCase) {
            @keys = sort { fc($b) cmp fc($a) } keys(%seqByKeyID);
        } else {
            @keys = sort { $b cmp $a } keys(%seqByKeyID);
        }
    } else {
        if ($ignoreCase) {
            @keys = sort { fc($a) cmp fc($b) } keys(%seqByKeyID);
        } else {
            @keys = sort { $a cmp $b } keys(%seqByKeyID);
        }
    }
}

foreach $key (@keys) {
    foreach $header ( keys( %{ $seqByKeyID{$key} } ) ) {
        print '>', $header, "\n", $seqByKeyID{$key}{$header}, "\n";
    }
}
