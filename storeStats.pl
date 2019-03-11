#!/usr/bin/env perl
# Samuel Shepard - 3.8.2019
# Stores alignment statistics for other uses.

use Storable qw(store_fd store);
use Getopt::Long;
GetOptions(
		'output|O=s' => \$outfile,
		'delim|D=s' => \$delim,
		'field|F=s' => \$fieldSet
	);

if ( -t STDIN && scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\tperl $0 <nts.fasta> [options]\n";
	$message .= "\t\t--output|-O <file.sto>\tOutput file for storable object. Default: STDOUT\n";
	$message .= "\t\t--delim|-D <CHAR>\tDelimiter for header fields.\n";
	$message .= "\t\t--field|-F <STR>\tComma-delimited set of fields to use for group. Default: no group.\n";
	die($message."\n");
}

$numberSelected = 0;
if ( defined($fieldSet) ) {
	@fields = split(',', $fieldSet);
	$numberSelected = scalar(@fields);
	foreach $x (@fields ) {
		if ( $x == 0 ) {
			die("$0 ERROR: field must be specified.\n");
		} elsif( $x < 0 ) {
			die("$0 ERROR: field must be a positive number.\n");
		}
	}
	for($x = 0; $x < $numberSelected; $x++ ) { $fields[$x]--; }
}

if ( !defined($delim) ) {
	$delim = '|';
} elsif( $delim eq '' ) {
	die("$0 ERROR: No delimiter argument detected.\n");
} elsif( length($delim) > 1 ) {
	die("$0 ERROR: single character delimiter expected instead of '$delim'.\n");
}

@patterns = (); $/ = ">";
while ( $record = <> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$id = shift(@lines);
	$sequence = lc(join('',@lines));
	$length = length($sequence);

	if ( $length == 0 ) {
		next;
	} elsif ( $length % 3 != 0 ) {
		die("$PROG:\tNot a codon alignment! See:\n\t$firstFile\n");
	} else {
		if ( $numberSelected > 0 ) {
			@values = split( /\Q$delim\E/, $id );
			$numberFound = scalar(@values);
			if ( $numberFound < $numberSelected ) {
				die("$0 ERROR: non-existant field specified. Wanted $numberSelected but found $numberFound\n");
			}
			$id = join($delim, (@values[@fields]) );
		}

		$group = $numberSelected > 0 ? $id : '__NIL__';
		$counts{$group}{$sequence}++;
	}
}


@groups = keys(%counts);
%codonData = ();
if ( scalar(@groups) == 1 && $groups[0] eq '__NIL__' ) {
	$group = '__NIL__';
	@patterns = keys(%{$counts{$group}});
	foreach my $pattern ( @patterns ) {
		$length = length($pattern);
		$count = $counts{$group}{$pattern};
		for ( $pos = 0; $pos < $length; $pos +=3 ) {
			# zero based codon number
			$codonNumber = int( $pos/3 );
			$codonStart = $codonNumber * 3;
			$codon = substr($pattern, $codonStart, 3);
			$codonData{$codonNumber}{$codon} += $count;
		}

	}
} else {
	foreach $group ( keys(%counts) ) {
		@patterns = keys(%{$counts{$group}});
		foreach my $pattern ( @patterns ) {
			$length = length($pattern);
			$count = $counts{$group}{$pattern};
			for ( $pos = 0; $pos < $length; $pos +=3 ) {
				# zero based codon number
				$codonNumber = int( $pos/3 );
				$codonStart = $codonNumber * 3;
				$codon = substr($pattern, $codonStart, 3);
				$codonData{$group}{$codonNumber}{$codon} += $count;
			}

		}
	}

}

#$g = 'CALI07|HA';
#foreach $p ( (0..3,545..548) ) {
#	print STDERR $p,"\n";
#	foreach $codon ( keys(%{$codonData{$g}{$p}}) ) {
#		print STDERR $g,"\t",$p,"\t",$codon,"\t",$codonData{$g}{$p}{$codon},"\n";
#	}
#}

if ( defined($outfile) ) {
	store(\%codonData, $outfile) or die("Cannot write to '$outfile'.\n");
} else {
	store_fd(\%codonData, *STDOUT) or die("Can't write to STDOUT.\n");
}
