#!/usr/bin/env perl
# Sam Shepard - 2014

#<<< augmented translation table
my %gc = (
	'TAA'=>'*','TAG'=>'*','TAR'=>'*','TGA'=>'*','TRA'=>'*','GCA'=>'A','GCB'=>'A','GCC'=>'A','GCD'=>'A','GCG'=>'A','GCH'=>'A',
	'GCK'=>'A','GCM'=>'A','GCN'=>'A','GCR'=>'A','GCS'=>'A','GCT'=>'A','GCV'=>'A','GCW'=>'A','GCY'=>'A','TGC'=>'C','TGT'=>'C',
	'TGY'=>'C','GAC'=>'D','GAT'=>'D','GAY'=>'D','GAA'=>'E','GAG'=>'E','GAR'=>'E','TTC'=>'F','TTT'=>'F','TTY'=>'F','GGA'=>'G',
	'GGB'=>'G','GGC'=>'G','GGD'=>'G','GGG'=>'G','GGH'=>'G','GGK'=>'G','GGM'=>'G','GGN'=>'G','GGR'=>'G','GGS'=>'G','GGT'=>'G',
	'GGV'=>'G','GGW'=>'G','GGY'=>'G','CAC'=>'H','CAT'=>'H','CAY'=>'H','ATA'=>'I','ATC'=>'I','ATH'=>'I','ATM'=>'I','ATT'=>'I',
	'ATW'=>'I','ATY'=>'I','AAA'=>'K','AAG'=>'K','AAR'=>'K','CTA'=>'L','CTB'=>'L','CTC'=>'L','CTD'=>'L','CTG'=>'L','CTH'=>'L',
	'CTK'=>'L','CTM'=>'L','CTN'=>'L','CTR'=>'L','CTS'=>'L','CTT'=>'L','CTV'=>'L','CTW'=>'L','CTY'=>'L','TTA'=>'L','TTG'=>'L',
	'TTR'=>'L','YTA'=>'L','YTG'=>'L','YTR'=>'L','ATG'=>'M','AAC'=>'N','AAT'=>'N','AAY'=>'N','CCA'=>'P','CCB'=>'P','CCC'=>'P',
	'CCD'=>'P','CCG'=>'P','CCH'=>'P','CCK'=>'P','CCM'=>'P','CCN'=>'P','CCR'=>'P','CCS'=>'P','CCT'=>'P','CCV'=>'P','CCW'=>'P',
	'CCY'=>'P','CAA'=>'Q','CAG'=>'Q','CAR'=>'Q','AGA'=>'R','AGG'=>'R','AGR'=>'R','CGA'=>'R','CGB'=>'R','CGC'=>'R','CGD'=>'R',
	'CGG'=>'R','CGH'=>'R','CGK'=>'R','CGM'=>'R','CGN'=>'R','CGR'=>'R','CGS'=>'R','CGT'=>'R','CGV'=>'R','CGW'=>'R','CGY'=>'R',
	'MGA'=>'R','MGG'=>'R','MGR'=>'R','AGC'=>'S','AGT'=>'S','AGY'=>'S','TCA'=>'S','TCB'=>'S','TCC'=>'S','TCD'=>'S','TCG'=>'S',
	'TCH'=>'S','TCK'=>'S','TCM'=>'S','TCN'=>'S','TCR'=>'S','TCS'=>'S','TCT'=>'S','TCV'=>'S','TCW'=>'S','TCY'=>'S','ACA'=>'T',
	'ACB'=>'T','ACC'=>'T','ACD'=>'T','ACG'=>'T','ACH'=>'T','ACK'=>'T','ACM'=>'T','ACN'=>'T','ACR'=>'T','ACS'=>'T','ACT'=>'T',
	'ACV'=>'T','ACW'=>'T','ACY'=>'T','GTA'=>'V','GTB'=>'V','GTC'=>'V','GTD'=>'V','GTG'=>'V','GTH'=>'V','GTK'=>'V','GTM'=>'V',
	'GTN'=>'V','GTR'=>'V','GTS'=>'V','GTT'=>'V','GTV'=>'V','GTW'=>'V','GTY'=>'V','TGG'=>'W','TAC'=>'Y','TAT'=>'Y','TAY'=>'Y'
);
#>>>

use Getopt::Long;
GetOptions( 'begin|B=i' => \$begin, 'end|E=i' => \$end, 'shift-right|R' => \$right, 'translate|T' => \$translate );

if ( scalar(@ARGV) != 1 ) {
    $message = "Usage:\n\tperl $0 <in.fasta> <...>\n";
    $message .= "\t\t-B|--begin <INT>\tBegin coordinate.\n";
    $message .= "\t\t-E|--end <INT>\t\tEnd coordinate.\n";
    $message .= "\t\t-R|--shift-right\tShift alignment right.\n";
    die( $message . "\n" );
}

if ( !defined($begin) ) {
    $begin = 0;
} else {
    $begin--;
}

if ( !defined($end) ) {
    $end = 'L';
} else {
    $end--;
}

@fillSeqs = @fillLens = ();
$/        = ">";
$i        = 0;
while ( $record = <> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/, $record );
    $header   = shift(@lines);
    $sequence = uc( join( '', @lines ) );
    $length   = length($sequence);

    if ( $length == 0 ) {
        next;
    }

    #	print '>',$header,"\n";

    if ( $end eq 'L' ) {
        $stop = $length;
    } else {
        $stop = $end;
    }

    $subLen = $stop - $begin + 1;
    $subSeq = lc( substr( $sequence, $begin, $subLen ) );
    $subSeq =~ tr/-.~//d;
    $fill = $subLen - length($subSeq);

    if ($translateSort) {
        $fillSeqs[$i] = translate($subSeq);
    } else {
        $fillSeqs[$i] = $subSeq;
    }
    $fillLens[$i] = length($subSeq);

    if ($right) {
        $subSeq = ( '-' x $fill ) . $subSeq;
    } else {
        $subSeq = $subSeq . ( '-' x $fill );
    }
    substr( $sequence, $begin, $subLen ) = $subSeq;

    #	print $sequence,"\n";

    @records[$i] = '>' . $header . "\n" . $sequence, "\n";
    $i++;
}
$N     = $i;
@order = sort { $fillLens[$a] <=> $fillLens[$b] || $fillSeqs[$a] cmp $fillSeqs[$b] } 0 .. $#fillSeqs;

for ( $i = 0; $i < $N; $i++ ) {
    print $records[$order[$i]], "\n";
}

sub translate($) {
    my $nt = uc( $_[0] );
    my $aa = '';
    my $i  = 0;
    for ( $i = 0; $i < length($nt); $i += 3 ) {
        $aa .= $gc{ substr( $nt, $i, 3 ) };
    }
    return $aa;
}
