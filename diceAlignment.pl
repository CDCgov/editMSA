#!/usr/bin/env perl
# Dices alignments into easier to manage smaller regions.
#
# Samuel S. Shepard (vfn4@cdc.gov)
# 2016, Centers for Disease Control & Prevention

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

use File::Basename;
use Getopt::Long;
GetOptions(
            'begin|B=i'       => \$begin,
            'end|E=i'         => \$end,
            'shift-left|L'    => \$shiftLeft,
            'fix-coords|F'    => \$fixCoords,
            'shift-right|R'   => \$shiftRight,
            'translate|T'     => \$translate,
            'prefix|P=s'      => \$prefix,
            'sort-by-count|C' => \$sortByCount
);

if ( scalar(@ARGV) != 1 ) {
    $message = "Usage:\n\tperl $0 <in.fasta> <...>\n";
    $message .= "\t\t-P|--prefix <STR>\tPrefix used for output.\n";
    $message .= "\t\t-B|--begin <INT>\tBegin coordinate.\n";
    $message .= "\t\t-E|--end <INT>\t\tEnd coordinate.\n";
    $message .= "\t\t-R|--shift-right\tShift alignment right.\n";
    $message .= "\t\t-T|--translate\t\tTranslate data as well.\n";
    $message .= "\t\t-C|--sort-by-count\tSort data by their cluster counts.\n";
    die( $message . "\n" );
}

if ( !defined($begin) ) {
    $begin = 0;
} else {
    $begin--;
}

if ($fixCoords) {
    if ( ( $begin % 3 ) != 0 ) {
        $tmp = $begin - ( $begin % 3 );
        if ( $tmp >= 0 ) {
            print STDERR "Moving starting coordinate up to ", ( $tmp + 1 ), ".\n";
            $begin = $tmp;
        }
    }
}

if ( !defined($end) ) {
    $end = 'L';
} else {
    $end--;
}

if ( !defined($prefix) ) {
    $filename = basename( $ARGV[0] );
    @pieces   = split( '\.', $filename );
    if ( scalar(@pieces) > 1 ) {
        pop(@pieces);
    }
    $prefix = join( '.', @pieces );
}

open( DICED, '>', $prefix . '.diced.txt' )    or die("Cannot write to $prefix.diced.txt\n");
open( MINI,  '>', $prefix . '.miniNT.fasta' ) or die("Cannot write to $prefix.mini.fasta\n");
%clusterCount = %cluster = ();
$/            = ">";
$cluster      = 0;
while ( $record = <> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/, $record );
    $header   = shift(@lines);
    $sequence = join( '', @lines );
    $length   = length($sequence);

    if ( $length == 0 ) { next; }
    if ( $end eq 'L' ) {
        $stop = $length;
    } else {
        $stop = $end;
    }
    $subLen = $stop - $begin + 1;
    if ($fixCoords) {
        $fixCoords = 0;
        if ( ( $subLen % 3 ) != 0 ) {
            $tmp = $stop - ( $subLen % 3 );
            if ( $stop > $begin ) {
                print STDERR "Moving stopping coordinate up to ", ( $tmp + 1 ), ".\n";
                $stop   = $end = $tmp;
                $subLen = $stop - $begin + 1;
            }
        }
    }
    $pattern = $mid = lc( substr( $sequence, $begin, $subLen ) );
    $pattern =~ tr/-.~//d;

    $left = $right = "";
    if ( $begin > 0 ) {
        $left = uc( substr( $sequence, 0, $begin ) );
    }

    if ( $stop < $length ) {
        $right = uc( substr( $sequence, $stop + 1, ( $length - $stop - 1 ) ) );
    }

    if ( !defined( $clusters{$pattern} ) ) {
        $cluster++;
        $clusters{$pattern}     = $cluster;
        $clusterCount{$cluster} = 1;
        $originals{$pattern}    = $mid;
    } else {
        $clusterCount{ $clusters{$pattern} }++;
    }
    print DICED $header, "\t", $left, "\t", $clusters{$pattern}, "\t", $right, "\n";
}

@fillLens = @fillPats = ();
@patterns = keys(%originals);
if ($shiftRight) {
    foreach $i ( 0 .. $#patterns ) {
        $pattern      = $patterns[$i];
        $fillLens[$i] = length($pattern);
        $fill         = $subLen - $fillLens[$i];
        $fillPats[$i] = ( '-' x $fill ) . $pattern;
    }
} elsif ($shiftLeft) {
    foreach $i ( 0 .. $#patterns ) {
        $pattern      = $patterns[$i];
        $fillLens[$i] = length($pattern);
        $fill         = $subLen - $fillLens[$i];
        $fillPats[$i] = $pattern . ( '-' x $fill );
    }
} else {
    foreach $i ( 0 .. $#patterns ) {
        $fillLens[$i] = length( $patterns[$i] );
        $fillPats[$i] = $originals{ $patterns[$i] };
    }
}

if ($sortByCount) {
    @order = sort {
             $clusterCount{ $clusters{ $patterns[$b] } } <=> $clusterCount{ $clusters{ $patterns[$a] } }
          || $fillPats[$a] cmp $fillPats[$b]
    } 0 .. $#patterns;
} else {
    @order = sort { $fillLens[$a] <=> $fillLens[$b] || $fillPats[$a] cmp $fillPats[$b] } 0 .. $#patterns;
}

for $i ( 0 .. $#order ) {
    $index   = $order[$i];
    $pattern = $patterns[$index];
    $cluster = $clusters{$pattern};
    print MINI '>', $cluster, '|', $clusterCount{$cluster}, "\n";
    print MINI $fillPats[$index], "\n";
}
close(MINI);
close(DICED);
print STDERR "Created '$prefix.miniNT.fasta' for editing.\n";

if ($translate) {
    open( AA, '>', $prefix . '.miniAA.fasta' ) or die("Cannot write to $prefix.miniAA.fasta\n");
    $aaLength = 0;
    %aaFill   = %ntClusterList = ();
    for $pattern (@patterns) {
        $aa = translate($pattern);
        if ( length($aa) > $aaLength ) {
            $aaLength = length($aa);
        }

        if ( !defined( $aaClusterCount{$aa} ) ) {
            $aaClusterCount{$aa} = $clusterCount{ $clusters{$pattern} };
        } else {
            $aaClusterCount{$aa} += $clusterCount{ $clusters{$pattern} };
        }
        push( @{ $ntClusterList{$aa} }, $clusters{$pattern} );
    }

    for $aa ( keys(%aaClusterCount) ) {
        $gap = $aaLength - length($aa);
        if ($shiftLeft) {
            $fill = $aa . ( '-' x $gap );
        } elsif ($shiftRight) {
            $fill = ( '-' x $gap ) . $aa;
        } else {
            if ( $gap == 0 ) {
                $fill = $aa;
            } elsif ( $gap % 2 == 0 ) {
                $gap /= 2;
                $fill = ( '-' x $gap ) . $aa . ( '-' x $gap );
            } else {
                $gap  = int( $gap / 2 );
                $fill = ( '-' x $gap ) . $aa . ( '-' x ( $gap + 1 ) );
            }
        }
        $aaFill{$fill} = $aa;
    }

    if ($sortByCount) {
        @peptides = sort { $aaClusterCount{ $aaFill{$b} } <=> $aaClusterCount{ $aaFill{$a} } } keys(%aaFill);
    } else {
        @peptides = sort { $a cmp $b } keys(%aaFill);
    }
    for $fill (@peptides) {
        $pat = $aaFill{$fill};
        print AA '>', join( '-', @{ $ntClusterList{$pat} } ), '|', $aaClusterCount{$pat}, "\n", $fill, "\n";
    }
    close(AA);
    print STDERR "Created '$prefix.miniAA.fasta' for editing.\n";
}

sub translate($) {
    my $nt    = uc( $_[0] );
    my $aa    = '';
    my $i     = 0;
    my $codon = '';
    for ( $i = 0; $i < length($nt); $i += 3 ) {
        $codon = substr( $nt, $i, 3 );
        if ( !defined( $gc{$codon} ) ) {
            $aa .= '?';
        } else {
            $aa .= $gc{$codon};
        }
    }
    return $aa;
}
