#!/usr/bin/env perl
# Samuel Shepard - 2019.03.11
# Version 1.0
# Codon corrects sequences from user-supplied statistics.

use warnings;
use Storable;
use File::Basename;
use Getopt::Long;
use constant { NIL => '__NIL__' };
GetOptions(
            'stats|S=s'           => \$statsFile,
            'insertion-table|I=s' => \$insertionTable,
            'output-table|O=s'    => \$outputTable,
            'delim|D=s'           => \$delim,
            'field|F=s'           => \$fieldSet
);

if ( -t STDIN && scalar(@ARGV) != 1 ) {
    $message = "Usage:\n\tperl $0 <nts.fasta> [options]\n";
    $message .= "\t\t--stats|-S <in.sto>\t\tInput file for storable object. Default: none\n";
    $message .= "\t\t--delim|-D <CHAR>\t\tDelimiter for header fields.\n";
    $message .= "\t\t--field|-F <STR>\t\tComma-delimited set of fields to use for group. Default: no group.\n";
    $message .= "\t\t--insertion-table|-I <STR>\tInsertion table for insertion corrections.\n";
    $message .= "\t\t--output-table|-O <STR>\t\tOutput file for the insertion table.\n";
    die( $message . "\n" );
}

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

my %inserts = ();
if ($insertionTable) {
    $/ = "\n";
    open( INS, '<', $insertionTable ) or die("Cannot open $insertionTable for reading.\n");
    @lines = <INS>;
    chomp(@lines);
    foreach $line (@lines) {
        ( $id, $pos, $insert ) = split( "\t", $line );
        $inserts{$id}{$pos} = lc($insert);
    }
    close(INS);
}

if ($outputTable) {
    open( TABL, '>', $outputTable ) or die("Cannot open $outputTable for writing.\n");
} else {
    *TABL = *STDERR;
}

# From: http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=465364&aa=1&style=GCG
# Stop codons are changed to pseudo-counts to avoid pre-mature stop codons.
# Old data:	'tag'=>  '8905','taa'=> '19221','tga'=> '12569'
#<<< Influenza A codon table
my %defaultCodonStats = (
    'ggg'=>'298542','gga'=>'491304','ggt'=>'162838','ggc'=>'148927','gag'=>'513739','gaa'=>'667159','gat'=>'429380','gac'=>'332513',
    'gtg'=>'323574','gta'=>'205254','gtt'=>'220590','gtc'=>'183609','gcg'=> '74308','gca'=>'398856','gct'=>'254219','gcc'=>'207403',
    'agg'=>'299901','aga'=>'521693','agt'=>'237994','agc'=>'252010','aag'=>'365188','aaa'=>'591384','aat'=>'515196','aac'=>'387110',
    'atg'=>'622324','ata'=>'402080','att'=>'389594','atc'=>'298586','acg'=> '80019','aca'=>'455971','act'=>'301672','acc'=>'217249',
    'tgg'=>'270035','tgt'=>'125361','tgc'=>'184026','tat'=>'233240','tac'=>'207199','ttg'=>'240030','tta'=>'133967','ttt'=>'293396',
    'ttc'=>'332322','tcg'=> '61962','tca'=>'312168','tct'=>'210832','tcc'=>'181613','cgg'=> '91717','cga'=> '99775','cgt'=> '33533',
    'cgc'=> '52445','cag'=>'302804','caa'=>'371086','cat'=>'162471','cac'=>'118729','ctg'=>'243631','cta'=>'219434','ctt'=>'276221',
    'ctc'=>'200260','ccg'=> '76633','cca'=>'242107','cct'=>'185417','ccc'=>'119819','tag'=>     '1','taa'=>     '3','tga'=>     '2'
);
#>>>

# Obtain alignment-specific codon statistics if available
my %codonStats = ();
if ( defined($statsFile) ) {
    %codonStats = %{ retrieve($statsFile) } or die("Cannot open statistics file '$statsFile'.\n");
    $stats      = 1;
} else {
    $stats = 0;
}

# is LEFT >= RIGHT
sub compareLeftRightGE($$$$) {
    my ( $left, $right, $group, $pos ) = ( @_[0 .. 3] );
    my ( $x, $y ) = ( 0, 0 );
    if ( !defined($group) || $group eq NIL || $group eq '' ) {
        $x = defined( $codonStats{$pos}{$left} )  ? $codonStats{$pos}{$left}  : 0;
        $y = defined( $codonStats{$pos}{$right} ) ? $codonStats{$pos}{$right} : 0;
    } else {
        $x = defined( $codonStats{$group}{$pos}{$left} )  ? $codonStats{$group}{$pos}{$left}  : 0;
        $y = defined( $codonStats{$group}{$pos}{$right} ) ? $codonStats{$group}{$pos}{$right} : 0;
    }

    if ( $x == $y && $x == 0 ) {
        $x = defined( $defaultCodonStats{$pos}{$left} )  ? $defaultCodonStats{$pos}{$left}  : 0;
        $y = defined( $defaultCodonStats{$pos}{$right} ) ? $defaultCodonStats{$pos}{$right} : 0;
    }

    return $x >= $y;
}

sub comparePosLeftRightGE($$$$$) {
    my ( $posL, $posR, $group, $codon, $def ) = ( @_[0 .. 4] );
    my ( $x, $y ) = ( 0, 0 );
    if ( !defined($group) || $group eq NIL || $group eq '' ) {
        $x = defined( $codonStats{$posL}{$codon} ) ? $codonStats{$posL}{$codon} : 0;
        $y = defined( $codonStats{$posR}{$codon} ) ? $codonStats{$posR}{$codon} : 0;
    } else {
        $x = defined( $codonStats{$group}{$posL}{$codon} ) ? $codonStats{$group}{$posL}{$codon} : 0;
        $y = defined( $codonStats{$group}{$posR}{$codon} ) ? $codonStats{$group}{$posR}{$codon} : 0;
    }

    if ( $x == $y && $x == 0 ) {
        $x = $def eq 'L' ? 1 : 0;
        $y = $def eq 'L' ? 0 : 1;
    }
    return $x >= $y;
}

$PROG      = basename( $0, '.pl' );
$firstFile = $ARGV[0];
$group     = '';
%groups    = ();
%sequences = ();
$/         = ">";
while ( $record = <> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/, $record );
    $id       = shift(@lines);
    $sequence = lc( join( '', @lines ) );
    $length   = length($sequence);

    if ( $length == 0 ) {
        next;
    } elsif ( $length % 3 != 0 ) {
        die("$PROG:\tNot a codon alignment! See:\n\t$firstFile\n");
    } else {
        $sequences{$id} = $sequence;

        if ( $numberSelected > 0 ) {
            @values      = split( /\Q$delim\E/, $id );
            $numberFound = scalar(@values);
            if ( $numberFound < $numberSelected ) {
                die("$0 ERROR: non-existant field specified. Wanted $numberSelected but found $numberFound\n");
            }
            $group = join( $delim, ( @values[@fields] ) );
        }

        $group = $numberSelected > 0 ? $group : NIL;
        $groups{$id} = $group;
    }
}

foreach $id ( keys(%sequences) ) {
    $sequence = $sequences{$id};
    $seqLimit = length($sequence) - 2;                          # TO-DO: what about the second to last insertion opportunity?
    $group    = defined( $groups{$id} ) ? $groups{$id} : NIL;

    if ( defined( $inserts{$id} ) ) {

        # position is 1 based
        foreach $pos ( keys( %{ $inserts{$id} } ) ) {
            if ( $pos > $seqLimit ) {

                # should be in frame if at length
                if ( $pos == length($sequence) ) { print TABL $id, "\t", $pos, "\t", uc( $inserts{$id}{$pos} ), "\n"; }
                next;
            }

            $insert = $inserts{$id}{$pos};
            $iDivis = length($insert) % 3;
            $iFrame = $pos % 3;
            if ( ( ( $iDivis == 0 ) && ( $iFrame != 0 ) ) ) {

                # zero based codon number
                $codonNumber = int( ( $pos - 1 ) / 3 );
                $codonStart  = $codonNumber * 3;

                # A2 insertion
                if ( $iFrame == 2 ) {
                    $A2L2 = substr( $insert,   -2 ) . substr( $sequence, $codonStart + 2, 1 );
                    $A2R1 = substr( $sequence, $codonStart, 2 ) . substr( $insert, 0, 1 );

                    if ( compareLeftRightGE( $A2R1, $A2L2, $group, $codonNumber ) ) {

                        #if ( $codonStats{$codonNumber}{$A2R1} >= $codonStats{$codonNumber}{$A2L2} ) {
                        #A2::R1 shift
                        $newInsert = substr( $insert, 1 ) . substr( $sequence, $pos, 1 );
                        $newCodon  = $A2R1;
                        $newPos    = $codonStart + 3;                                       # +2 +1 => +3 given 1-based
                    } else {

                        #A2::L2 shift
                        $newInsert = substr( $sequence, $codonStart, 2 ) . substr( $insert, 0, -2 );
                        $newCodon  = $A2L2;
                        $newPos    = $codonStart;                                           # -1 + 1 => + 0 given 1-based
                    }

                    # A1 insertion
                } else {
                    $A1L1 = substr( $insert,   -1 ) . substr( $sequence, $codonStart + 1, 2 );
                    $A1R2 = substr( $sequence, $codonStart, 1 ) . substr( $insert, 0, 2 );
                    if ( compareLeftRightGE( $A1L1, $A1R2, $group, $codonNumber ) ) {

                        #if ( $codonStats{$codonNumber}{$A1L1} >= $codonStats{$codonNumber}{$A1R2} ) {
                        #A1::L1 shift
                        $newInsert = substr( $sequence, $codonStart, 1 ) . substr( $insert, 0, -1 );
                        $newCodon  = $A1L1;
                        $newPos    = $codonStart;    # -1 + 1 => + 0 given 1-based
                    } else {

                        #A1::R2 shift
                        $newInsert = substr( $insert, 2 ) . substr( $sequence, $pos, 2 );
                        $newCodon  = $A1R2;
                        $newPos    = $codonStart + 3;                                       # +2 +1 => +3 given 1-based
                    }
                }

                print TABL $id, "\t", $newPos, "\t", uc($newInsert), "\n";
                substr( $sequence, $codonStart, 3 ) = $newCodon;
            } else {
                print TABL $id, "\t", $pos, "\t", uc( $inserts{$id}{$pos} ), "\n";
            }
        }
    }

    # perform after insertion corrections
    while ( $sequence =~ /([A-Za-z]{3})((---)+)([A-Za-z]{3})/g ) {
        ( $left, $gaps, $right ) = ( $1, $2, $4 );

        $frame = $-[1] % 3;

        if ( $frame == 1 ) {
            $pivotCodon = substr( $left, -1 ) . substr( $right, 0, 2 );
            $leftCN     = int( $-[2] / 3 );
            $rightCN    = int( ( $+[2] - 1 ) / 3 );

            if ( comparePosLeftRightGE( $leftCN, $rightCN, $group, $pivotCodon, 'R' ) ) {

                # LEFT SHIFT, move 2 to the left
                $replacement = $left . substr( $right, 0, 2 ) . $gaps;
                substr( $sequence, $-[1], length($replacement) ) = $replacement;
            } else {

                # RIGHT SHIFT, move 1 to right
                $replacement = substr( $left, 0, -1 ) . $gaps . substr( $left, -1 );
                substr( $sequence, $-[1], length($replacement) ) = $replacement;
            }
        } elsif ( $frame == 2 ) {
            $pivotCodon = substr( $left, -2 ) . substr( $right, 0, 1 );
            $leftCN     = int( $-[2] / 3 );
            $rightCN    = int( ( $+[2] - 1 ) / 3 );

            if ( comparePosLeftRightGE( $leftCN, $rightCN, $group, $pivotCodon, 'L' ) ) {

                # LEFT SHIFT, move 1 to the left
                $replacement = $left . substr( $right, 0, 1 ) . $gaps;
                substr( $sequence, $-[1], length($replacement) ) = $replacement;
            } else {

                # RIGHT SHIFT, move 2 to right
                $replacement = substr( $left, 0, -2 ) . $gaps . substr( $left, -2 );
                substr( $sequence, $-[1], length($replacement) ) = $replacement;
            }
        }
    }

    print '>', $id, "\n", $sequence, "\n";
}

if ($outputTable) {
    close(TABL);
}
