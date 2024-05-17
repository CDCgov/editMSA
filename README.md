# Edit MSA

A collection of scripts for editing multiple sequence alignments, usually in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format. Some functionality may exist for downstream processes and may not be generally useful.

POC: Sam Shepard (<vfn4@cdc.gov>), CDC/NCIRD/ID/OD

## Overview

| Script | Author | Description | Notes |
| ---------------------- | ------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------- |
| *codonCorrect.pl* | S. S. Shepard | Corrects insertions (given by the [TSV](https://en.wikipedia.org/wiki/Tab-separated_values)) and deletions to be in-frame using the MSA stats. | Deprecated. |
| *codonCorrectStats.pl* | S. S. Shepard | Corrects insertions (given by the TSV) and deletions to be in-frame using user-supplied statistics in Perl [storable format](https://perldoc.perl.org/Storable). | Used by dais-ribosome (GitHub pending). |
| *displayStats.pl* | S. S. Shepard | Useful for displaying the codon-weight matrix used by `codonCorrectStats.pl`. | Uses Perl storable format. |
| *storeStats.pl* | S. S. Shepard | Takes an MSA and produces a codon-weight matrix in Perl storable format. | Format used by `codonCorrectStats.pl` |
| *fillSequences.pl* | S. S. Shepard | Re-inserts insertions from a TSV (ID / Position / Insert) back into the MSA. | |
| *ordinalHeaders.pl* | S. S. Shepard | Changes FASTA headers to an ordinal `S#`. Has various options for annotating. | |
| *removeGapColumns.pl* | S. S. Shepard | Various options for dealing with gap columns in an MSA, including K-tons. | Overwrites in-place |
| *reviseTaxa.pl* | S. S. Shepard | Revises annotation taxa in an MSA. | A component of the [LABEL MS](https://wonder.cdc.gov/amd/flu/label) |
| *shiftAlignment.pl* | S. S. Shepard | Shifts alignments to the left or right at a particular region. Useful for gaps before manual editing. ||
| *sortFASTA.pl* | S. S. Shepard | Sorts the FASTA file by the header or a designatd set of fields. | |
| *stripSequences.pl* | S. S. Shepard | Strips unwanted characters from a FASTA file. | FASTA may be unaligned after this operation. |
| *treeOrderedFasta.pl* | S. S. Shepard | Orders the FASTA file using a provided Newick where IDs or headers match. ||
| *diceAlignment.pl* | S. S. Shepard | Chops alignment at specified region into a de-duplicated MSA for editing and downstream splicing. | Used in tandem with `spliceAlignment.pl` |
| *spliceAlignment.pl* | S. S. Shepard | Takes the manually edited miniNT or miniAA file and splices it back into a full nucleotide alignment. | Used in tandem with `diceAlignment.pl` |

## Usage

For up-to-date usages, run `perl <script.pl>` with no arguments. Example usage for each script:

```bash
# codonCorrect.pl
Usage:
        perl codonCorrect.pl <nts.fasta> [options]
                --insertion-table|-I <STR>      Insertion table for insertion corrections.
                --output-table|-O <STR>         Output file for the insertion table.

# codonCorrectStats.pl
Usage:
        perl codonCorrectStats.pl <nts.fasta> [options]
                --stats|-S <in.sto>             Input file for storable object. Default: none
                --delim|-D <CHAR>               Delimiter for header fields.
                --field|-F <STR>                Comma-delimited set of fields to use for group. Default: no group.
                --insertion-table|-I <STR>      Insertion table for insertion corrections.
                --output-table|-O <STR>         Output file for the insertion table.

# diceAlignment.pl
Usage:
        perl diceAlignment.pl <in.fasta> <...>
                -P|--prefix <STR>       Prefix used for output.
                -B|--begin <INT>        Begin coordinate.
                -E|--end <INT>          End coordinate.
                -R|--shift-right        Shift alignment right.
                -T|--translate          Translate data as well.
                -C|--sort-by-count      Sort data by their cluster counts.

# displayStats.pl
Usage:
        perl displayStats.pl <STATS.sto> [options]
                --verbose|V     Make display verbose.

# fillSequences.pl
Usage:
        perl fillSequences.pl <fasta> <ins.txt> [-F|--fill-triplets]

# ordinalHeaders.pl
Usage:
        perl ordinalHeaders.pl <input.fasta> [-A] [-H] [-N <STR>] [-O <STR>]
        
# removeGapColumns.pl
Usage:
        perl removeGapColumns.pl <file.fasta> [options]
                -D|-display-k-tons      Displays sequences where for columns with 0 < bases <= K.
                -T|-tab-formatted       Print out tab-delimited data.
                -S|-show-column-support Show column support.

# reviseTaxa.pl
Usage:
        perl reviseTaxa.pl <input.fasta> [OPTIONS] [-M <file1 file2 ...>]
                --confirm-prediction|-C         Confirm predicted annotations.
                --delete-prev|-D                Delete previous annotations where there are two.
                --delete-single|-S              Delete previous annotations where there is one, filters AFTER delete-prev.
                --find|-F <TEXT>                Select sequences including this annotation.
                --replace|-R <TEXT>             Replace the annotation with TEXT.
                --add-annot|-A <FILE>           Add annotations based on tab delimited file (ID ANNOT).
                --ignore-fasta-annot|-G <FILE>  Ignores previous annotation on FASTA headers vs. annotation file.
                --fuzzy-match|-Z                Searches for IDs in FASTA (-A option), matching if header contains the ID.
                --order-mode|-O <out.file>      Annotation with ordinals for truncated names.
                --match-files|-M <file1 ...>    Output filenames containing annotation names in the input file.
                --previous-infix|N              Suffix and prefix only apply to a 'previous' or first of a doublet annotation.
                --prefix|-P <TEXT>              Prefix for matching filenames with annotations OR adds prefix to header.
                --suffix|-X <TEXT>              Suffix for matching filenames with annotations OR adds suffix to header.
                --in-place|-I                   Revise files in-place (could be dangerous).
                --join-to-end|-J <TEXT>         Join annotation to end of the header (similar to add but without a file).
                --last-field|-L <delim>         Clips the last field of the header and turns it into an annotation. Uses the specified delimiter to determine fields.
                --append-pipe-annot|-p          Appends annotation as a header field with pipe delim.

# shiftAlignment.pl
Usage:
        perl shiftAlignment.pl <in.fasta> <...>
                -B|--begin <INT>        Begin coordinate.
                -E|--end <INT>          End coordinate.
                -R|--shift-right        Shift alignment right.

# sortFASTA.pl
Usage:
        perl sortFASTA.pl [FASTA ...] [OPTIONS]
                -C|--case-sensitive     Ignore case.
        <STDIN> read if no FASTA given.

# spliceAlignment.pl
Usage:
        perl spliceAlignment.pl <diced> <miniNT> [options]
                -F|--fill-missing               Fill missing data with gaps rather than excluding.
                -A|--amino-acid-msa <FILE>      Amino acid alignment.
                -P|--prefix <STR>               Prefix used for AAtoNT file output.
                -S|--splice-now                 Output the splice without further editing of the AAtoNT file.

# storeStats.pl
Usage:
        perl storeStats.pl <nts.fasta> [options]
                --output|-O <file.sto>  Output file for storable object. Default: STDOUT
                --delim|-D <CHAR>       Delimiter for header fields. Default delim is '|'.
                --field|-F <STR>        Comma-delimited set of fields to use for group. Default: no group.

# stripSequences.pl
Usage:
        stripSequences.pl <file.fas> {-N|-L|<quoted_characters_to_delete>}
                -N|--remove-bad-bases   Removes invalid nucleotide characters from the sequence.
                -L|--strip-lower        Removes lowercase letters from the sequence.
                -F|--fix-header         Removes and replaces troublesome characters from the FASTA header.

# treeOrderedFasta.pl
Usage:
        treeOrderedFasta.pl <tree.nwk> <sequences.fasta> <output.fasta> [-Q|--quoted-ids]
```

If you spot a bug in the usage, please file an issue or make a PR to this repo.

## Installation

1. Install [Perl 5](https://www.perl.org/get.html). If you are running MacOS or Linux *you likely already have it*. No, you don't need a container, package manager, *venv*, or Conda. Perl is aggressively backwards compatible and ubiquitous.
2. [Clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) this repo.

## Notices

Copyright is public domain but [attributions to the original authors are welcomed](CITATION.bib). The software license uses [ASL2](http://www.apache.org/licenses/LICENSE-2.0.html).

### Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not subject to domestic copyright protection under 17 USC § 105. This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the CC0 1.0 Universal public domain dedication. All contributions to this repository will be released under the CC0 dedication. By submitting a pull request you are agreeing to comply with this waiver of copyright interest.

### License Standard Notice

The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see <http://www.apache.org/licenses/LICENSE-2.0.html>

The source code forked from other open source projects will inherit its license.

### Contributing Standard Notice

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

### Records Management Standard Notice

This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).
