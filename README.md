**As a first step, this document is under governance review. When the review completes as appropriate per local and agency processes, the project team will be allowed to remove this notice. This material is draft.**

POC: Sam Shepard (vfn4@cdc.gov), CDC/DDID/NCIRD/ID/OD

# Edit MSA

A collection of scripts for editing multiple sequence alignments, usually in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format. Some functionality may exist for downstream processes and may not be generally useful.

## Overview

| Script                 | Author        | Description                                                                                                                                                      | Notes                                                               |
| ---------------------- | ------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------- |
| *codonCorrect.pl*      | S. S. Shpeard | Corrects insertions (given by the [TSV](https://en.wikipedia.org/wiki/Tab-separated_values)) and deletions to be in-frame using the MSA stats.                   | Deprecated.                                                         |
| *codonCorrectStats.pl* | S. S. Shepard | Corrects insertions (given by the TSV) and deletions to be in-frame using user-supplied statistics in Perl [storable format](https://perldoc.perl.org/Storable). | Used by dais-ribosome (GitHub pending).                             |
| *displayStats.pl*      | S. S. Shepard | Useful for displaying the codon-weight matrix used by `codonCorrectStats.pl`.                                                                                    | Uses Perl storable format.                                          |
| *storeStats.pl*        | S. S. Shepard | Takes an MSA and produces a codon-weight matrix in Perl storable format.                                                                                         | Format used by `codonCorrectStats.pl`                               |
| *fillSequences.pl*     | S. S. Shepard | Re-inserts insertions from a TSV (ID / Position / Insert) back into the MSA.                                                                                     |                                                                     |
| *ordinalHeaders.pl*    | S. S. Shepard | Changes FASTA headers to an ordinal `S#`. Has various options for annotating.                                                                                    |                                                                     |
| *removeGapColumns.pl*  | S. S. Shepard | Various options for dealing with gap columns in an MSA, including K-tons.                                                                                        | Overwrites in-place                                                 |
| *reviseTaxa.pl*        | S. S. Shepard | Revises annotation taxa in an MSA.                                                                                                                               | A component of the [LABEL MS](https://wonder.cdc.gov/amd/flu/label) |
| *shiftAlignment.pl* | S. S. Shepard | Shifts alignments to the left or right at a particular region. Useful for gaps before manual editing. ||
| *sortFASTA.pl*         | S. S. Shepard | Sorts the FASTA file by the header or a designatd set of fields.                                                                                                                              |                                                                     |
| *stripSequences.pl* | S. S. Shepard | Strips unwanted characters from a FASTA file. | FASTA may be unaligned after this operation. |
| *treeOrderedFasta.pl* | S. S. Shepard | Orders the FASTA file using a provided Newick where IDs or headers match. ||
| *diceAlignment.pl*     | S. S. Shepard | Chops alignment at specified region into a de-duplicated MSA for editing and downstream splicing.                                                                | Used in tandem with `spliceAlignment.pl`                            |
| *spliceAlignment.pl*   | S. S. Shepard | Takes the manually edited miniNT or miniAA file and splices it back into a full nucleotide alignment.                                                            | Used in tandem with `diceAlignment.pl`                              |

## Notices

Copyright is public domain but [attributions to the original authors are welcomed](CITATION.bib). The software license uses ASF2 with additional notices and a [standard DISCLAIMER](DISCLAIMER.md).

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
