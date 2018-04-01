[![Build Status](https://travis-ci.org/tseemann/samclip.svg?branch=master)](https://travis-ci.org/tseemann/samclip) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# samclip

Filter SAM file for soft and hard clipped alignments

## Introduction

Most short read aligners perform 
[*local alignment*](https://en.wikipedia.org/wiki/Sequence_alignment#Global_and_local_alignments)
of reads to the reference genome.
Examples includes `bwa mem`, `minimap2`, and `bowtie2` (unless in `--end-to-end` mode).
This means the *ends* of the read may not be part of the best alignment.

This can be caused by:
* adapter sequences (aren't in the reference)
* poor quality bases (mismatches only make the alignment score worse)
* structural variation in your sample compared to the reference
* reads overlapping the start and end of contigs (includign circular genomes)

Read aligners output a [SAM file](https://en.wikipedia.org/wiki/SAM_(file_format)).
Column 6 in this format stores the 
[CIGAR string](https://www.drive5.com/usearch/manual/cigar.html).
which describes which parts of the read aligned and which didn't.
The unaligned ends of the read can be "soft" or "hard" clipped, 
denoted with `S` and `H` at each end of the CIGAR string. 
It is possible for both types to be present, but that is not common.
Soft and hard don't mean anything biologically, they just refer
to whether the full read sequence is in the SAM file or not.

Some examples of a 100bp aligned read:
* `100M` - the read is fully aligned
* `30S70M` - soft clipped 30 bp on left end, followed by 70 aligned bp.
* `95M5H` - 95 aligned bp, followed by 5 bp unaligned (hard clipped).
* `25S40M35S` - only the middle 40bp of the read aligned

## Motivation

One may wish to remove these alignments to avoid downstream problems.
In particular, `samclip` was designed to remove clipped alignments
to improve variant calling, by removing suspicious local aligments
causing false positives near structural variation. However, it does
keep them if they hit the ends of contigs, which is particularly
important given the lower coverage often observed at those locations.

## Installation

`samclip` has no dependencies except [Perl 5](https://www.perl.org/). 
It uses no non-core modules. Hooray!

### Direct script download
```
% cd /usr/local/bin  # choose a folder in your $PATH
% wget https://raw.githubusercontent.com/tseemann/samclip/master/samclip
% chmod +x samclip
```
### Homebrew
```
% brew install brewsci/bio/samclip
```
### Conda
```
% conda install -c bioconda -c conda-forge samclip
```
### Github
```
% git clone https://github.com/tseemann/samclip.git
% cp samclip/samclip /usr/local/bin # choose a folder in your $PATH
```

## Test Installation

```
% ./samclip --version
samclip 0.1.0

% ./samclip --help
Usage: samclip.pl --ref ref.fasta [--max=5] [-v/--invert] < in.sam > out.sam
```

## Usage
```
% samclip --ref ref.fa < in.sam > out.sam

% samtools view in.bam | samclip --ref ref.fa | samtools sort > out.bam

% bwa mem ref.fa R1.fq R2.fq | samclip --ref ref.fa | samtools sort > out.bam 
```

## Options

* `--ref FILE` should be a FASTA file indexed with `samtools faidx FILE`
* `--max INTEGER` is the maximum soft+hard clipping to allow. Set to 0 to allow none.
* `-v` or `--invert` will output the records that would have clipped and discard the good ones
* `--debug` is verbose debugging information for testing purposes

## Issues

Submit feedback to the [Issue Tracker](https://github.com/tseemann/samclip/issues).

## License

[GPL v3](https://raw.githubusercontent.com/tseemann/samclip/master/LICENSE)

## Author

[Torsten Seemann](http://tseemann.github.io/)
