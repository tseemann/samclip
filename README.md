# samclip

Filter SAM file for soft and hard clipped alignments

## Introduction

TODO: description of the problem

## Install
`samclip` has no dependencies except Perl 5.  It uses no non-core modules.
```
% wget https://raw.githubusercontent.com/tseemann/samclip/master/samclip

% chmod +x samclip

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

TODO: explain these options

* `--ref FILE`
* `--max INTEGER`

## Issues

Post to https://github.com/tseemann/samclip/issues

## License

GPLv3

## Author

Torsten Seemann

