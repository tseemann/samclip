# samclip

Filter SAM file for soft and hard clipped alignments

## Introduction

TODO: description of the problem

## Install
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

## License

GPLv3

## Author

Torsten Seemann
