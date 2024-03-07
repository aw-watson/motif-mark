# Motif-mark

## Purpose

Given a list of sequences (FASTA format), and a line-by-line list of motifs, draw all sequences to scale, and annotate all occurrences of each motif. Expects sequences to be formatted with lowercase letters denoting intron sequence,and uppercase letters denoting exon sequence. Motifs can use ambiguous nucleotide notation, and all possible matches will be annotated in the final drawing. 

## Environment

```motif-mark``` is built with python 3.11, and uses ```pycairo``` 1.23 to draw its annotations. An environment containing those installations should be able to run motif-mark.

## Usage

```motif-mark``` is useable from the command line, and expects two arguments:
* -f, followed by the path to the FASTA file containing the sequences to annotate.
* -m, followed by the path to a text file containing the motifs to annotate on all sequences under consideration.
