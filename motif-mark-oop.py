#!/usr/bin/env python

import argparse
import re

base_mapping = {"A":"A", "C":"C","G":"G","T":"T",
                "W":"[AT]","S":"[CG]","M":"[AC]","K":"[GT]",
                "R":"[AG]","Y":"[CT]","B":"[CGT]","D":"[AGT]",
                "H":"[ACT]","V":"[ACG]","N":"[AGCT]"}

def oneline_fasta(infile: str, outfile: str):
    """Takes a FASTA file (infile) and puts all sequence lines on one line. Writes the result to outfile"""
    with open(infile, 'rt') as rf, open(outfile, 'wt') as wf:
        seq: str = ''
        first_record: bool = True
        while True:
            line = rf.readline().strip()
            if not line:
                break
            if line.startswith(">"):
                if first_record:
                    wf.write(line + "\n")
                    first_record = False
                else:
                    wf.write(seq + "\n")
                    seq = ""
                    wf.write(line + "\n")
            else:
                seq += line
        wf.write(seq)

def regex_maker(motif_seq: str) -> str: 
    motif_regex = ""
    for b in motif_seq:
        motif_regex += base_mapping[b]
    return motif_regex

class _Motif:
    def __init__(self, motif_name:str, motif_start:int, motif_end:int):
        self.name = motif_name
        self.start = motif_start
        self.end = motif_end

class Sequence:
    def __init__(self, nucleotides:str):
        seq_match = re.fullmatch("[actg]*([ACTG]+)[actg]*", nucleotides)
        self.exon_start:int = 0
        self.exon_end:int = 0
        if seq_match is None:
            raise ValueError("Innapropriately formatted input: expects a single stretch of uppercase nucleotides \
                             (denoting an exon) \
                             flanked by surrounding intron sequence (in lowercase).")
        else:
            self.exon_start, self.exon_end = seq_match.span(1)
        self.sequence:str = nucleotides.upper()
        self.motifs:list[_Motif] = []
    def find_motif(self, motif_seq:str, motif_name:str):
        #regex_maker accounts for degenerate nucleotide code
        motif_pattern:str = "(?=(" + regex_maker(motif_seq) + "))"
        #looking for all overlapping matches
        motif_iter = re.finditer(motif_pattern, self.sequence)
        for m in motif_iter:
            #add all matches found to the motifs associated with this sequence
            self.motifs.append(_Motif(motif_name, m.start(), m.end()))


if __name__ == "__main__":
    raise NotImplementedError()