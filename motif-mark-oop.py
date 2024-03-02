#!/usr/bin/env python

import argparse
import re
import os
import cairo
import math
from queue import PriorityQueue

base_mapping = {"A":"A", "C":"C","G":"G","T":"T",
                "W":"[AT]","S":"[CG]","M":"[AC]","K":"[GT]",
                "R":"[AG]","Y":"[CT]","B":"[CGT]","D":"[AGT]",
                "H":"[ACT]","V":"[ACG]","N":"[AGCT]", "U":"T"}

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

class Canvas:
    def __init__(self, n_seq:int):
        self._seq_amt = float(n_seq)
        self._surface = cairo.ImageSurface(cairo.FORMAT_RGB16_565, 
                                           2300, 
                                           50+150*math.ceil(self._seq_amt/2.0))
        self._cx = cairo.Context(self._surface)
        #white background
        self._cx.set_line_width(1)
        self._cx.set_source_rgba(1,1,1,1)
        self._cx.rectangle(0,0,2300,50+150*math.ceil(self._seq_amt/2.0))
        self._cx.fill()
        self._xoffset = 0
        self._yoffset = 0
    def draw_seq(self):
        #base position (leftmost position we want to draw at)
        base_position = (100+ self._xoffset*1100, 50 + self._yoffset*150) #(x,y)
        #draw intron/exon

        #draw motifs
        #increment offset to draw next sequence
        if self._xoffset == 0 and self._yoffset == math.ceil(self._seq_amt/2.0):
            self.y_offset = 0
            self.x_offset = 1
        else:
            self._yoffset += 1
        raise NotImplementedError
    
class _Motif:
    def __init__(self, motif_name:str, motif_start:int, motif_end:int):
        self.name = motif_name
        self.start = motif_start
        self.end = motif_end

class _Exon:
    def __init__(self, exon_start:int, exon_end:int):
        self.start = exon_start
        self.end = exon_end

class Sequence:
    def __init__(self, name:str,  nucleotides:str):
        self.name = name
        seq_match = re.fullmatch("[actg]*([ACTG]+)[actg]*", nucleotides)
        self.exon = None
        if seq_match is None:
            raise ValueError("Innapropriately formatted reference sequence: expects a single stretch of \
                             uppercase DNA nucleotides (denoting an exon) \
                             flanked by surrounding intron sequence (in lowercase).")
        else:
            self.exon = _Exon(seq_match.start(1), seq_match.end(1))
        self.sequence:str = nucleotides.upper()
        self.motifs:list[_Motif] = []
    def find_motif(self, motif_seq:str):
        #regex_maker accounts for degenerate nucleotide code
        motif_pattern:str = "(?=(" + regex_maker(motif_seq) + "))"
        #looking for all overlapping matches
        motif_iter = re.finditer(motif_pattern, self.sequence)
        for m in motif_iter:
            #add all matches found to the motifs associated with this sequence
            self.motifs.append(_Motif(motif_seq, m.start(1), m.end(1)))


def get_args():
    parser = argparse.ArgumentParser(description = "Find and mark motifs in a given genomic sequence")
    parser.add_argument("-f", "--fastafile", help="Path to FASTA file containing sequences of interest", required = True)
    parser.add_argument("-m", "--motifs", help="Path to file of motifs to annotate (supports ambiguous nucleotide notation)", required = True)
    return parser.parse_args()

if __name__ == "__main__":
    #Read in arguments
    args = get_args()
    #Put sequences of interest on one line
    temp_filename = "oneline_" + args.fastafile
    oneline_fasta(args.fastafile, temp_filename)
    #parse fasta file into sequence objects
    seq_list:list[Sequence] = []
    with open(temp_filename, mode = 'rt') as f:
        for line in f:
            seq = Sequence(line.strip()[1:], f.readline().strip())
            seq_list.append(seq)
    os.remove(temp_filename)

    #Find all motifs in all sequences
    with open(args.motifs, mode='rt') as f:
        for line in f:
            for s in seq_list:
                s.find_motif(line.strip().upper())


    raise NotImplementedError()