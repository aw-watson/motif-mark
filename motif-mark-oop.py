#!/usr/bin/env python

import argparse
import re
import os
import cairo
from queue import PriorityQueue
from typing import Iterator

base_mapping:dict[str,str] = {"A":"A", "C":"C","G":"G","T":"T",
                "W":"[AT]","S":"[CG]","M":"[AC]","K":"[GT]",
                "R":"[AG]","Y":"[CT]","B":"[CGT]","D":"[AGT]",
                "H":"[ACT]","V":"[ACG]","N":"[AGCT]", "U":"T"}

colors:list[tuple[float,float,float,float]] = [(1,38/255,0,1),
          (1,159/255,0,1),
          (1,221/255,0,1),
          (159/255,180/255,0,1),
          (50/255,130/255,0,1),
          (0,130/255,90/255,1),
          (0,77/255,132/255,1),
          (81/255,0,32/255,1),
          (72/255,38/255,13/255,1),
          (210/255,180/255,140/255,1)]

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
    motif_regex:str = ""
    for b in motif_seq:
        motif_regex += base_mapping[b]
    return motif_regex
    
class _Motif:
    def __init__(self, motif_name:str, motif_start:int, motif_end:int):
        self.name:str = motif_name
        self.start:int = motif_start
        self.end:int = motif_end

class _Exon:
    def __init__(self, exon_start:int, exon_end:int):
        self.start:int = exon_start
        self.end:int = exon_end

class Sequence:
    def __init__(self, name:str,  nucleotides:str):
        self.name:str = name
        seq_match: None | re.Match[str] = re.fullmatch("[actg]*([ACTG]+)[actg]*", nucleotides)
        self.exon:_Exon = _Exon(0,0)
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
        motif_iter:Iterator = re.finditer(motif_pattern, self.sequence)
        for m in motif_iter:
            #add all matches found to the motifs associated with this sequence
            self.motifs.append(_Motif(motif_seq, m.start(1), m.end(1)))

class Canvas:
    def __init__(self, n_seq:int):
        self._surface:cairo.ImageSurface = cairo.ImageSurface(cairo.FORMAT_RGB16_565, 
                                           1200, 
                                           50+250*n_seq)
        self._cx:cairo.Context = cairo.Context(self._surface)
        #white background
        self._cx.set_line_width(1)
        self._cx.set_source_rgba(1,1,1,1)
        #imagine a 1-column grid of 1000x200 rectangles, with 100 pixels of horizontal padding and 50 pixels of vertical padding
        self._cx.rectangle(0,0,1200,50+250*n_seq)
        self._cx.fill()
        self._cx.set_source_rgba(0,0,0,1)
        self._yoffset:int = 0
    def output(self, name:str):
        self._surface.write_to_png(f"{name}.png")
    def draw_seq(self, seq:Sequence):
        #base position (top left corner of the "rectangle" we're drawing our sequence in)
        base_position:list[int] = [100, 50 + self._yoffset*250] #(x,y)
        #write out sequence name from FASTA header
        self._cx.set_source_rgba(0,0,0,1)
        self._cx.move_to(base_position[0], base_position[1])
        self._cx.show_text(seq.name)
        #draw intron
        self._cx.move_to(base_position[0], base_position[1]+100) #halfway down the left side
        self._cx.set_line_width(2)

        self._cx.line_to(base_position[0] + len(seq.sequence), base_position[1] + 100) #straight line across
        self._cx.stroke()
        #draw exon
        self._cx.move_to(base_position[0] + seq.exon.start, base_position[1]+100) #halfway down the left side, skip intron bases
        self._cx.set_line_width(10)
        self._cx.line_to(base_position[0] + seq.exon.end + 1, base_position[1] + 100) #straight line
        self._cx.stroke()
        #draw motifs
        self.draw_motifs(seq, base_position)
        #increment offset to draw next sequence
        self._yoffset += 1
    def draw_motifs(self, seq:Sequence, base_position:list[int]):
        #staggers to account for overlaps

        #setup
        motifs_present:set[str]  = set()
        m_todraw:PriorityQueue[tuple[int,_Motif]] = PriorityQueue()

        for m in seq.motifs:
            motifs_present.add(m.name) #track all motif types in sequence
            m_todraw.put((m.end,m)) #place motifs in priority queue based on end position
        
        #start drawing
        #set line size
        self._cx.set_line_width(10)
        stagger_ctr = 0
        #draw until all motifs are drawn
        while m_todraw.qsize() != 0: #not multithread safe
            m_tostagger:PriorityQueue[tuple[int,_Motif]] = PriorityQueue()
            while m_todraw.qsize() != 0: #run through queue once
                #get motif with earliest end position
                m = m_todraw.get()[1]
                #put all overlapping motifs in m_tostagger
                while m_todraw.qsize() != 0 and m_todraw.queue[0][1].start < m.end:
                    m_tostagger.put(m_todraw.get())
                #draw motif
                self._cx.set_source_rgba(*motif_color_map[m.name])
                self._cx.move_to(base_position[0] +  m.start,base_position[1] + 100 + 12*stagger_ctr)
                self._cx.line_to(base_position[0] + m.end + 1, base_position[1] + 100 + 12*stagger_ctr)
                self._cx.stroke()
            #add offset
            stagger_ctr += 1
            #feed all undrawn motifs back into m_todraw
            m_todraw = m_tostagger

        #legend
        self._cx.set_source_rgba(0,0,0,1)
        self._cx.set_line_width(2)
        for i, m in enumerate(motifs_present):
            legend_pos = [base_position[0] + 50, base_position[1] + 20 + 10*i]
            self._cx.set_source_rgba(0,0,0,1)
            self._cx.rectangle(legend_pos[0],legend_pos[1], 6, 6)
            self._cx.stroke()
            self._cx.set_source_rgba(*motif_color_map[m])
            self._cx.rectangle(legend_pos[0], legend_pos[1], 5 ,5)
            self._cx.fill()
            self._cx.set_source_rgba(0,0,0,1)
            self._cx.move_to(legend_pos[0] + 10,legend_pos[1] + 6)
            self._cx.show_text(f": {m}")


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
    n_motifs = 0
    motif_color_map = {}
    with open(args.motifs, mode='rt') as f:
        for line in f:
            motif = line.strip().upper()
            if motif not in motif_color_map:
                motif_color_map[motif] = colors[n_motifs]
                n_motifs += 1
            for s in seq_list:
                s.find_motif(motif)

    #draw sequences
    cnv: Canvas = Canvas(len(seq_list))
    for s in seq_list:
        cnv.draw_seq(s)
                
    #output file
    cnv.output("testing")
