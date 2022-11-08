# -*- coding: utf-8 -*-
"""
    utils.genomic_region
    ~~~~~~~~~~~~

    Genomic region definition.

    @Copyright: (c) 2017 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""
import copy
import pybedtools
import os, sys

from tenxtools.utils import sequence
from tenxtools.utils import constants

class GenomicRegion(object): # [start, end]; length = end - start + 1

    def __init__(self, chrom=None, pos1=None, pos2=None, strand=None, name=None, code_name=None, full_name=None):
        if pos1 and pos1 == pos2: raise Exception('In class GenomicRegion, pos1 %s cannot equal to pos2 %s' % (str(pos1), str(pos2)))
        self.chrom = chrom
        try: # [start, end] closed interval
            self.start = min(int(pos1), int(pos2))
            self.end = max(int(pos1), int(pos2))
        except Exception as e:
            self.start = min(pos1, pos2)
            self.end = max(pos1, pos2)
        self.strand = strand
        self.name = name
        self.code_name = code_name
        self.full_name = full_name


    def slice(self, pos1, pos2, new_strand, code_name, full_name):
        new_region = copy.copy(self)
        new_region.start = min(pos1, pos2)
        new_region.end = max(pos1, pos2)
        new_region.strand = new_strand
        new_region.code_name = code_name
        new_region.full_name = full_name
        return new_region

    def get_3prime_pos(self):
        if self.strand == '+':      return self.end
        if self.strand == '-':      return self.start
    def get_5prime_pos(self):
        if self.strand == '+':      return self.start
        if self.strand == '-':      return self.end

    @property
    def bp_length(self):
        return self.end - self.start + 1

    def __eq__(self, other):
        return self.chrom == other.chrom and self.start == other.start and self.end == other.end and self.name == other.name and self.strand == other.strand and self.code_name == other.code_name and self.full_name == other.full_name

    ####
    # representation #
    ####
    def to_fasta_header(self):
        return '>%s::%s::%s:%s-%s' % (self.name, self.strand, self.chrom, str(self.start), str(int(self.end)+1))

    def to_bed(self):
        if self.strand == '+':  start, end = self.start, int(self.end) + 1
        else:                   start, end = int(self.start) - 1, self.end
        return '\t'.join(str(v) if v else '' for v in [self.chrom, start, end, self.name])
    def __repr__(self):
        if not self.start: return '%s %s %s ' % (' ','  ','  ')
        code_name = str(self.code_name) if self.code_name else ''
        if self.strand == '+':
            return '%s[%2s,%2s]%s' % (str(self.strand), str(self.start), str(self.end), code_name)
        if self.strand == '-':
            return '%s[%2s,%2s]%s' % (str(self.strand), str(self.end), str(self.start), code_name)
    def __str__(self):
        return '%6s' % (str(self.code_name) if self.code_name else '')

    def to_stat(self):
        return '[SEG]\t{}\t{}\t{}'.format(self.code_name, self.start, self.end)
    ####
    # get region fasta sequence #
    ####
    def get_fasta(self, ref_fasta, out_fasta=False):
        tmp = 'tmp_fasta'
        if out_fasta: tmp=out_fasta
        bedtool = pybedtools.BedTool(self.to_bed(),from_string=True)
        bedtool.sequence(fi=ref_fasta, fo=tmp, name=True)

        header, seq = open(bedtool.seqfn, 'rb').read().split()
        header, seq = header.strip(), seq.strip()
        if self.strand == '-':
            seq = sequence.reverse_complement(seq)
        if out_fasta:   # return file
            with open(out_fasta, 'wb') as out:
                out.write(header + '\n')
                out.write(seq + '\n')
            return
        else:           # return string
            os.remove(tmp)
            return header, seq

    @staticmethod
    def safe_genomic_region(chrom, start, end):
        if start < 0: start = 0   # start from 0???
        max_end = constants.hg19_fai_bp[chrom]    # hg19 !!!????
        if end > max_end: end = max_end
        return chrom, start, end

if __name__ == "__main__":
    ref_fasta = '/home/BIOINFO_DATABASE/reference/genome_DNA/Homo_sapiens/hg19/BWA_GATK_index/hg19.fa'
    ####
    # region1 and region2 are the same!!! #
    ####
    region1 = GenomicRegion('chr1','100010','100020','+')
    region2 = GenomicRegion('chr1','100010','100020','-')
    '''
    print region1
    header1, seq1 = region1.get_fasta(ref_fasta)
    print seq1, sequence.reverse_complement(seq1)
    print region2
    header2, seq2 = region2.get_fasta(ref_fasta)
    print seq2, sequence.reverse_complement(seq2)
    '''

