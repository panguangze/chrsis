# -*- coding: utf-8 -*-
"""
    tenxtools.utils.psl_tran
    ~~~~~~~~~~~~~~

    PSL transcript class

    @Copyright: (c) 2018-04 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""
import pysam
from tenxtools.utils import genomic_region
from tenxtools.utils import gene


class Reader(object):

    def __init__(self, fn, gene_list=None):
        self._fn = fn
        self.gene_list = gene.read_gene_list_from_dir(gene_list)

    def fetch(self, chrom, pos1, pos2):
        if pos1 > pos2:
            start, end = pos2, pos1
        else:
            start, end = pos1, pos2
        tbx = pysam.TabixFile(self._fn)
        try:
            for i, row in enumerate(tbx.fetch(chrom, start, end)):
                yield Transcript(row)
        except Exception:
            return

    def fetch_genes_trans(self, chrom, pos1, pos2):
        groups = []
        gene = None
        for i, tran in enumerate(self.fetch(chrom, pos1, pos2)):
            if gene != tran.gene:
                if groups:
                    yield groups
                groups = []
            groups.append(tran)
            gene = tran.gene
        # yield last group
        if groups:
            yield groups

    def fetch_genes_tran_by_point(self, chrom, point):
        return self.fetch_genes_tran_by_region(chrom, point, point+1)

    def fetch_genes_tran_by_region(self, chrom, pos1, pos2):
        genes_trans = self.fetch_genes_trans(chrom, pos1, pos2)
        for i, gene_trans in enumerate(genes_trans):
            tran = self.choose_tran_from_group(gene_trans)
            if not tran:
                continue
            for group, genes in self.gene_list.items():
                if tran.gene in genes:
                    tran.group.append(group)
            yield tran

    def fetch_genes_tran_by_seg(self, chrom, pos, strand, seg_extend=500000):
        if strand == '+':
            start, end = pos - seg_extend, pos
        else:
            start, end = pos, pos + seg_extend
        chrom, start, end = genomic_region.GenomicRegion.safe_genomic_region(chrom, start, end)
        return self.fetch_genes_tran_by_region(chrom, start, end)

    @staticmethod
    def choose_tran_from_group(trans,
                               default_biotypes=['protein_coding', 'lincRNA'],  # ,'miRNA'],
                               other_biotypes=[]):
        # choose one transcript to represent the gene transcript
        accepted_biotypes = default_biotypes + other_biotypes

        protein_001, protein_201, protein_min = None, None, None
        other_001, other_201, other_min = None, None, None
        protein_name_id_min, other_name_id_min = '999', '999'
        for t in trans:
            if t.biotype not in accepted_biotypes:
                continue
            name_id = t.name.split('-')[-1]
            if 'protein_encoding' in t.biotype:
                if name_id == '001':
                    protein_001 = t
                if name_id == '201':
                    protein_201 = t
                if name_id < protein_name_id_min:
                    protein_min = t
                    protein_name_id_min = name_id
            else:
                if name_id == '001':
                    other_001 = t
                if name_id == '201':
                    other_201 = t
                if name_id < other_name_id_min:
                    other_min = t
                    other_name_id_min = name_id

        if protein_001:
            return protein_001
        if protein_201:
            return protein_201
        if protein_min:
            return protein_min
        if other_001:
            return other_001
        if other_201:
            return other_201
        if other_min:
            return other_min
        return None


class Transcript(object):

    def __init__(self, line):

        self._fields = 'source,version,col3,col4,start_condon,end_condon,protein_ens_id,col8,strand,name,ens_id,gtf_name,biotype,chrom,cytoband,left,right,exon_num,exon_lens,cds_info,exon_starts,gene'.split(',')

        self._parse(line)

    def _set(self, source=None, version=None, col3=None, col4=None,
             start_condon=None, end_condon=None, protein_ens_id=None, col8=None,
             strand=None, name=None, ens_id=None, gtf_name=None,
             biotype=None, chrom=None, cytoband=None, left=None,
             right=None, exon_num=None, exon_lens=None, cds_info=None,
             exon_starts=None, gene=None,
             extend=20000,
             group=None):
        self.source = source
        self.version = version
        self.start_condon = tuple(int(x) for x in start_condon.split(',') if x)
        self.end_condon = tuple(int(x) for x in end_condon.split(',') if x)
        self.protein_ens_id = protein_ens_id
        self.col8 = col8
        self.strand = strand
        self.name = name
        self.ens_id = ens_id
        self.gtf_name = gtf_name
        self.biotype = biotype
        self.chrom = chrom
        self.cytoband = cytoband
        self.left = int(left) + 1   # first exon start
        self.right = int(right)     # last exon end
        self.exon_num = exon_num
        self.exon_lens = list(int(x) for x in exon_lens.split(',') if x)
        self.cds_info = cds_info
        self.exon_starts = list(int(x) for x in exon_starts.split(',') if x)
        self.gene = gene

        self.extend = int(extend)
        self.left_extend = self.left - self.extend
        self.right_extend = self.right + self.extend

        self.group = group or []

    def _parse(self, line):  # parse line to a record object
        args = dict(zip(self._fields, line.strip().split()))
        self._set(**args)

    def __str__(self):
        return '{} {} {} {} {} {} {}'.format(self.chrom, self.left, self.right, self.name, self.biotype, self.gene, self.strand)

    @staticmethod
    def _in_region(point, region):  # left and right are all in region
        if not region:
            return False
        left, right = region
        return left <= point and right >= point

    def in_exons(self, point):
        for i, exon in enumerate(self.exons):
            if self._in_region(point, exon):
                return i+1
        return 0

    def in_introns(self, point):
        for i, intron in enumerate(self.introns):
            if self._in_region(point, intron):
                return i+1
        return 0

    def in_utr5(self, point):
        return self._in_region(point, self.utr5)

    def in_utr3(self, point):
        return self._in_region(point, self.utr3)

    def in_promoter(self, point):
        return self._in_region(point, self.promoter)

    def in_terminator(self, point):
        return self._in_region(point, self.terminator)

    def in_upstream(self, point):  # 5' upstream
        if self.strand == '+':
            return point < self.left_extend
        if self.strand == '-':
            return point > self.right_extend

    def in_downstream(self, point):  # 3' downstream
        if self.strand == '+':
            return point > self.right_extend
        if self.strand == '-':
            return point < self.left_extend

    def in_genetic(self, point):
        return self.in_exons(point) \
                or self.in_introns(point) \
                or self.in_utr5(point) \
                or self.in_utr3(point)

    def mark_functional_region(self, point):
        g = 'GENETIC'
        ig = 'INTERGENETIC'
        reg = 'REGULATION'

        if self.in_utr5(point):
            return g, "utr5"
        if self.in_utr3(point):
            return g, "utr3"
        # utr5 and utr3 are inside exon
        exon = self.in_exons(point)
        if exon:
            return g, 'e{}'.format(exon)
        intron = self.in_introns(point)
        if intron:
            return g, 'i{}'.format(intron)
        if self.in_promoter(point):
            return reg, "promoter"
        if self.in_terminator(point):
            return reg, "terminator"
        if self.in_upstream(point):
            return ig, "upstream"
        if self.in_downstream(point):
            return ig, "downstream"
        return None, None

    def mark_exon_distance(self, point):
        if self.in_exons(point):
            return 0
        intron_i = self.in_introns(point)
        if intron_i:
            l, r = self.get_intron(intron_i)
            l, r = l-1, r+1
            return min(point-l, r-point)
            ''' in intergenetic or regulation area '''
        else:
            return min(abs(point-self.left), abs(point-self.right))
        return None

    def get_exon(self, i):
        exons = self.exons
        if abs(i) > len(exons) or abs(i) < 1:
            return None, None
        if i > 0:
            return exons[i-1]
        if i < 0:
            return exons[i]

    def get_intron(self, i):
        introns = self.introns
        if abs(i) > len(introns) or abs(i) < 1:
            return None, None
        if i > 0:
            return introns[i-1]
        if i < 0:
            return introns[i]

    @property
    def exons(self):  # start = start + 1, end = start + len, len is the number of base pairs
        exon_starts = [start + 1 for start in self.exon_starts]
        exon_ends = [start + l for start, l in zip(self.exon_starts, self.exon_lens)]
        if self.strand == '-':
            exon_starts = exon_starts[::-1]
            exon_ends = exon_ends[::-1]
        return list(zip(exon_starts, exon_ends))

    @property
    def introns(self):
        if self.strand == '+':
            return [(left[-1]+1, right[0]-1) for left, right in zip(self.exons[:-1], self.exons[1:])]
        if self.strand == '-':
            return [(left[-1]+1, right[0]-1) for right, left in zip(self.exons[:-1], self.exons[1:])]

    @property
    def utr5(self):
        if not self.start_condon:
            return tuple()
        if self.strand == '+':
            return (self.left, self.start_condon[0] - 1)
        if self.strand == '-':
            return (self.start_condon[-1] + 1, self.right)

    @property
    def utr3(self):
        if not self.end_condon:
            return tuple()
        if self.strand == '+':
            return (self.end_condon[-1] + 1, self.right)
        if self.strand == '-':
            return (self.left, self.end_condon[0] - 1)

    @property
    def promoter(self):
        if self.strand == '+':
            return (self.left_extend, self.left - 1)
        if self.strand == '-':
            return (self.right + 1, self.right_extend)

    @property
    def terminator(self):
        if self.strand == '-':
            return (self.left_extend, self.left - 1)
        if self.strand == '+':
            return (self.right + 1, self.right_extend)


'''
    def test(self, start=100, end=200):
        print 'left', self.left
        print 'right', self.right
        print self.right - self.left
        print '--- exons'
        print 'exon_num', self.exon_num
        print 'exon_starts', self.exon_starts
        print 'exon_lens', self.exon_lens
        print 'exons', self.exons
        print '--- intron'
        print self.introns

        print '--- utr'
        print 'start condon', self.start_condon
        print 'end condon', self.end_condon
        print 'cds_info', self.cds_info   #????
        print '3utr', self.utr3
        print '5utr', self.utr5
        print 'promoter', self.promoter
        print 'terminator', self.terminator
        if self.start_condon and self.end_condon:
            print 'start condon - end condon', self.start_condon[0] - self.end_condon[-1]
        print '-----distance----'
        point = self.left - 1000
        print 'point', point, 'distance', self.mark_exon_distance(point)
        point = self.right + 1000
        print 'point', point, 'distance', self.mark_exon_distance(point)
        point = 1167272
        print 'point', point, 'distance', self.mark_exon_distance(point)
        point = 1159349 + 99
        print 'point', point, 'distance', self.mark_exon_distance(point)
        point = 1163847 - 99
        print 'point', point, 'distance', self.mark_exon_distance(point)

if __name__ == "__main__":
    line = 'ensembl NA      0       0       1164171,1164172,1164173,        1156079,1156080,1156081,        ENSP00000444451 0       -       SDF4-201        ENST00000545427 SDF4-201        protein_coding  chr1    p36.33  1156078 1167370 5       32,114,137,479,99,   1156081(29),1158623(114),1159211(137),1163847(326),     1156078,1158623,1159211,1163847,1167271,        SDF4'
    print line
    tran = Transcript(line)
    tran.test()

'''
