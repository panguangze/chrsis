# -*- coding: utf-8 -*-

"""
    utils.complex_haplotype
    ~~~~~~~~~~~~~~~

    Local complex haplotype definition.

    @Copyright: (c) 2017 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""
from genomic_region import GenomicRegion

class Interval(object):

    def __init__(self, pos1, pos2, index, notation='H'):
        self.start = min(pos1, pos2)
        self.end = max(pos1, pos2)
        self.index = index # start from 0
        self.notation = notation
        self.name = self.notation + str(self.index+1)

    def __cmp__(self, other):
        return self.index > other.index

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end and self.name == other.na0
    def __repr__(self):
        return '%s[%s, %s]' % (self.name, str(self.start), str(self.end))
    def __str__(self):
        return str(self.name)
    def to_stat(self):
        return '[SEG]\t{}\t{}\t{}'.format(self.name, self.start, self.end)



class SVLink(object):

    def __init__(self, sv_name=None, up_pos=None, down_pos=None):
        self.sv_name = sv_name # 'DEL', 'INV', 'TD'
        self.up_pos = up_pos
        self.down_pos = down_pos

    def __eq__(self, other):
        return self.sv_name == other.sv_name and self.up_pos == other.up_pos and self.down_pos == other.down_pos

    def __repr__(self):
        if not self.sv_name and not self.up_pos and not self.down_pos: return ' %2s%3s%2s ' % ('','','')
        return ' %2s%3s%2s ' % (str(self.up_pos), self.sv_name, str(self.down_pos))
    def __str__(self):
        return ' %3s ' % (str(self.sv_name) if self.sv_name else '')

    @staticmethod
    def render_sv_link(new_region, last_region, sv_name=None):
        up_pos, down_pos = last_region.get_3prime_pos(), new_region.get_5prime_pos()
        if sv_name: return SVLink(sv_name, up_pos, down_pos)
        # infer sv name
        if new_region.strand != last_region.strand:                             sv_name = 'INV'
        elif new_region.strand == '+':
            if last_region.get_3prime_pos() < new_region.get_5prime_pos()-1:    sv_name = 'DEL'
            if last_region.get_3prime_pos() == new_region.get_5prime_pos()-1:   sv_name = '   '
            if last_region.get_3prime_pos() > new_region.get_5prime_pos()-1:    sv_name = 'TD'
        elif new_region.strand == '-':
            if last_region.get_3prime_pos() > new_region.get_5prime_pos()-1:    sv_name = 'DEL'
            if last_region.get_3prime_pos() == new_region.get_5prime_pos()-1:   sv_name = '   '
            if last_region.get_3prime_pos() < new_region.get_5prime_pos()-1:    sv_name = 'TD'
        return SVLink(sv_name, up_pos, down_pos)


class Contig(object):

    def __init__(self, contig=None):
        self.contig = contig or []

    def append_connection(self, region, sv_link=SVLink()):
        self.contig.append((sv_link, region))

    def __add__(self, other):
        contig = self.contig + other.contig
        return Contig(contig)
    def __len__(self):
        return self.contig.__len__()

    def get_contig_length(self):    return len(self.contig)
    def get_first_region(self):     return self.contig[0][1]
    def get_last_region(self):      return self.contig[-1][1]
    def get_first_sv_link(self):    return self.contig[0][0]
    def get_last_sv_link(self):     return self.contig[-1][0]

    def __getitem__(self, key): # index or slice object
        item = self.contig.__getitem__(key)
        if type(item) == tuple: item = [item]
        return Contig(item)

    def __delitem__(self, key): # index or slice object
        contig = self.__getitem__(key)
        self.contig.__delitem__(key)

    def pop(self, key): # index or slice object
        contig = self.__getitem__(key)
        self.__delitem__(key)
        return contig

    def __repr__(self):
        return ''.join([repr(sv_link) + repr(region) for sv_link, region in self.contig])

    def __str__(self):
        return ''.join([str(sv_link) + str(region) for sv_link, region in self.contig])

    def to_stat(self):
        return ','.join(region.full_name for _, region in self.contig)

    ####
    # transform contig to fasta #
    ####
    def get_fasta(self, ref_fasta, out_fasta=False):
        contig_header, contig_seq = '', ''
        for sv_link, region in self.contig:
            header, seq = region.get_fasta(ref_fasta, out_fasta)
            contig_header += header
            contig_seq += seq
        return contig_header, contig_seq

    @property
    def bp_length(self):
        return sum(region.bp_length for svlink, region in self.contig)

class TandemDuplicationContig(object):

    def __init__(self, repeat_contig, repeat_time=5, generate_type='forward', base_contig=Contig(), margin=0):
        self.repeat_contig = repeat_contig
        self.repeat_time = repeat_time
        self.generate_type = generate_type # 'backword' or 'forward'
        self.base_contig = base_contig # Empty if 'forward'
        self.margin = margin # graphical representation margin, 0 if 'forward'


    def get_last_region(self):
        return self.repeat_contig.get_last_region()

    def __repr__(self):
        base_str = repr(self.base_contig) + '\n' if self.generate_type == 'backward' else ''
        repeat_str = (repr(GenomicRegion())+repr(SVLink()))*self.margin + repr(self.repeat_contig)
        repeat_str_list = [repeat_str]*self.repeat_time
        return base_str + '\n'.join(repeat_str_list)

    def __str__(self):
        base_str = str(self.base_contig) + '\n' if self.generate_type == 'backward' else ''
        repeat_str = (' '*(9+1+3))*self.margin + ' '*self.margin + str(self.repeat_contig)
        repeat_str_list = [repeat_str]*self.repeat_time
        return base_str + '\n'.join(repeat_str_list)

    def to_stat(self):
        if self.generate_type == 'backward':
            res = '{},'.format(self.base_contig.to_stat())
        else:
            res = ''
        res += '({}):{}'.format(self.repeat_contig.to_stat(), self.repeat_time)
        return res

    ####
    # transform tandemduplication contig to fasta #
    ####
    def get_fasta(self, ref_fasta, out_fasta=False):
        base_contig_header, base_contig_seq = self.base_contig.get_fasta(ref_fasta)
        repeat_contig_header, repeat_contig_seq = self.repeat_contig.get_fasta(ref_fasta)
        tdc_header = base_contig_header
        tdc_seq = base_contig_seq + ''.join([repeat_contig_seq]*self.repeat_time)

        return tdc_header, tdc_seq

    @property
    def bp_length(self):
        bp_length = self.repeat_contig.bp_length*self.repeat_time
        if self.generate_type == 'backward':
            bp_length += self.base_contig.bp_length
        return bp_length

class ComplexHaplotype(object):

    def __init__(self, region, interval_list, is_local=True, is_inter=True,
            elongated_length=500000,
            ):
        self.region = region
        self.interval_list = interval_list

        self.contig_list = [Contig()]
        self.is_local = is_local
        self.is_inter = is_inter
        self._in_dup = False
        self.elongated_length = elongated_length

        self.up_elongated_region = self.region.slice(self.region.start-1, self.region.start-elongated_length, '+', 'UP','UP')
        self.down_elongated_region = self.region.slice(self.region.end+1, self.region.end+elongated_length, '+', 'DOWN','DOWN')
        self.elongated_region = self.region.slice(self.region.start-elongated_length, self.region.end+elongated_length, '+', 'elongated', 'elongated')

    ####
    # append methods #
    ####
    def append_contig(self, contig):
        self.contig_list.append(contig)

    def append_connection(self, region, sv_link=SVLink()):
        if type(self.get_last_contig()) != Contig:
            self.contig_list.append(Contig())
        self.get_last_contig().append_connection(region, sv_link)
    ####
    # get methods #
    ####
    def get_last_contig(self): return self.contig_list[-1]
    def get_last_region(self): return self.get_last_contig().get_last_region()
    def pop_subcontig_from_last_contig(self, key):
        if type(self.get_last_contig()) != Contig: return Contig()
        subcontig = self.get_last_contig().pop(key)
        return subcontig

    def __eq__(self, other):
        return self.haplotype == other.haplotype
    def __repr__(self):
        res = ''
        for contig1, contig2 in zip(['']+self.contig_list, self.contig_list+[Contig()]):
            delimiter = '' if type(contig1) == Contig and type(contig2) == TandemDuplicationContig and contig2.generate_type == 'backward' else '\n'
            res += delimiter + repr(contig2)
        return res
    def __str__(self):
        res = self.up_elongated_region.full_name
        for contig1, contig2 in zip(['']+self.contig_list, self.contig_list+[Contig()]):
            delimiter = '' if type(contig1) == Contig and type(contig2) == TandemDuplicationContig and contig2.generate_type == 'backward' else '\n'
            res += delimiter + str(contig2)
        res += self.down_elongated_region.full_name
        return res
    def to_stat(self):
        res = []
        # segment_list
        res.append(self.up_elongated_region.to_stat())
        res.append(self.interval_list.to_stat())
        res.append(self.down_elongated_region.to_stat())
        # length
        res.append('[LENGTH]\t{}\n'.format(self.bp_length))
        # hap
        res.append('[COMPLEX_HAPLOTYPE]\n'+
                self.up_elongated_region.full_name + ',' +
                ','.join(contig.to_stat() for contig in self.contig_list) +
                self.down_elongated_region.full_name
                )
        res.append('[GRAPH]\n'+str(self))
        return '\n'.join(res)

    @property
    def bp_length(self):
        return sum(contig.bp_length for contig in self.contig_list) + self.elongated_length*2


    ####
    # write haplotype to fasta #
    ####
    def write_to_fasta(self, ref_fasta, out_fasta, mod='wb'):

        hap_seq = ''

        header, seq = self.up_elongated_region.get_fasta(ref_fasta)
        hap_seq += seq

        for contig in self.contig_list:
            header, seq = contig.get_fasta(ref_fasta)
            hap_seq += seq
        hap_header = self.elongated_region.to_fasta_header()

        header, seq = self.up_elongated_region.get_fasta(ref_fasta)
        hap_seq += seq


        with open(out_fasta, mod) as out:   out.write(hap_header + '\n' + hap_seq + '\n')
        return hap_header, hap_seq

    def write_to_stat(self, stat_fn, mod='wb'):
        with open(stat_fn, mod) as out:   out.write(self.to_stat()+'\n')

if __name__ == "__main__":
    h = ComplexHaplotype(GenomicRegion('chr1','0','10000000','+','test'))
    '''
    h.append_connection(GenomicRegion('chr1','a','b','+'))
    h.append_connection(GenomicRegion('chr1','c','d','+'), SVLink('DEL','b','c'))

    base_c = Contig()
    base_c.append_connection(GenomicRegion('chr1','e','f','+'), SVLink('DEL','d','e'))
    base_c.append_connection(GenomicRegion('chr1','d','f','-'), SVLink('INV','f','f'))
    base_c.append_connection(GenomicRegion('chr1','b','c','-'), SVLink('DEL','d','c'))
    base_c.append_connection(GenomicRegion('chr1','g','h','+'), SVLink('INV','b','g'))
    repeat_c = Contig()
    repeat_c.append_connection(GenomicRegion('chr1','e','f','+'), SVLink('TD','h','e'))
    repeat_c.append_connection(GenomicRegion('chr1','d','f','-'), SVLink('INV','f','f'))
    repeat_c.append_connection(GenomicRegion('chr1','b','c','-'), SVLink('DEL','d','c'))
    repeat_c.append_connection(GenomicRegion('chr1','g','h','+'), SVLink('INV','b','g'))
    tdc = TandemDuplicationContig(repeat_c,generate_type='backward',base_contig=base_c, margin=2)
    h.append_contig(tdc)

    repeat_c = Contig()
    repeat_c.append_connection(GenomicRegion('chr1','i','j','+'), SVLink('DEL','h','i'))
    repeat_c.append_connection(GenomicRegion('chr1','b','h','+'), SVLink('TD','j','b'))
    tdc=TandemDuplicationContig(repeat_c, repeat_time=3)
    h.append_contig(tdc)
    h.append_connection(GenomicRegion('chr1','h0','k','+'), SVLink('','h','h'))
    h.append_connection(GenomicRegion('chr1','h1','k','+'), SVLink('','h','h'))
    h.append_connection(GenomicRegion('chr1','h2','k','+'), SVLink('','h','h'))
    h.append_connection(GenomicRegion('chr1','h3','k','+'), SVLink('','h','h'))
    '''
    # test for get_fasta
    h.append_connection(GenomicRegion('chr1','100010','100020','+'))
    base_c = Contig()
    base_c.append_connection(GenomicRegion('chr1','100030','100040','+'), SVLink('DEL','100020','100030'))
    base_c.append_connection(GenomicRegion('chr1','100025','100035','+'), SVLink('TD','100040','100025'))
    base_c.append_connection(GenomicRegion('chr1','100035','100025','-'), SVLink('INVS','100035','100035'))
    tdc = TandemDuplicationContig(base_c, repeat_time=2)
    h.append_contig(tdc)


    repeat_c = Contig()
    repeat_c.append_connection(GenomicRegion('chr1','100030','100040','+'), SVLink('INV','100025','100030'))
    repeat_c.append_connection(GenomicRegion('chr1','100025','100035','+'), SVLink('TD','100040','100025'))
    repeat_c.append_connection(GenomicRegion('chr1','100035','100025','-'), SVLink('INVS','100035','100035'))
    tdc = TandemDuplicationContig(repeat_c, generate_type='backward', base_contig=base_c, margin=0, repeat_time=3)
    h.append_contig(tdc)
    #print h
    #print '-'*20

    ref_fasta = '/home/BIOINFO_DATABASE/reference/genome_DNA/Homo_sapiens/hg19/BWA_GATK_index/hg19.fa'
    out_fasta = 'test'
    header, seq = h.get_fasta(ref_fasta, out_fasta)
    #print header
    #print seq
