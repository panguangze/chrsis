# -*- coding: utf-8 -*-
"""
    sv.vcf
    ~~~~~~

    @Copyright: (c) 2018-03 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""

import gzip
import vcf
import copy

import structure_variant
from tenxtools.utils import basic


class Record(vcf.model._Record, structure_variant.SVRecord):

    def __init__(self, *args, **kwargs):
        vcf.model._Record.__init__(self, *args, **kwargs)
        structure_variant.SVRecord.__init__(self)

    @property
    def to_sv(self):
        if self.var_type == 'INV':
            return '\n'.join(str(x) for x in self.bnd_records)
        else:
            return structure_variant.SVRecord.__str__(self)

    @property
    def var_type(self):
        if self.is_snp:
            return "snp"
        elif self.is_indel:
            return "indel"
        elif self.is_sv:
            return self.INFO.get('SVTYPE', 'sv')
        else:
            return "unknown"

    @property
    def secondary(self):
        return self.INFO.get('SECONDARY', False)

    @property
    def sv_end(self):
        return self.INFO.get('END', None)

    def set_sv_fields(self):
        pass

    @property
    def _get_general_sv_fields(self):
        args = {}
        args['id'] = self.ID
        args['qual'] = self.QUAL
        args['filter'] = self.FILTER
        args['span_reads'] = self.INFO.get('PE', [None])[0]
        args['junc_reads'] = self.INFO.get('SR', [None])[0]
        return args

    @property
    def _get_inv_bnd_records(self):  # transform INV format to BNDs format
        '''
        pos:            ab            cd
                        AC            GT
        +: 5' ----------oo------------oo---------- 3'
        -: 3' ----------oo------------oo---------- 5'

        INV Format:
        #CHROM  POS ID      REF ALT         QUAL    FILTER  INFO
        chr1    b   1       N   <INV>       10      PASS    SVTYPE=INV;END=c

        BNDs Format:
        - BDD1
        #CHROM  POS ID      REF ALT         QUAL    FILTER  INFO
        chr1    a   1.1_1   N   A]chr1:c]   10      PASS    SVTYPE=BND;MATEID=.;EVENT=1
        chr1    c   1.1_2   N   G]chr1:a]   10      PASS    SVTYPE=BND;MATEID=.;EVENT=1;SECONDARY

                        a
                        A
        +: 5' ----------o 3' ==================>
        -: 3' ----------o 5'                   |
                                               |
                                      c        |
                                      G        |
        +:            5' o------------o 3'     |
        -:            3' o------------o 5' <===

        - BND2
        chr1    b   1.2_1   N   [chr1:d[C   10      PASS    SVTYPE=BND;MATEID=.;EVENT=1
        chr1    d   1.2_2   N   [chr1:b[T   10      PASS    SVTYPE=BND;MATEID=.;EVENT=1:SECONDARY

                         b
                         C
        +:            5' o------------o 3'
        -:   <======= 3' o------------o 5'
             |
             |                         d
             |                         T
        +:   =====================> 5' o---------- 3'
        -:                          3' o---------- 5'

        '''
        args = {}
        args['chrom_5p'] = self.CHROM
        args['chrom_3p'] = self.CHROM
        args.update(self._get_general_sv_fields)
        args['meta_info'] = self._get_sv_meta_info   # no need to change for INV
        args['anno_info'] = self._get_sv_anno_info

        id = args['id']
        # BND 1
        args['id'] = id + '.1_1'
        args['bkpos_5p'] = self.POS - 1
        args['bkpos_3p'] = self.sv_end
        args['strand_5p'] = '+'     # ++ version
        args['strand_3p'] = '-'
        record1 = structure_variant.SVRecord(args=args)

        # BND2
        args['id'] = id + '.2_1'
        args['bkpos_5p'] = self.POS
        args['bkpos_3p'] = self.sv_end + 1
        args['strand_5p'] = '-'     # -- verson
        args['strand_3p'] = '+'
        record2 = structure_variant.SVRecord(args=args)
        return record1, record2

    @property
    def _get_inv_sv_fields(self):  # just for INV type comparision
        args = {}
        args['chrom_5p'] = self.CHROM
        args['bkpos_5p'] = self.POS
        args['strand_5p'] = '+'

        args['chrom_3p'] = self.CHROM
        args['bkpos_3p'] = self.sv_end
        args['strand_3p'] = '+'
        return args

    @property
    def _get_del_sv_fields(self):
        args = {}
        args['chrom_5p'] = self.CHROM
        args['bkpos_5p'] = self.POS
        args['strand_5p'] = '+'

        args['chrom_3p'] = self.CHROM
        args['bkpos_3p'] = self.sv_end
        args['strand_3p'] = '+'

        args['inner_ins'] = None
        return args

    @property
    def _get_dup_sv_fields(self):   # +: larger, small, switch start and end
        args = {}
        args['chrom_5p'] = self.CHROM
        args['bkpos_5p'] = self.sv_end
        args['strand_5p'] = '+'

        args['chrom_3p'] = self.CHROM
        args['bkpos_3p'] = self.POS
        args['strand_3p'] = '+'

        args['inner_ins'] = None
        return args

    @property
    def _get_unk_sv_fields(self):
        args = {}
        args['chrom_5p'] = self.CHROM
        args['bkpos_5p'] = self.POS
        args['strand_5p'] = '+'

        args['chrom_3p'] = self.CHROM
        args['bkpos_3p'] = self.sv_end
        args['strand_3p'] = '+'

        args['inner_ins'] = None
        return args

    @property
    def _get_bnd_sv_fields(self):
        '''
        VCF Format 1
        #CHROM  POS ID      REF ALT         QUAL    FILTER  INFO
        chr1    a   1_1     N   A]chr2:b]   10      PASS    SVTYPE=BND;MATEID=.;EVENT=1
        chr2    b   1_2     N   [chr1:a[C   10      PASS    SVTYPE=BND;MATEID=.;EVENT=1;SECONDARY

        SV Format
        #chrom_5p   bkpos_5p    strand_5p   chrom_3p    bkpos_3p    strand_3p
        chr1        a           +           chr2        b           +
                        a
                        A
        chr1 +: 5' ----------o 3' =========
        chr1 -: 3' ----------o 5'          |
                                           |
                                           |             b
                                           |             C
        chr2 +:                             ========> 5' o------------o 3'
        chr2 -:                                       3' o------------o 5'

        VCF Format 2
        #CHROM  POS ID      REF ALT         QUAL    FILTER  INFO
        chr1    a   1_1     N   A[chr2:b[   10      PASS    SVTYPE=BND;MATEID=.;EVENT=1
        chr2    b   1_2     N   ]chr1:a]C   10      PASS    SVTYPE=BND;MATEID=.;EVENT=1;SECONDARY

        SV Format
        #chrom_5p   bkpos_5p    strand_5p   chrom_3p    bkpos_3p    strand_3p
        chr1        a           +           chr2        b           +
                                                         b
                                                         C
        chr2 +:                             ========> 5' o------------o 3'
        chr2 -:                            |          3' o------------o 5'
                        a                  |
                        A                  |
        chr1 +: 5' ----------o 3' =========
        chr1 -: 3' ----------o 5'

        '''
        args = {}
        alt = self.ALT[0]

        # if strand1 == '-' and strand2 == '-':   # inversion
        if alt.orientation and alt.remoteOrientation:   # inversion
            args['strand_5p'] = '-'
            args['strand_3p'] = '+'
        # elif strand1 == '+' and strand2 == '+':  # inversion
        elif not alt.orientation and not alt.remoteOrientation:  # inversion
            args['strand_5p'] = '+'
            args['strand_3p'] = '-'
        else:                                   # not inversion
            args['strand_5p'] = '+'
            args['strand_3p'] = '+'

        if alt.orientation:  # If the breakend is connected to sequence 3', extending to 5'
            args['chrom_5p'] = self.CHROM
            args['bkpos_5p'] = self.POS
            args['chrom_3p'] = alt.chr
            args['bkpos_3p'] = alt.pos

            args['inner_ins'] = alt.connectingSequence[:-1]

        else:               # If the breakend is connected to sequence 5', extending to 3'
            args['chrom_5p'] = alt.chr
            args['bkpos_5p'] = alt.pos
            args['chrom_3p'] = self.CHROM
            args['bkpos_3p'] = self.POS

            args['inner_ins'] = alt.connectingSequence[1:]

        return args

    @property
    def _get_sv_meta_info(self):
        return self.INFO

    @property
    def _get_sv_anno_info(self):
        self.anno_info = None


class LumpyRecord(Record):

    def __init__(self, *args, **kwargs):
        super(LumpyRecord, self).__init__(*args, **kwargs)

    def __str__(self):
        res = "Record(CHROM=%(CHROM)s, POS=%(POS)s, REF=%(REF)s, ALT=%(ALT)s), " % self.__dict__
        res += 'END={}, TYPE={}, STRANDS={})'.format(self.sv_end, self.var_type, self.strands)
        return res

    ''' newly add '''
    @property
    def strands(self):
        return self.INFO.get('STRANDS', None)

    def set_sv_fields(self):
        args = {}
        if self.var_type == 'INV':
            self.bnd_records = self._get_inv_bnd_records
            args = self._get_inv_sv_fields
            self.set(**args)
            return

        if self.var_type == 'DEL':
            args = self._get_del_sv_fields
        if self.var_type == 'DUP':
            args = self._get_dup_sv_fields
        if self.var_type == 'BND':
            args = self._get_bnd_sv_fields
        if self.var_type == 'DUP:TANDEM':
            args = self._get_dup_sv_fields
        if self.var_type == 'CNV':
            args = None         # ???

        args.update(self._get_general_sv_fields)

        args['meta_info'] = self._get_sv_meta_info
        args['anno_info'] = self._get_sv_anno_info

        self.set(**args)

    @property
    def _get_sv_meta_info(self):
        meta_info = copy.deepcopy(self.INFO)

        # delete tags
        del_tags = 'END,STRANDS,PE,SR'.split(',')
        for tag in del_tags:
            if tag in meta_info:
                del meta_info[tag]

        # switch tags
        switch_tags = [
            ('CIPOS', 'CIEND'),
            ('CIPOS95', 'CIEND95'),
            ('PRPOS', 'PREND')
        ]
        # longranger, only record CIPOS
        # If the breakend is connected to sequence 5', extending to 3', switch, lumpy didn't need to switch???
        if (self.var_type in ['DUP', 'DUP:TANDEM']) or (self.var_type == 'BND' and not self.ALT[0].orientation):
            for l_tag, r_tag in switch_tags:
                if l_tag in meta_info and r_tag in meta_info:
                    temp = meta_info[l_tag]
                    meta_info[l_tag] = meta_info[r_tag]
                    meta_info[r_tag] = temp

        # add new tags
        meta_info['GT'] = ','.join([sample['GT'] for sample in self.samples])
        meta_info['TOOL'] = 'Lumpy'

        # update mateid for INV
        if self.var_type == 'INV':
            meta_info['MATEID'] = None

        return meta_info


class LongrangerRecord(Record):

    def __init__(self, *args, **kwargs):
        super(LongrangerRecord, self).__init__(*args, **kwargs)

    def __str__(self):
        res = "Record(CHROM=%(CHROM)s, POS=%(POS)s, REF=%(REF)s, ALT=%(ALT)s), " % self.__dict__
        res += 'END={}, TYPE={})'.format(self.sv_end, self.var_type)
        return res

    @property
    def sv_type(self):
        if self.is_snp:
            return "snp"
        elif self.is_indel:
            return "indel"
        elif self.is_sv:
            return self.INFO.get('SVTYPE', 'sv')
        else:
            return "unknown"

    @property
    def sv_type2(self):
        return self.INFO.get('SVTYPE2', None)

    @property
    def var_type(self):
        if self.sv_type == 'BND' and not self.sv_type2:
            return self.sv_type2
        else:
            return self.sv_type

    def set_sv_fields(self):
        args = {}
        if self.var_type == 'INV':
            self.bnd_records = self._get_inv_bnd_records
            args = self._get_inv_sv_fields
            self.set(**args)
            return

        if self.var_type == 'DEL':
            args = self._get_del_sv_fields
        if self.var_type == 'DUP':
            args = self._get_dup_sv_fields
        if self.var_type == 'BND':
            args = self._get_bnd_sv_fields
        if self.var_type == 'DUP:TANDEM':
            args = self._get_dup_sv_fields

        if self.var_type == 'UNK':
            args = self._get_unk_sv_fields  # ???

        args.update(self._get_general_sv_fields)

        args['meta_info'] = self._get_sv_meta_info
        args['anno_info'] = self._get_sv_anno_info

        self.set(**args)

    @property
    def _get_sv_meta_info(self):
        meta_info = copy.deepcopy(self.INFO)

        # switch tags
        switch_tags = [
            ('CIPOS', 'CIEND'),
        ]
        if self.var_type in ['DUP', 'DUP:TANDEM']:
            for l_tag, r_tag in switch_tags:
                temp = meta_info[l_tag]
                meta_info[l_tag] = meta_info[r_tag]
                meta_info[r_tag] = temp

        # add new tags
        meta_info['GT'] = ','.join([sample['GT'] for sample in self.samples])
        meta_info['TOOL'] = 'Longranger'

        return meta_info


class Reader(vcf.parser.Reader):

    def __init__(self, record_cls=Record, *args, **kwargs):
        self._record_cls = record_cls
        super(Reader, self).__init__(*args, **kwargs)

    def next(self):
        '''Return the next record in the file.'''
        line = next(self.reader)
        row = self._row_pattern.split(line.rstrip())
        chrom = row[0]
        if self._prepend_chr:
            chrom = 'chr' + chrom
        pos = int(row[1])

        if row[2] != '.':
            ID = row[2]
        else:
            ID = None

        ref = row[3]
        alt = self._map(self._parse_alt, row[4].split(','))

        try:
            qual = int(row[5])
        except ValueError:
            try:
                qual = float(row[5])
            except ValueError:
                qual = None

        filt = self._parse_filter(row[6])
        info = self._parse_info(row[7])

        try:
            fmt = row[8]
        except IndexError:
            fmt = None
        else:
            if fmt == '.':
                fmt = None

        record = self._record_cls(chrom, pos, ID, ref, alt, qual, filt, info, fmt, self._sample_indexes)  # change the _Record type to Record

        if fmt is not None:
            samples = self._parse_samples(row[9:], fmt, record)
            record.samples = samples

        return record

    def _parse_filter(self, filt_str):
        '''Parse the FILTER field of a VCF entry into a Python list

        NOTE: this method has a cython equivalent and care must be taken
        to keep the two methods equivalent
        '''
        if filt_str == '.':
            return None
        elif filt_str == 'PASS':
            return []
        else:
            return filt_str.split(';')


class VCF(object):

    def __init__(self, vcf_fn, allowance=10, record_cls=Record):
        self._vcf_fn = vcf_fn
        self._allowance = allowance
        self._record_cls = record_cls

    @property
    def reader(self):
        # return Reader(basic.safe_open(self._vcf_fn, 'r'))
        compressed = self._vcf_fn.endswith('gz')
        return Reader(filename=self._vcf_fn, compressed=compressed, record_cls=self._record_cls)

    def read(self, ex_filter_funcs=[], in_filter_funcs=[]):
        for record in self.reader:
            if (not any(filter_func(record) for filter_func in ex_filter_funcs)) \
                    and (all(filter_func(record) for filter_func in in_filter_funcs)):
                yield record

    @property
    def to_sv(self):  # vcf to sv string
        for record in self.reader:
            record.set_sv_fields()
            yield 'whole', record.to_sv

    @property
    def set_sv_fields(self):  # vcf to sv record
        for record in self.reader:
            record.set_sv_fields()
            yield 'whole', record

    def write(self, out_fn, records, fmt='vcf'):
        if fmt == 'vcf':
            writer = vcf.Writer(open(out_fn, 'wb'), self.reader)
        if fmt == 'sv':
            writer = structure_variant.Writer(out_fn)
        for record in records:
            writer.write_record(record)

    def write_inter_diff(self, shared_fn, diff_fn, shared_or_diff, fmt='vcf'):
        if fmt == 'vcf':
            inter_writer = vcf.Writer(open(shared_fn, 'wb'), self.reader)
            diff_writer = vcf.Writer(open(diff_fn, 'wb'), self.reader)
        if fmt == 'sv':
            inter_writer = structure_variant.Writer(shared_fn)
            diff_writer = structure_variant.Writer(diff_fn)
        for tag, record in shared_or_diff:
            if tag == 'intersection':
                inter_writer.write_records([record])
            else:
                diff_writer.write_records([record])

    def fetch(self, chrom, start, end, fmt='vcf', ex_filter_funcs=[], in_filter_funcs=[]):  # check
        try:
            records = self.reader.fetch(chrom, start, end)
        except:
            return

        for record in records:
            if (not any(filter_func(record) for filter_func in ex_filter_funcs)) \
                    and (all(filter_func(record) for filter_func in in_filter_funcs)):

                if fmt == 'sv':
                    record.set_sv_fields()
                yield record

    def inter_diff(self, other, ex_filter_funcs=[], in_filter_funcs=[]):
        for i, self_record in enumerate(self.read(ex_filter_funcs=ex_filter_funcs, in_filter_funcs=in_filter_funcs)):
            self_record.set_sv_fields()
            chrom = self_record.chrom_5p
            pos = self_record.bkpos_5p
            start, end = pos - self._allowance, pos + self._allowance

            other_records = other.fetch(chrom, start, end, fmt='sv', ex_filter_funcs=ex_filter_funcs, in_filter_funcs=in_filter_funcs)
            if self_record in other_records:  # variants in self that also in other, self record is written
                yield 'intersection', self_record
            else:  # variants in self but not in other
                yield 'difference', self_record

    def intersection(self, other):
        for tag, record in self.inter_diff(other):
            if tag == 'intersection':
                yield self_record

    def difference(self, other):
        for tag, record in self.inter_diff(other):
            if tag == 'difference':
                yield self_record
