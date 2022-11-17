# -*- coding: utf-8 -*-
"""
    cnv.cnv
    ~~~~~~~~~~~~~

    copy number variant

    @Copyright: (c) 2018-07 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""
import pandas as pd

from bioutensil import myio
from bioutensil import module


def read_cnv_fn(cnv_fn):  # original csv file
    df = pd.read_csv(cnv_fn, index_col=0)
    return df


class CNVRecord(myio.Record):

    fields = 'id,chr,start,end,np,mean,arm,snvs,ai,median,Cn,mCn,fullCN,meanCn,purity'.split(',')  # tumor_percent

    def __init__(self, *args, **kwargs):
        super(CNVRecord, self).__init__(*args, **kwargs)

    def set(self,
            id=None,
            chr=None,
            start=None,
            end=None,
            np=None,
            mean=None,
            arm=None,
            snvs=None,
            ai=None,
            median=None,
            Cn=None,
            mCn=None,
            fullCN=None,
            meanCn=None,
            purity=None,
            args=None,
            **kwargs):

        self.id = self._validate(id, int)
        self.chr = self._validate(chr)
        self.start = self._validate(start, int)
        self.end = self._validate(end, int)
        self.np = self._validate(np, int)
        self.mean = self._validate(mean, float)
        self.arm = self._validate(arm)
        self.snvs = self._validate(snvs, int)
        self.ai = self._validate(ai, float)
        self.median = self._validate(median, float)
        self.Cn = self._validate(Cn, int)
        self.mCn = self._validate(mCn, int)
        self.fullCN = self._validate(fullCN)
        self.meanCn = self._validate(meanCn, float)
        self.purity = self._validate(purity, float)

    def parse(self, line):  # parse a stream_or_string to cnv object
        args = dict(zip(self.fields, line.strip().split(',')))
        self.set(**args)

    def __str__(self):
        return ','.join(str(getattr(self, field)) for field in self.fields)


class CNVReader(myio.Reader):

    def __init__(self, cnv_cls=CNVRecord, sep='\t', has_header=True, *args, **kwargs):
        super(CNVReader, self).__init__(record_cls=cnv_cls, sep=sep, has_header=has_header, *args, **kwargs)


class CNVWriter(myio.Writer):

    def __init__(self, cnv_cls=CNVRecord, sep='\t', *args, **kwargs):
        super(CNVWriter, self).__init__(record_cls=cnv_cls, sep=sep, *args, **kwargs)


class CNVGenerator(module.Module):

    def __init__(self, reader_cls=CNVReader, writer_cls=CNVWriter, *args, **kwargs):
        super(CNVGenerator, self).__init__(reader_cls=reader_cls, writer_cls=writer_cls,
                                           *args, **kwargs)
    def __iter__(self):
        for i, record in enumerate(self.in_iter):
            print(record)
            meta, record = record
            yield meta, record


# if __name__ == "__main__":
#     fn = '/home/chenlingxi/mnt/chenlingxi/workspace/FS_Projects/WHOC.WGS_10X.batch01/CNACalling/patchwork/OC001/OC001/oc001t_Copynumbers.csv'
#     CNVGenerator(in_data=fn).evaluate()
