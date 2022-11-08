# -*- coding: utf-8 -*-
"""
    tenxtools.io
    ~~~~~~~~~~~~


    @Copyright: (c) 2018-06 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""
import os
import gzip
import csv
import vcf
import yaml
import pysam


class Record(object):

    fields = []
    sep = '\t'

    def __init__(self, stream_or_str=None, parent=None, args=None, sep='\t',
                 **kwargs):
        self.parent = parent
        self.sep = sep
        args = args or {}

        if parent:
            self.from_parent()
        elif stream_or_str:
            self.parse(stream_or_str)
        elif args:
            self.set(**args)
        else:
            self.set(**kwargs)

    def set(self, **kwargs):

        for key in kwargs:
            setattr(self, key, kwargs[key])

        for k in kwargs.keys():
            if k not in self.fields:
                self.fields.append(k)

    def from_parent(self, **kwargs):
        pass

    def parse(self, stream_or_str):
        # parse a stream_or_string to record object
        pass

    def _validate(self, x, func=None):
        if func in [None, int, float]:
            if not x:
                return x
            if x == '.':
                return None
        if func is str:
            if not x:
                return ''
        if func:
            return func(x)
        return x

    def _parse_list(self, line, sep=','):
        if isinstance(line, list):
            return line
        if not line or line == '.':
            return []
        return line.split(sep)

    def _format_list(self, x, sep=','):
        return sep.join(map(self._format_value, x))

    def _format_value(self, x):  # None, int, float, list
        if x == 0:
            return '0'
        if not x:
            return '.'   # notice not 0 is True
        if type(x) == list:
            return self._format_list(x)
        return str(x)

    def _parse_value(self, x):  # None, int, float, list
        if x == '.':
            return None
        splits = x.split(',')
        if splits:
            return map(self._parse_value, splits)
        return float(x)

    def __repr__(self):
        res_list = []
        for field in self.fields:
            try:
                x = self._format_value(getattr(self, field))
                res_list.append(x)
            except Exception:
                pass
        return self.sep.join(res_list)

    def __str__(self):
        return self.__repr__()

    def __eq__(self):
        pass


def safe_open(fn, mode):
    if not fn:
        return None
    if fn.endswith('gz'):
        return gzip.open(fn, mode)
    else:
        return open(fn, mode)


class Writer(object):

    def __init__(self, fn=None, record_cls=Record,
                 fields_prefix='#',
                 *args, **kwargs):
        self.fn = fn

        self.fields_prefix = fields_prefix
        self.sep = record_cls.sep
        self.fields = record_cls.fields
        self.write('', 'w')

    def write_headers(self, mode='a'):
        pass

    def write_fields(self, fields=None, mode='a'):
        if fields is None:
            fields = self.fields

        line = self.fields_prefix + self.sep.join(fields) + '\n'
        with safe_open(self.fn, mode) as f:
            f.write(line)
        return True

    def write_chunk(self, chunk, mode='a', dynamic_fields=False):
        with safe_open(self.fn, mode) as f:
            for i, record in enumerate(chunk):
                if i == 0 and dynamic_fields:
                    self.write_fields(fields=record.fields)
                f.write(str(record) + '\n')

    def write(self, content, mode='a'):
        with safe_open(self.fn, mode) as f:
            f.write(content)


class Reader(object):  # currently just csv, tsv file

    record_cls = Record

    def __init__(self,
                 in_data=None, in_fn=None,
                 record_cls=Record,
                 file_type='csv', sep=',', has_header=False,
                 sample=None):
        self.in_data = in_data
        self.in_fn = in_fn
        self.record_cls = record_cls
        self.fields = record_cls.fields
        self.file_type = file_type
        self.sep = sep
        self.has_header = has_header

        self.sample = sample
        if isinstance(self.in_fn, str) and self.in_fn.endswith('gz'):
            self.tbx = pysam.TabixFile(self.in_fn)

        if isinstance(self.in_fn, list) or isinstance(self.in_data, list):
            self.in_iter = self._multiple_data(in_datas=in_data, in_fns=in_fn)
        else:
            self.in_iter = self._single_data(in_data=in_data, in_fn=in_fn)

    def _read(self, in_fn):

        for line in safe_open(in_fn, 'rt'):
            if self.file_type != 'vcf' and line.startswith('#'):
                continue

            yield line

    def _single_data(self, in_data=None, in_fn=None):
        if in_data:
            return in_data
        if in_fn:
            if self.file_type == 'csv':
                in_iter = self._csv_fn(in_fn)
            elif self.file_type == 'vcf':
                in_iter = self._vcf_fn(in_fn)
            elif self.file_type == 'yaml':
                in_iter = self._yaml_fn(in_fn)
            else:
                in_iter = None
        else:
            in_iter = None
        return in_iter

    def _csv_fn(self, in_data):
        if self.has_header:
            in_iter = csv.DictReader(
                self._read(in_data),
                delimiter=self.sep, skipinitialspace=True)
        else:
            in_iter = csv.DictReader(
                self._read(in_data),
                fieldnames=self.fields,
                delimiter=self.sep, skipinitialspace=True)
        return in_iter

    def _vcf_fn(self, in_data):
        in_iter = vcf.Reader(self._read(in_data))
        return in_iter

    def _yaml_fn(self, in_fn):
        docs = yaml.load_all(safe_open(in_fn, 'r'))
        for doc in docs:
            if not doc:
                continue
            yield doc

    def _multiple_data(self, in_datas=None, in_fns=None):
        if in_datas:
            for in_data in in_datas:
                for x in self._single_data(in_data=in_data):
                    yield x
        elif in_fns:
            for in_fn in in_fns:
                for x in self._single_data(in_fn=in_fn):
                    yield os.path.split(in_fn)[1], x

    def __iter__(self):
        return self.read_record()

    def conditional_read(self, condition_func=None):
        for meta, record in self:
            if condition_func and condition_func(record):
                yield meta, record

    def fetch(self, chrom, start, end, condition_func=None):
        for line in self.tbx.fetch(chrom, start, end):
            row = dict(zip(self.record_cls.fields, line.split(self.sep)))
            record = self.record_cls(args=row)
            if condition_func is None:
                yield None, record
            elif condition_func and condition_func(record):
                yield None, record

    def read_record(self):
        for i, row_or_record in enumerate(self.in_iter):
            if isinstance(row_or_record, dict):
                row_or_record['sample'] = self.sample
                record = self.record_cls(args=row_or_record)
                meta = None

            elif isinstance(row_or_record, tuple):
                meta, record = row_or_record
                if isinstance(record, dict):
                    record['sample'] = self.sample
                    record = self.record_cls(args=record)
                elif not isinstance(record, self.record_cls):
                    record.sample = self.sample
                    record = self.record_cls(parent=record)
            elif not isinstance(row_or_record, self.record_cls):
                meta = None
                row_or_record.sample = self.sample
                record = self.record_cls(parent=row_or_record)
            else:
                meta = None
                record.sample = self.sample
                record = row_or_record

            yield meta, record

    def read_chunks(self, field):
        chunks = []
        prev_value = None

        for i, record in enumerate(self.read_record()):
            meta, record = record

            value = getattr(record, field)

            if prev_value is None or prev_value == value:
                chunks.append(record)
            else:  # update
                if i:
                    yield prev_value, chunks
                chunks = [record]
            prev_value = value

        if chunks:
            yield prev_value, chunks
