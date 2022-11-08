# -*- coding: utf-8 -*-
"""
    tenxtools.module
    ~~~~~~~~~~~~~~~~

    @Copyright: (c) 2018-04 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.
"""
import os
from itertools import tee


class Module(object):

    def __init__(self, in_data=None, in_fn=None,
                 reader_cls=None, file_type=None,
                 record_cls=None,
                 write=False, write_stats=True, write_headers=False,
                 write_fields=True, dynamic_fields=False,
                 writer_cls=None, write_mode='whole', out_dir='',
                 prefix='', suffix='', chunk_size=1000,
                 sample=None,
                 **extra):
        self.sample = sample

        self.in_fn = in_fn
        self.in_data = in_data
        self.reader_cls = reader_cls
        self.file_type = file_type
        self.record_cls = record_cls

        if in_data:
            if file_type:
                self.in_iter = self.reader_cls(in_data=in_data, sample=sample, file_type=file_type)
            else:
                self.in_iter = self.reader_cls(in_data=in_data, sample=sample)
        elif in_fn:
            if file_type:
                self.in_iter = self.reader_cls(in_fn=in_fn, sample=sample, file_type=file_type)
            else:
                self.in_iter = self.reader_cls(in_fn=in_fn, sample=sample)
        else:
            self.in_iter = in_data

        self.write = write
        self.write_stats = write_stats
        self.write_headers = write_headers
        self.write_fields = write_fields
        self.dynamic_fields = dynamic_fields
        self.writer_cls = writer_cls
        self.out_fn = None
        self.out_fns = []

        self.write_mode = write_mode
        self.out_dir = out_dir
        self.prefix = prefix
        self.suffix = suffix
        self.chunk_size = chunk_size

        self.meta_dict = {}
        self.writer_dict = {}
        self.meta_count = {}

    def evaluate(self):
        gen1, gen2 = tee(self)

        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        for meta, record in gen1:
            if self.write:
                self._write_line(meta, record)

        if self.write:
            self._write_after_iter()

        return gen2

    def __iter__(self):
        pass

    def _write_line(self, meta, line):
        # update count
        count = self.meta_count.get(meta, 0)
        self.meta_count[meta] = count + 1

        # write
        if self.write_mode == 'whole':
            meta = 'ALL'
        chunk = self.meta_dict.get(meta, [])

        if len(chunk) >= self.chunk_size:
            self._write_chunk(meta, chunk)
            self.meta_dict[meta] = [line]
        else:
            chunk.append(line)
            self.meta_dict[meta] = chunk

    def _write_chunk(self, meta, chunk):
        dynamic_fields = False
        if meta not in self.writer_dict:
            if meta == 'ALL':
                fns = [self.prefix, self.suffix]
            else:
                fns = [self.prefix, meta, self.suffix]
            fns = [x for x in fns if x]
            self.out_fn = os.path.join(self.out_dir, '.'.join(fns))
            if self.out_fn not in self.out_fns:
                self.out_fns.append(self.out_fn)
            writer = self.writer_cls(fn=self.out_fn)
            if self.write_headers:
                writer.write_headers()
            if self.write_fields:
                if not self.dynamic_fields:
                    writer.write_fields()
                else:
                    dynamic_fields = True
            self.writer_dict[meta] = writer
        else:
            writer = self.writer_dict[meta]

        writer.write_chunk(chunk, dynamic_fields=dynamic_fields)

    def _write_after_iter(self):
        self._write_last_chunk()
        if self.write_stats:
            self._write_stat()

    def _write_last_chunk(self):
        for meta, chunk in self.meta_dict.items():
            self._write_chunk(meta, chunk)

    def _write_stat(self):
        stat_fn = os.path.join(
            self.out_dir, self.prefix + '.' + self.suffix + '.stat')
        with open(stat_fn, 'w') as f:
            f.write('\t'.join(map(str, self.meta_count.keys())) + '\n')
            f.write('\t'.join(map(str, self.meta_count.values())) + '\n')
