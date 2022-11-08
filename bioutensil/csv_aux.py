import csv
import gzip
import sys

csv.field_size_limit(sys.maxsize)

# --- load mapping function --- #
def str2int(x):
    try: x = int(x)
    except: x = 0
    return x

def str2float(x):
    try: x = float(x)
    except: x = 0.0
    return x

def str2list(sep=','):
    def func(x):
        return [] if x == '.' else x.split(sep)
    return func

# --- dump mapping function --- #
list2str = lambda l: ','.join([ str(x) for x in l]) if l else '.'

# --- reduce mapping function --- #
right = lambda x, y: y
left = lambda x, y: x
add = lambda x, y: x + y
list_append = lambda l, x: [x] if l == 'empty_list' else l + [x]

def assign_map(func2headers):
    _map = {}
    for func, headers in func2headers.items():
        for h in headers: _map[h] = func
    return _map

def assign_reduce_map(func_init2headers):
    _map = {}
    for func_init, headers in func_init2headers.items():
        func, init = func_init
        for h in headers: _map[(h,init)] = func
    return _map


class csvReader(object):

    def __init__(self, in_fn, load_map, headers=None, sep=',', reduce_map={}, has_header=False):
        self.in_fn = in_fn
        self.load_map = load_map
        self.headers = headers
        self.sep = sep
        self.reduce_map = reduce_map
        self.has_header = has_header

    def _load_row(self, row, filtered_funcs=[], required_headers='all'):
        error = []
        for header, filtered_func in filtered_funcs:
            if header in row:
                value = self.load_map.get(header, lambda x:x)(row[header])
                if not filtered_func(value):
                    error.append('%s=%s' % (header, str(value)))
        if error:
            return 'Error: '+','.join(error)

        if required_headers == 'all': required_headers = self.headers
        if not required_headers: required_headers = row.keys()

        data = {}
        for header in required_headers:
            value = row[header]
            value = self.load_map.get(header, lambda x: x)(value)
            data[header] = value
        return data

    def read_rows(self, filtered_funcs=[], required_headers='all'):

        if self.in_fn.endswith('.gz'):
            f = gzip.open(self.in_fn, 'rb')
        else:
            f = open(self.in_fn, 'rb')

        collection = csv.DictReader(f, fieldnames=self.headers, delimiter=self.sep, skipinitialspace=True)
        for i, row in enumerate(collection):
            if i == 0 and self.has_header: pass
            else:
                row = self._load_row(row, filtered_funcs, required_headers)
                if 'Error' not in row:
                    if self.has_header: i-=1
                    yield (i, row)

    def read_chunks(self, header, filtered_funcs=[], required_headers='all', reduce=False):
        chunks = []
        prev_value = None
        chunk_i = 0
        for i, row in self.read_rows(filtered_funcs, required_headers):
            value = row[header]
            if prev_value is None or prev_value == value: chunks.append(row)
            else: # update
                if i:
                    yield (prev_value, self.reduce_chunks(chunks, reduce))
                    chunk_i += 1
                chunks = [row]
            prev_value = value

        if chunks: yield (prev_value, self.reduce_chunks(chunks, reduce))

    def reduce_chunks(self, chunks, reduce):
        if not reduce: return chunks

        reduce_map = {  head_init:func for head_init, func in self.reduce_map.items() if head_init[0] in chunks[0].keys()}

        res = { header:init for header, init in reduce_map.keys() }
        for row in chunks:
            for header_init, func in reduce_map.items():
                header, _ = header_init
                res[header] = func(res[header], row[header])

        return [res]

class csvWriter():

    def __init__(self, out_fn, dump_map, headers, sep=','):
        self.out_fn = out_fn
        self.dump_map = dump_map
        self.headers = headers
        self.sep = sep

        self._write_header()

    def _write_header(self):
        if self.out_fn.endswith('.gz'):
            f = gzip.open(self.out_fn, 'wb')
        else:
            f = open(self.out_fn, 'wb')

        line = self.sep.join(self.headers) + '\n'
        f.write(line)
        return True

    def _dump_row(self, row):
        data = {}
        for header, value in row.items():
            value = self.dump_map.get(header, str)(value)
            data[header] = value
        return data

    def write_row(self, row, mod='ab'):


        if self.out_fn.endswith('.gz'):
            f = gzip.open(self.out_fn, mod)
        else:
            f = open(self.out_fn, mod)

        row = self._dump_row(row)
        # write
        line = self.sep.join([ row[header] for header in self.headers]) + '\n'
        f.write(line)
        f.close()
        return True

    def write_rows(self, rows, mod='ab'):

        if self.out_fn.endswith('.gz'):
            f = gzip.open(self.out_fn, mod)
        else:
            f = open(self.out_fn, mod)

        for row in rows:
            row = self._dump_row(row)
            # write
            line = self.sep.join([ row[header] for header in self.headers]) + '\n'
            f.write(line)

        f.close()
        return True
