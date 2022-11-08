import os, sys, subprocess

def sort_files(fns, sorted_col, header=True):
    for fn in fns:
        sort_file(fn, sorted_col, header=header)
    return True

def sort_file(fn, sorted_col, header=True):
    path, name = os.path.split(fn)
    sorted_fn = os.path.join(path, name.replace('tsv','sortBy.%s.tsv' % str(sorted_col)))
    cmd = '(zcat %s | head -1; zcat %s | sed "1d" | sort ) | gzip -c > %s' % (fn, fn, sorted_fn)
    #print cmd
    subprocess.call(cmd, shell=True)
    #cmd = 'rm -f %s' % (fn)
    #subprocess.call(cmd, shell=True)
    return True

