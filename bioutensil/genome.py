import os
import pandas as pd
from pybedtools import BedTool


def read_chrom_size(assembly):
    fn_dir = os.path.dirname(os.path.abspath(__file__))
    fn = os.path.join(fn_dir, 'bio_database/{}.chrom.size'.format(assembly))
    df = pd.read_csv(fn, sep='\t')
    return df.set_index('chrom')['size'].to_dict()


def read_cytoband(assembly):
    fn_dir = os.path.dirname(os.path.abspath(__file__))
    fn = os.path.join(fn_dir, 'bio_database/{}.cytoBand.txt'.format(assembly))
    df = pd.read_csv(fn, sep='\t')
    df['cytoband'] = df['#chr'].str.replace('chr', '') + df['cytoband']
    df['startPos'] = df['startPos'] + 1
    df['region'] = df['#chr'].astype(str) + ':' + df['startPos'].astype(str) + '-' + df['endPos'].astype(str)
    return df.set_index('cytoband')[['region', 'color']].to_dict()




def get_cytoband_from_pos(cytoband_dict, chrom, pos):
    bin_list = ['{}\t{}'.format(v, k) for k, v in cytoband_dict['region'].items()]
    bed_str = '\n'.join(bin_list).replace(':', '\t').replace('-', '\t')
    bin_bed = BedTool(bed_str, from_string=True)

    pos_bed = BedTool('{}\t{}\t{}'.format(chrom, pos, pos+1), from_string=True)

    for i, hit in enumerate(pos_bed.window(bin_bed).overlap(cols=[2, 3, 5, 6])):
        return hit[6]

def get_cytoband_from_bin_str(cytoband_dict, bin_str):
    chrom = bin_str.split(':')[0]
    start, end = map(int, bin_str.split(':')[1].split('-'))
    pos = int((start+end)/2)
    return get_cytoband_from_pos(cytoband_dict, chrom, pos)
