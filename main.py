#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    SVAS.complexsv
    ~~~~~~~~~~~~~~~~~~~~~~~

    @Input: DESCRIPTION
    @Return:

    @Copyright: (c) 2020-05 by Lingxi Chen (chanlingxi@gmail.com).
    @License: LICENSE_NAME, see LICENSE for more details.

Usage:
    complexsv.py call --sv_fn=IN_FILE --out_dir=OUT_DIR --tool=TOOL --max_dis=max_dis --max_range=max_range --cnv_fn=CNV_FN
    complexsv.py cnv --sv_fn=IN_FILE --out_dir=OUT_DIR --tool=TOOL --max_dis=max_dis --max_range=max_range --cnv_fn=CNV_FN
    complexsv.py -h | --help

Options:
    -h --help           Show this screen.
    --version           Show version.
    --sv_fn=IN_FILE     Path of SvABA sv vcf file.
    --out_dir=OUT_DIR   Path of output directory.
"""

import os
from cmath import inf
import argparse
import sys

from bioutensil import constants

import sv

import base
import  bplot_complexsv
import chromothripsis

MIN_SV_NUM = 6

class ComplexSVRegionGroupGenerator():

    def __init__(self, sv_records=None, groups=None, bam_fn=None, *args, **kwargs):
        self.sv_records = sv_records or []
        self.groups = groups or []
        self.bam_fn = bam_fn
        self.args = kwargs
        # print(kwargs)

    def compare(self, sv1, sv2):
        if int(sv1.bkpos_5p) < int(sv2.bkpos_5p):
            return -1
        elif int(sv1.bkpos_5p) > int(sv2.bkpos_5p):
            return 1
        if int(sv1.bkpos_3p) < int(sv2.bkpos_3p):
            return -1
        elif int(sv1.bkpos_3p) > int(sv2.bkpos_3p):
            return 1
        return 0

    def setRange(self, chr_range, sv):
        if sv.chrom_5p in chr_range.keys():
            chr_range[sv.chrom_5p][0] = min(chr_range[sv.chrom_5p][0], int(sv.bkpos_5p))
            chr_range[sv.chrom_5p][1] = max(chr_range[sv.chrom_5p][1], int(sv.bkpos_5p))
        else:
            chr_range[sv.chrom_5p] = [int(sv.bkpos_5p), int(sv.bkpos_5p)]
        if sv.chrom_3p in chr_range.keys():
            chr_range[sv.chrom_3p][0] = min(chr_range[sv.chrom_3p][0], int(sv.bkpos_3p))
            chr_range[sv.chrom_3p][1] = max(chr_range[sv.chrom_3p][1], int(sv.bkpos_3p))
        else:
            chr_range[sv.chrom_3p] = [int(sv.bkpos_3p), int(sv.bkpos_3p)]
        return chr_range

    def min_dis(self, sv1, sv2):
        diff1 = abs(int(sv1.bkpos_5p)-int(sv2.bkpos_5p))
        if sv1.chrom_5p != sv2.chrom_5p:
            diff1 = inf
        diff2 = abs(int(sv1.bkpos_5p)-int(sv2.bkpos_3p))
        if sv1.chrom_5p != sv2.chrom_3p:
            diff2 = inf
        diff3 = abs(int(sv1.bkpos_3p)-int(sv2.bkpos_5p))
        if sv1.chrom_3p != sv2.chrom_5p:
            diff3 = inf
        diff4 = abs(int(sv1.bkpos_3p)-int(sv2.bkpos_3p))
        if sv1.chrom_3p != sv2.chrom_3p:
            diff4 = inf
        return min(diff1,diff2,diff3,diff4)

    def check_range(self, chr_range, max_range):
        for val in chr_range.values():
            if val[1]-val[0] > max_range:
                return False
        return True

    def _group_sv(self):
        juncs = []
        for i, record in enumerate(self.sv_records):
            if record.chrom_5p not in constants.chrs or record.chrom_3p not in constants.chrs:
                continue
            # info = [record.chrom_5p, record.bkpos_5p, record.strand_5p, record.chrom_3p, record.bkpos_3p, record.strand_3p]
            if record.strand_5p == '-' and record.strand_3p == '-':
                record.strand_5p, record.strand_3p = '+','+'
                record.chrom_5p, record.chrom_3p = record.chrom_3p, record.chrom_5p
                record.bkpos_5p, record.bkpos_3p = record.bkpos_3p, record.bkpos_5p
            # if info[2] == '-' and info[5] == '-':
            #     info[2], info[5] = '+', '+'
            #     info[0], info[3] = info[3], info[0]
            #     info[1], info[4] = info[4], info[1]
            juncs.append(record)
        # sort juncs in chromosome
        juncs = sorted(juncs, key=lambda x: (x.bkpos_5p, self.compare))
        # cluster SV junctions based on distance
        # cluster = []
        svIdx = list(range(0, len(juncs))) # sv index for selection
        while len(svIdx) > 0:
            subcluster = base.RegionGroup(linkage_distance = int(self.args['max_dis']))
            subcluster.append_sv(juncs[svIdx[0]])
            queue = [svIdx[0]] # index of sv_info
            sv, chr_range = juncs[svIdx[0]], {}
            self.setRange(chr_range, sv)
            svIdx.pop(0)
            while len(queue) > 0:
                idx = queue[0]
                queue.pop(0)
                for i in svIdx:
                    if self.min_dis(juncs[i], juncs[idx]) < int(self.args['max_dis']):
                        temp_range = chr_range.copy()
                        self.setRange(temp_range, juncs[i])
                        if self.check_range(temp_range, int(self.args['max_range'])):
                            self.setRange(chr_range, juncs[i])
                            queue.append(i)
                            subcluster.append_sv(juncs[i])
                            svIdx.remove(i)
                            # print(subcluster)
            # if self.hasFBI(subcluster, juncs) == True:
            self.groups.append(subcluster)
            for r in subcluster.region_list:
                self._assign_cn_list(self.args['cn_fn'],r)
        # print(self.groups)
    def call(self):
        self._group_sv()
        # self._cluster()
        cnt_dict = {}
        for i, group in enumerate(self.groups):
            # filter
            group_type = group.group_type

            if group_type in cnt_dict:
                cnt_dict[group_type] += 1
            else:
                cnt_dict[group_type] = 1

            group.group_name = '{}{}'.format(group_type, cnt_dict[group_type])
            group.expand_region()

            # linkage information
            length = 0
            for region in group.region_list:
                length += region.end - region.start

            '''
            if length <= 10000:
                linkage_heatmap.run_call(bam_fn=self.bam_fn,
                                         regions=group.region_list,
                                         out_dir=out_dir)
            '''
            # self._write_group(group)
            yield group

    def _write_group(self, group):
        out_dir = os.path.join(self.out_dir, group.group_type, group.group_name)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        sv_ids = set()
        out_prefix = os.path.join(out_dir, '{}.{}'.format(self.sample, group.group_name))
        with open(out_prefix + '.sv', 'w') as sv_f, open(out_prefix + '.region', 'w') as region_f:
            for region in group.region_list:
                region_f.write('{}:{}-{}\n'.format(region.chrom, region.start, region.end))
                # write sv
                for sv_re in region.sv_list:
                    if sv_re.id not in sv_ids:
                        sv_f.write('{}\n'.format(sv_re))
                        sv_ids.add(sv_re.id)

    def _assign_cn_list(self, cn_fn, region):
        if not os.path.isfile(cn_fn):
            return
        for line in open(cn_fn, 'r'):
            region_str, cn = line.strip().split('\t')
            cn = float(cn)
            chrom = region_str.split(':')[0]
            start = int(region_str.split(':')[1].split('-')[0])
            end = int(region_str.split(':')[1].split('-')[1])

            if region.chrom == chrom and start >= region.start and end <= region.end and end - start > 1:
                region.cn_list.append(base.CN(chrom=chrom, start=start, end=end, cn=cn))

    def read(self):
        for group_type in os.listdir(os.path.join(self.out_dir)):
            if not os.path.isdir(group_type):
                continue
            for group_id in os.listdir(os.path.join(self.out_dir, group_type)):
                in_dir = os.path.join(self.out_dir, group_type, group_id)
                prefix_fn = os.path.join(in_dir, '{}.{}'.format(self.sample, group_id))
                sv_fn = prefix_fn + '.sv'
                cn_fn = prefix_fn + '.cn'
                enhancer_fn = prefix_fn + '.enhancer'
                super_enhancer_fn = prefix_fn + '.super_enhancer'
                region_fn = os.path.join(in_dir, '{}.{}.region'.format(self.sample, group_id))
                region_list = []
                for line in open(region_fn, 'r'):
                    chrom = line.split(':')[0]
                    start = int(line.split(':')[1].split('-')[0])
                    end = int(line.split(':')[1].split('-')[1])
                    region = base.Region(chrom=chrom, start=start, end=end, sv_fn=sv_fn, sample=self.sample)
                    self._assign_cn_list(cn_fn, region)
                    self._assign_enhancer_list(enhancer_fn, region)
                    self._assign_super_enhancer_list(super_enhancer_fn, region)
                    region_list.append(region)

                group = base.RegionGroup(region_list=region_list, group_name=group_id)
                yield group_type, group

    # def _group_sv(self, minimal_distance=1000):  # 1k
    #     for i, record in enumerate(self.sv_records):
    #         if i > 50:
    #             break
    #         if record.chrom_5p not in constants.chrs or record.chrom_3p not in constants.chrs:
    #             continue
    #         if record.chrom_5p == record.chrom_3p and \
    #                 abs(record.bkpos_5p - record.bkpos_3p) < minimal_distance:
    #             continue
    #         self._append_sv_to_groups(record)

    def _append_sv_to_groups(self, record):  # code review
        merge_groups = [group for group in self.groups if True in group.has_sv(record)[0]]
        groups = [group for group in self.groups if True not in group.has_sv(record)[0]]
        if merge_groups:
            group = merge_groups[0]
            group.append_sv(record)
            for g in merge_groups[1:]:
                g.append_sv(record)
                group.merge(g)
            groups.append(group)
        else:  # new group
            group = base.RegionGroup()
            group.append_sv(record)
            groups.append(group)
        self.groups = groups

    def _bkp_in_regions(self, bkp):
        bs = [bkp.inside_region(region) for region in self.regions]
        return True in bs, bs

    def _sv_in_regions(self, record):
        bkp_5p = base.Breakpoint(record.chrom_5p, record.bkpos_5p,
                                 record.bkpos_5p_orientation, record.orientation)
        bkp_3p = base.Breakpoint(record.chrom_3p, record.bkpos_3p,
                                 record.bkpos_3p_orientation, record.orientation)
        _, bs_5p = self._bkp_in_regions(bkp_5p)
        _, bs_3p = self._bkp_in_regions(bkp_3p)
        bs = [b1 or b2 for b1, b2 in zip(bs_5p, bs_3p)]
        return True in bs, bs

    def _append_sv_to_regions(self, record, bs):
        if True in bs:
            for i, b in enumerate(bs):
                if b:
                    self.regions[i].append_sv(record, self.tran_psl_reader)

def run_cnv(**args):
    chromothripsis.run_cnv_state(**args)

def run_csis(sv_file_type,**args):
    '''
    regions = [base.Region(chrom=chrom, start=0, end=constants.hg19_fai_bp[chrom]) for chrom in constants.chrs]
    print(regions)
    region = base.Region(chrom='chr7', start=13500000, end=15000000)
    groups = [base.RegionGroup(region_list=[region])]
    '''
    groups = []
    # sv_records = sv.parse(args['sv_fn'])
    if sv_file_type == "tsv":
        sv_records = sv.read_txt(args['sv_fn'])
    else:
        sv_records = sv.read_vcf(args['sv_fn'],tool=args['tool'])
    # print(sv_records)
    groups = ComplexSVRegionGroupGenerator(
        groups=groups,
        sv_records=sv_records,
        **args
    ).call()

    for g in groups:
        print([region.sv_ids for region in g.region_list])
        # if g.sv_num <= MIN_SV_NUM:
            # print()
            # continue
        # print('group', g.group_type, g.region_list,chromo_metrics["chromothripsis"])
        chromo_metrics = chromothripsis.evaluate_region(regions=g.region_list, search=False, **args)
        for item in chromo_metrics:
            if item["chromothripsis"]:
                print(g.region_list,item["chromothripsis"])
        # if chromo_metrics['cluster'] and chromo_metrics['random_walk']:
        # print(chromo_metrics["chromothripsis"])

    # groups = list(groups)
    # data = {}
    # for meta, group in groups:
    #     if meta in data:
    #         data[meta].append((meta, group))
    #     else:
    #         data[meta] = [(meta, group)]

    # for meta, group_list in data.items():
    #     out_prefix = os.path.join(args['out_dir'], '{}.{}.region'.format(args['sample'], meta))
    #     bplot_complexsv.run(draw=True,
    #                         groups=group_list,
    #                         out_prefix=out_prefix,
    #                         **args)
    # return groups, None


def run_draw(**args):
    reader_cls = base.RegionGroupReader
    groups = ComplexSVRegionGroupGenerator(
        reader_cls=reader_cls,
        write=False,  # True,
        **args
    ).evaluate()

    groups = list(groups)
    data = {}
    for meta, group in groups:
        if meta in data:
            data[meta].append((meta, group))
        else:
            data[meta] = [(meta, group)]

    for meta, group_list in data.items():
        out_fn = os.path.join(args['out_dir'], '{}.{}.region.svg'.format(args['sample'], meta))
        bplot_complexsv.run(draw=True,
                            sample=args['sample'],
                            groups=group_list,
                            out_fn=out_fn)
    return groups, None




if __name__ == "__main__":
    # complexsv.py call --sv_fn=IN_FILE --out_dir=OUT_DIR --tool=TOOL --max_dis=max_dis --max_range=max_range --cnv_fn=CNV_FN
    parser = argparse.ArgumentParser(description='Calling chromothricsis and chromoplexy from SV and CN files')
    parser.add_argument('--func', dest='func', required=True, default='csis', choices=['csis','other'], help='Sub functions')
    parser.add_argument('--sv_fn', dest='sv_fn', required=True, help='vcf or tsv (must end with .vcf or .tsv)')
    parser.add_argument('--cn_fn', dest='cn_fn', required=False, help='Segment CN file')
    parser.add_argument('--tool', dest='tool', required=False, help='SV calling tool, manta or svaba')
    parser.add_argument('--out_dir', dest='out_dir', required=False, help='Output dir')
    parser.add_argument('--max_dis', dest='max_dis', required=False, type=int, default=10000000, help='Maximum distance of two SVs grouped in a cluster')
    parser.add_argument('--max_range', dest='max_range', required=False, type=int, default=500000000, help='Maximum range of a cluster')

    args = parser.parse_args()
    
    if args.sv_fn.lower().endswith('tsv'):
        sv_fn_extention = "tsv"
    elif args.sv_fn.lower().endswith('vcf'):
        sv_fn_extention = "vcf"
    else:
        sys.exit("Error: sv_fn file should end with vcf or tsv file")
    run_csis(sv_file_type= sv_fn_extention, **vars(args))
