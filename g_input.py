import pysam
import argparse
import sv
def depth2cn(sampleDepth, wgsDepth, purity):
    ploidy = 2
    haploDepth = (wgsDepth*purity/ploidy)
    return sampleDepth/haploDepth
def generate_seg():
    import pysam
    parser = argparse.ArgumentParser(description='Generate SEG file from SV and BAM files.')
    parser.add_argument('-sv', '--sv_file', dest='svPath', required=True, help='Path to SV file')
    parser.add_argument('--tool', dest='tool', required=False, help='SV calling tool, manta or svaba')
    parser.add_argument('-bam', '--bam_file', dest='bamPath', required=False, help='Path to BAM file')
    parser.add_argument('-p', '--purity', required=False, type=int, default= 1,help='Tumor purity')
    parser.add_argument('-s', '--sample_name', dest='sampleName', required=False, default = 'sample', help='Sample name for output file')
    parser.add_argument('-d', '--wgs_depth', dest='wgsDepth', required=False, type=int, default = 0, help='The whole genome average depth')
    args = parser.parse_args()
    # find all breakpoints on each chromosome
    all_sv, pos = [], {}
    sv_records = sv.read_vcf(args.svPath, tool=args.tool, precise=False)
    for i, record in enumerate(sv_records):
        info = []
        # chrom_5p,bkpos_5p,strand_5p,chrom_3p,bkpos_3p,strand_3p
        info.append(record.chrom_5p)
        info.append(record.bkpos_5p)
        info.append(record.strand_5p)
        info.append(record.chrom_3p)
        info.append(record.bkpos_3p)
        info.append(record.strand_3p)

        # f.write(str(record) + '\n')
    # for line in open(args.svPath, 'r').readlines()[1:]:
        # info = line.strip('\n').split('\t')
        info[1], info[4] = int(info[1]), int(info[4])
        all_sv.append(info)
        if info[0] in pos.keys():
            pos[info[0]].append(info[1])
        else:
            pos[info[0]] = [info[1]]
        if info[3] in pos.keys():
            pos[info[3]].append(info[4])
        else:
            pos[info[3]] = [info[4]]
    for key, arr in pos.items():
        pos[key] = list(set(arr))
        pos[key].sort()
        pos[key].insert(0, max(1, pos[key][0]-1000))
        pos[key].append(pos[key][-1]+1000)

    # call coverage depth of segments from BAM
    segDepth = {}
    if args.bamPath == None:
        for chr, arr in pos.items():
            for n in range(1, len(arr)):
                key = chr+':'+str(arr[n-1])+'-'+str(arr[n])
                segDepth[key] = 100
    else:
        bam = pysam.AlignmentFile(args.bamPath, 'rb')
        for key, arr in pos.items():
            for n in range(1, len(arr)):
                cnt = bam.count_coverage(key, arr[n-1], arr[n], quality_threshold = 0)
                posDepth = []
                for i in range(len(cnt[0])):
                    temp = 0
                    for j in range(4):
                        temp += cnt[j][i]
                    posDepth.append(temp)
                name = key+':'+str(arr[n-1])+'-'+str(arr[n])
                segDepth[name] = sum(posDepth)/len(posDepth) # average on read depth of all positions
    # output seg.txt file
    if args.wgsDepth!=0 and args.purity!=0:
        for key, value in segDepth.items():
            segDepth[key] = depth2cn(value, args.wgsDepth, args.purity)
            # print(key, segDepth[key])
    resStr = ''
    for key, value in segDepth.items():
        resStr += key+'\t'+str(value)+'\n'
    segFile = open('{}_seg.txt'.format(args.sampleName), 'w')
    segFile.write(resStr)
    segFile.close()

if __name__ == "__main__":
    generate_seg()