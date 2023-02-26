import argparse
import sys

class MainArgParser:
    def __init__(self):
        parser = argparse.ArgumentParser(prog='prelocalhap')
        parser.add_argument(dest='subfunc', help='Subcommands: ')
        args = parser.parse_args(sys.argv[1:2])
        getattr(self, args.subfunc)()
   
    def mergeFASTA(self):
        parser = argparse.ArgumentParser(description='Merge sequences in the fasta file')
        parser.add_argument('-fa', '--fasta',
                            dest='Fasta',
                            required=True,
                            help='Input fasta file')
        args = parser.parse_args(sys.argv[2:])
        res = ">ecDNA\n"
        for line in open(args.Fasta, "r").readlines():
            if line[0] == '>':
                continue
            res += line.strip('\n')
        fasta = open(args.Fasta, "w")
        fasta.write(res)

    def vcf2sv(self):
        parser = argparse.ArgumentParser(description='Extract SV information from VCF file')
        parser.add_argument('-v', '--vcf',
                            dest='vcf',
                            required=True,
                            help='Input VCF file')
        parser.add_argument('-p', '--prefix',
                            dest='prefix',
                            required=True,
                            help='Prefix of output file')
        args = parser.parse_args(sys.argv[2:])
        sv = []
        vaf, A, D = [], [], []
        ssm = []
        for line in open(args.vcf, "r").readlines():
            if line.startswith('#'):
                continue
            info = line.strip('\n').split('\t')
            prop = info[7].split(';')
            if prop[-1].split('=')[1] == 'BND':
                end = info[4].split('[')
                str1, str2 = '+', '+'
                if ']' in info[4]:
                    end = info[4].split(']')
                    str2 = '-'
                if end[0] == '':
                    str1 = '-'
                chr1, bkp1 = info[0], info[1]
                chr2, bkp2 = end[1].split(':')
                key, num = info[8].split(':'), info[12].split(':')
                data = {}
                for i in range(len(key)):
                    data[key[i]] = num[i]
                ad, dp = int(data['AD']), int(data['DP'])
                A.append(data['AD'])
                D.append(str(ad+dp))
                sv.append([chr1, bkp1, str1, chr2, bkp2, str2, data['AD']])
                vaf.append([chr1[3:], bkp1, data['DP'], data['AD'], str(ad/(dp+ad)*100)])
                vaf.append([chr2[3:], bkp2, data['DP'], data['AD'], str(ad/(dp+ad)*100)])
                id = '{}:{}:{}:{}:{}:{}'.format(chr1, bkp1, str1, chr2, bkp2, str2)
                ssm.append(['s{}'.format(len(ssm)), id, data['AD'], data['DP'], '0.999', '0.5'])
        with open('{}_sv.txt'.format(args.prefix), 'w') as sv_file:
            sv_file.write('chr_3p	bkp_3p	str_3p	chr_5p	bkp_5p	str_5p	depth\n')
            for info in sv:
                sv_file.write('\t'.join(info)+'\n')
            sv_file.close()
        with open('{}.vafs'.format(args.prefix), 'w') as vafs_file:
            for info in vaf:
                vafs_file.write('\t'.join(info)+'\n')
            vafs_file.close()
        with open('{}.input'.format(args.prefix), 'w') as input_file:
            input_file.write('> A\n({}, {})\n'.format(len(A), 1))
            for a in A:
                input_file.write(a + '\n')
            input_file.write('\n> D\n({}, {})\n'.format(len(D), 1))
            for d in D:
                input_file.write(d + '\n')
            input_file.close()

        # with open('{}_ssm.txt'.format(args.prefix), 'w') as ssm_file:
        #     ssm_file.write('id	gene	a	d	mu_r	mu_v\n')
        #     for info in ssm:
        #         ssm_file.write('\t'.join(info)+'\n')
        #     ssm_file.close()

    def multivcf2sv(self):
        parser = argparse.ArgumentParser(description='Extract SV information from multi-sample VCF file')
        parser.add_argument('-v', '--vcf',
                            dest='vcf',
                            required=True,
                            help='Input VCF files')
        parser.add_argument('-p', '--prefix',
                            dest='prefix',
                            required=True,
                            help='Prefix of output file')
        args = parser.parse_args(sys.argv[2:])
        vcfs = args.vcf.split(',')
        vcf_num = len(vcfs)
        sv_key = []
        sv, A, D = [{} for i in range(vcf_num)], [{} for i in range(vcf_num)], [{} for i in range(vcf_num)]
        vaf = [{} for i in range(vcf_num)]
        for n in range(vcf_num):
            sv_key_temp = []
            for line in open(vcfs[n], "r").readlines():
                if line.startswith('#'):
                    continue
                info = line.strip('\n').split('\t')
                prop = info[7].split(';')
                if prop[-1].split('=')[1] == 'BND':
                    end = info[4].split('[')
                    str1, str2 = '+', '+'
                    if ']' in info[4]:
                        end = info[4].split(']')
                        str2 = '-'
                    if end[0] == '':
                        str1 = '-'
                    chr1, bkp1 = info[0], info[1]
                    chr2, bkp2 = end[1].split(':')
                    # data for single sample
                    key, num = info[8].split(':'), info[12].split(':')
                    data = {}
                    for i in range(len(key)):
                        data[key[i]] = num[i]
                    ad, dp = int(data['AD']), int(data['DP'])

                    key1 = ':'.join([chr1,bkp1,str1,chr2,bkp2,str2])
                    if len(sv_key) == 0 or key1 in sv_key:
                        if key1 not in sv_key_temp:
                            sv_key_temp.append(key1)
                    if key1 not in sv[n].keys() or int(sv[n][key1][-1]) < ad:
                        sv[n][key1] = [chr1, bkp1, str1, chr2, bkp2, str2, data['AD']]
                        A[n][key1], D[n][key1] = data['AD'], str(ad+dp)
                        vaf[n][key1] = [chr1[3:], bkp1, data['DP'], data['AD'], str(ad/(dp+ad)*100)]
                    # complementary sv
                    if str1 == str2:
                        if str1 == '+':
                            str1, str2 = '-', '-'
                        else:
                            str1, str2 = '+', '+'
                    key2 = ':'.join([chr2,bkp2,str1,chr1,bkp1,str2])
                    if len(sv_key) == 0 or key2 in sv_key:
                        if key2 not in sv_key_temp:
                            sv_key_temp.append(key2)
                    sv[n][key2] = sv[n][key1]
                    A[n][key2], D[n][key2] = A[n][key1], D[n][key1]
                    vaf[n][key2] = [chr2[3:], bkp2, data['DP'], data['AD'], str(ad/(dp+ad)*100)]
            sv_key = sv_key_temp.copy()

        sv_key = list(filter(lambda key: (sv_key.index(key)%2==0), sv_key))
        for n in range(vcf_num):
            with open('{}.vafs'.format(args.prefix+str(n)), 'w') as vafs_file:
                for key in sv_key:
                    vafs_file.write('\t'.join(vaf[n][key])+'\n')
                vafs_file.close()
        for n in range(vcf_num):
            with open('{}_sv.txt'.format(args.prefix+str(n)), 'w') as sv_file:
                sv_file.write('chr_3p	bkp_3p	str_3p	chr_5p	bkp_5p	str_5p	depth\n')
                for key in sv_key:
                    sv_file.write('\t'.join(sv[n][key])+'\n')
                sv_file.close()
        with open('{}.input'.format(args.prefix), 'w') as input_file:
            input_file.write('> A\n({}, {})\n'.format(len(sv_key), vcf_num))
            for key in sv_key:
                for n in range(vcf_num-1):
                    input_file.write(A[n][key]+'.\t')
                input_file.write(A[-1][key]+'\n')
            input_file.write('\n> D\n({}, {})\n'.format(len(sv_key), vcf_num))
            for key in sv_key:
                for n in range(vcf_num-1):
                    input_file.write(D[n][key]+'.\t')
                input_file.write(D[-1][key]+'\n')
            input_file.close()

if __name__ == '__main__':
    MainArgParser()

