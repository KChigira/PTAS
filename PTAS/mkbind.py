#!/usr/bin/env python3

import csv
import pandas as pd
import sys
import subprocess as sbp
from PTAS.utils import read_vcf, time_stamp, prepare_cmd
from PTAS.params import Params
from PTAS.visualize_marker import VisualizeMarker

pm = Params('mkbind')
args = pm.set_options()

class MKBind(object):
    def __init__(self, args):
        self.args = args
        self.fai = args.fai
        self.vcf = args.vcf
        self.out = args.output

        #VCF data
        self.header = [] #VCF header
        self.colnames = [] #VCF colnames
        self.data = [] #VCF content
        self.primers = [] #Primer information of selected variants

    def readvcf(self):
        #Read first VCF
        vcf_list = read_vcf(self.vcf[0])
        self.header = vcf_list[0]
        self.colnames = vcf_list[1]
        self.data = vcf_list[2]
        self.header.append('##mkbind.py_{},source_vcf={}'.format(time_stamp(), self.vcf))
        self.data = pd.DataFrame(self.data, columns=self.colnames)
        for i in range(len(self.data)):
            #column 7 is INFO
            #add the number of source VCF file
            self.data.iat[i, 7] += ';VCF=0'

        if len(self.vcf) == 1:
            return
        
        #Add 2nd ~ VCF
        for i in range(1, len(self.vcf)):
            vcf_list = read_vcf(self.vcf[i])
            sub_data = vcf_list[2]
            sub_data = pd.DataFrame(sub_data, columns=self.colnames)
            for j in range(len(sub_data)):
                #column 7 is INFO
                #add the number of source VCF file
                sub_data.iat[j, 7] += ';VCF={}'.format(i)
            #Bind rows
            self.data = pd.concat([self.data, sub_data], axis=0)

        #Sort binded VCF
        self.data['POS'] = self.data['POS'].astype(int)
        self.data = self.data.sort_values(['#CHROM', 'POS'],
                                          ascending=[True, True])
        self.data = self.data.reset_index(drop=True)

        #Correct infomation of primers
        self.primers = self.get_info_primers(self.data)

    def get_info_primers(self, vcfdata):
        primers_list = []
        for i in range(len(vcfdata)):
            info = str(vcfdata.at[vcfdata.index[i], 'INFO'])
            info_spl = info.split(';')
            for j in range(len(info_spl)):
                elem = info_spl[j].split('=')
                if elem[0] == 'PRIMER':
                    pri = elem[1].split('|')
                    primers_list.append(pri)
        primers_col = ['name','num','chr','L_seq','R_seq','L_pos','R_pos','L_TM','R_TM','product_size']
        return pd.DataFrame(primers_list, columns=primers_col)

    def maketable(self):
        name_line = self.data.columns[9:]
        out_table_col = ['Chr','Pos','ID','Left','Right','L_TM','R_TM']
        for str_n in name_line:
            out_table_col.append('{}_Size'.format(str_n))

        out_list = []
        for i in range(len(self.data)):
            out_row = []
            out_row.append(self.data.at[self.data.index[i], '#CHROM'])
            out_row.append(self.data.at[self.data.index[i], 'POS'])
            out_row.append(self.primers.at[self.primers.index[i], 'name'])
            out_row.append(self.primers.at[self.primers.index[i], 'L_seq'])
            out_row.append(self.primers.at[self.primers.index[i], 'R_seq'])
            out_row.append(self.primers.at[self.primers.index[i], 'L_TM'])
            out_row.append(self.primers.at[self.primers.index[i], 'R_TM'])
            for line in name_line:
                line_info = self.data.at[self.data.index[i], line]
                #ex. 1/1:0,6:6:18:270,18,0
                if line_info[0:3] == '0/0':
                    out_row.append(self.primers.at[self.primers.index[i], 'product_size'])
                elif line_info[0:3] == '1/1':
                    size = (int(self.primers.at[self.primers.index[i], 'product_size']) +
                            len(self.data.at[self.data.index[i], 'ALT']) - 
                            len(self.data.at[self.data.index[i], 'REF']))
                    out_row.append(size)
            out_list.append(out_row)
        self.out_table = pd.DataFrame(out_list, columns=out_table_col)

    def output(self):
        #1. output VCF
        out_vcf_name = '{}.vcf'.format(self.out)
        with open(out_vcf_name, 'w') as o:
            for h in self.header:
                o.write('{}\n'.format(h))
            writer = csv.writer(o, delimiter='\t')
            writer.writerow(self.colnames)

        out_vcf_name_tmp = '{}.tmp'.format(out_vcf_name)
        self.data.to_csv(out_vcf_name_tmp, sep='\t', header=False, index=False)
        cmd1 = 'cat {} >> {}'.format(out_vcf_name_tmp, out_vcf_name)
        cmd2 = 'rm {}'.format(out_vcf_name_tmp)
        cmd1 = prepare_cmd(cmd1)
        cmd2 = prepare_cmd(cmd2)
        try:
            sbp.run(cmd1,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True, check=True)
            sbp.run(cmd2,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True, check=True)
        except sbp.CalledProcessError:
            print(time_stamp(),
              'Error occured writing vcf data.',
              flush=True)
            sys.exit(1) 

        #2. output primer data
        self.out_table.to_csv('{}.txt'.format(self.out), sep='\t', header=True, index=False)

    def draw(self):
        out_vcf_name = '{}.vcf'.format(self.out)
        out_png_name = '{}.png'.format(self.out)

        vm = VisualizeMarker(out_vcf_name, out_png_name, self.fai, True)
        vm.run()
        

def main():
    print(time_stamp(), 'MKBind started.', flush=True)

    prog = MKBind(args)
    prog.readvcf()
    prog.maketable()
    prog.output()
    prog.draw()

    print(time_stamp(), 'MKBind successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()
    
