#!/usr/bin/env python


########################################################################################
#DESCRIPTION                                                                           #
#This script is designed to annotate splicing prediction based on ENIGMA cut-off point #
#The input needs to be annotated with MaxEntScan.pm and Intronic distance              # 
#https://github.com/Ensembl/VEP_plugins/blob/release/94/MaxEntScan.pm                  # 
#Contact: Jan.Shamsani@qimrberghofer.edu.au                                            #
#EXAMPLE USAGE                                                                         #
#module load python/3.6.1                                                              #
#python annotate_splicegenecity.py -i input.tsv.gz -o output.tsv                       #
########################################################################################


from textwrap import dedent

import argparse
import csv
import gzip
import os
import re
import sys
import platform
import itertools
import math
import hgvs.parser
csv.field_size_limit(sys.maxsize)
def get_option_parser():

    parser = argparse.ArgumentParser(
        prog = "summarize_VEP_TSV",
        description = dedent("""
        Summarize VEP output in TSV format.
        """.format(os.path.basename(__file__))))


    parser.add_argument('-i', '--inputfile',
        help = "Read VEP TAB output from FILE",
        metavar = "FILE")

    parser.add_argument('-o', '--output_basename',
        help = "Basename of output files to create",
        metavar = "STR")

    return parser


def main():

    parser = get_option_parser()
    args = parser.parse_args()

    Splicing_var_type = {}

    with gzip.open(args.inputfile, "rt") as infile, \
        gzip.open(args.output_basename + '.splicepredict.tsv.gz', 'wt') as outfile:
   
        tsv_reader = csv.reader(infile, delimiter='\t')

        hgvsparser = hgvs.parser.Parser()

        for row in tsv_reader:
            
            if row[0].startswith('##'):
                continue

            if row[0].startswith('#'):

                column_header = row
                outfile.write('\t'.join(column_header + ['Splicing_var_type','Native_loss','Donor_gain_native_site','Donor_gain']) + '\n')

                continue

            row_data = dict(zip(column_header, row))

            VARIANT_CLASS = row_data['VARIANT_CLASS']
            HGVSc = row_data['HGVSc']
            #REVEL = row_data['REVEL']
            Feature = row_data['Feature']
            Consequence = row_data['Consequence']
            IntronEndOffset = row_data['IntronEndOffset']
            IntronStartOffset = row_data['IntronStartOffset']
            MES_NCSS_downstream_acceptor = row_data['MES-NCSS_downstream_acceptor']
            MES_NCSS_upstream_acceptor = row_data['MES-NCSS_upstream_acceptor']
            MES_NCSS_downstream_donor = row_data['MES-NCSS_downstream_donor']
            MES_NCSS_upstream_donor = row_data['MES-NCSS_upstream_donor']
            MES_SWA_acceptor_alt = row_data['MES-SWA_acceptor_alt']
            MES_SWA_acceptor_ref = row_data['MES-SWA_acceptor_ref']
            MES_SWA_acceptor_ref_comp = row_data['MES-SWA_acceptor_ref_comp']
            MES_SWA_acceptor_diff = row_data['MES-SWA_acceptor_diff']
            MES_SWA_donor_alt = row_data['MES-SWA_donor_alt']
            MES_SWA_donor_ref = row_data['MES-SWA_donor_ref']
            MES_SWA_donor_ref_comp = row_data['MES-SWA_donor_ref_comp']
            MES_SWA_donor_diff = row_data['MES-SWA_donor_diff']
            MaxEntScan_alt = row_data['MaxEntScan_alt']
            MaxEntScan_diff = row_data['MaxEntScan_diff']
            MaxEntScan_ref = row_data['MaxEntScan_ref']
            
            #Splicing_var_type = row_data['Splicing_var_type']
            #IntronStartOffset = int(IntronStartOffset)
            #IntronEndOffset = int(IntronEndOffset)
            
            try: 
                IntronStartOffset = int(IntronStartOffset)
                IntronEndOffset = int(IntronEndOffset)
            except ValueError:
                continue
            #Annotate splice region according Burge et al
            
            if IntronStartOffset == 0 and IntronEndOffset ==0 and re.search('splice',Consequence):
                if MES_NCSS_downstream_donor == '-' or MES_NCSS_upstream_acceptor == '-' :
                    Splicing_var_type = 'last_exon'
                else :
                    MES_NCSS_downstream_acceptor = float(MES_NCSS_downstream_acceptor)
                    MES_NCSS_upstream_acceptor = float(MES_NCSS_upstream_acceptor)
                    MES_NCSS_downstream_donor = float(MES_NCSS_downstream_donor)
                    MES_NCSS_upstream_donor = float(MES_NCSS_upstream_donor)
                    MES_SWA_donor_ref = float(MES_SWA_donor_ref)
                    MES_SWA_acceptor_ref = float(MES_SWA_acceptor_ref)
                    if MES_SWA_donor_ref == MES_NCSS_downstream_donor:
                        Splicing_var_type = 'Exonic_donor_splice_region'
                    elif MES_SWA_acceptor_ref == MES_NCSS_upstream_acceptor:
                        Splicing_var_type = 'Exonic_acceptor_splice_region'
                    else:
                        Splicing_var_type = 'check'            
                     
            elif 0 < IntronStartOffset  <= 6:
                Splicing_var_type = 'Intronic_donor_splice_region'
            elif 0 > IntronStartOffset >= -20:
                Splicing_var_type = 'Intronic_acceptor_splice_region' 
            elif IntronStartOffset ==0 and IntronEndOffset >0:
                Splicing_var_type = 'Intronic_donor_splice_region'
            elif IntronStartOffset <0 and IntronEndOffset ==0:
                Splicing_var_type = 'Intronic_acceptor_splice_region'
            else:
                Splicing_var_type = 'Outside_native'
           

            #Annotate donor splice site 1st
            if Splicing_var_type == 'Intronic_donor_splice_region' or Splicing_var_type == 'Exonic_donor_splice_region' : 
                if re.search('SNV', VARIANT_CLASS):
                    if MaxEntScan_alt  == '-' or MaxEntScan_diff == '-':
                        Native_loss = '-'
                    else:
                        MaxEntScan_alt = float(MaxEntScan_alt)
                        MaxEntScan_diff = float(MaxEntScan_diff)

                        if MaxEntScan_diff > 0:    
                            if MaxEntScan_alt < 6.2:
                                if MaxEntScan_diff >= 1.15:
                                    Native_loss = 'HIGH'
                                else:
                                    Native_loss = 'MODERATE*'
                            elif MaxEntScan_alt > 8.5:
                                Native_loss = 'LOW'
                            else:
                                if MaxEntScan_diff >= 1.15:
                                    Native_loss = 'MODERATE'
                                else:
                                    Native_loss = 'LOW'   
                        else:
                            Native_loss = 'IMPROVED'
                else:
                    if MES_SWA_donor_diff == '-' or MES_SWA_donor_alt == '-':
                        Native_loss = '-'
                    else:  
                        MES_SWA_donor_alt = float(MES_SWA_donor_alt)
                        MES_SWA_donor_diff = float(MES_SWA_donor_diff)
         
                        if MES_SWA_donor_diff > 0 :
                            if MES_SWA_donor_alt < 6.2 :
                                if MES_SWA_donor_diff >= 1.15:
                                    Native_loss = 'HIGH'
                                else:
                                    Native_loss = 'MODERATE*'
                            elif MES_SWA_donor_alt > 8.5:
                                Native_loss = 'LOW'
                            else:
                                if MES_SWA_donor_diff >= 1.15: 
                                    Native_loss = 'MODERATE'
                                else: 
                                    Native_loss = 'LOW'
                        else:
                            Native_loss = 'IMPROVED'
            
            #Annotate acceptor splice site    
            elif Splicing_var_type == 'Exonic_acceptor_splice_region' or Splicing_var_type == 'Intronic_acceptor_splice_region':    
                if re.search('SNV', VARIANT_CLASS):
                    if MaxEntScan_alt  == '-' or MaxEntScan_diff == '-':
                        Native_loss = '-'
                    else:
                        MaxEntScan_alt = float(MaxEntScan_alt)
                        MaxEntScan_diff = float(MaxEntScan_diff)
                        if MaxEntScan_diff > 0:
                            if MaxEntScan_alt < 6.2:
                                if MaxEntScan_diff >= 1.15:
                                    Native_loss = 'HIGH'
                                else:
                                    Native_loss = 'MODERATE*'
                            elif MaxEntScan_alt > 8.5:
                                Native_loss = 'LOW'
                            else:
                                if MaxEntScan_diff >= 1.15:
                                    Native_loss = 'MODERATE'
                                else:
                                    Native_loss = 'LOW'
                        else:
                            Native_loss = 'IMPROVED'
                else:
                    if MES_SWA_acceptor_diff == '-' or MES_SWA_acceptor_alt =='-':
                        Native_loss = '-'
                    else:
                        MES_SWA_acceptor_alt = float(MES_SWA_acceptor_alt)
                        MES_SWA_acceptor_ref = float(MES_SWA_acceptor_ref)
                        MES_SWA_acceptor_ref_comp = float(MES_SWA_acceptor_ref_comp)
                        MES_SWA_acceptor_diff = float(MES_SWA_acceptor_diff)
  
                        if MES_SWA_acceptor_diff > 0 :
                            if MES_SWA_acceptor_alt < 6.2 :
                                if MES_SWA_acceptor_diff >= 1.15:
                                    Native_loss = 'HIGH'
                                else:
                                    Native_loss = 'MODERATE*'
                            elif MES_SWA_acceptor_alt > 8.5:
                                Native_loss = 'LOW'
                            else:
                                if MES_SWA_acceptor_diff >= 1.15:
                                    Native_loss = 'MODERATE'
                                else:
                                    Native_loss = 'LOW'
                        else:
                            Native_loss = 'IMPROVED'
                   
            else:
                Native_loss = '-'

            #Annoate donor gain within native splice
            if Splicing_var_type == 'Intronic_donor_splice_region' or Splicing_var_type == 'Exonic_donor_splice_region' or Splicing_var_type == 'Exonic_acceptor_splice_region':
                if not re.search('SNV', VARIANT_CLASS):
                    if MES_SWA_donor_alt == '-':
                        Donor_gain_native_site = '-'
                    else:
                        MES_SWA_donor_alt = float(MES_SWA_donor_alt)                        
                        
                        if MES_SWA_donor_alt > 8.5 :
                            Donor_gain_native_site = 'HIGH'
                        elif MES_SWA_donor_alt < 6.2:
                            Donor_gain_native_site = 'LOW'
                        else:
                            Donor_gain_native_site = 'MODERATE'
                else:
                    if Splicing_var_type == 'Intronic_donor_splice_region':
                        if MES_SWA_donor_ref_comp == '-' or MES_NCSS_upstream_donor == '-' or MES_SWA_donor_alt == '-':
                            Donor_gain_native_site = '-'
                        else:
                            MES_SWA_donor_ref_comp = float(MES_SWA_donor_ref_comp)
                            MES_NCSS_upstream_donor = float(MES_NCSS_upstream_donor)
                            MES_SWA_donor_alt = float(MES_SWA_donor_alt)
                           
                            if MES_SWA_donor_ref_comp != MES_NCSS_upstream_donor:
                           # if MES_SWA_donor_diff < 0:
                                if MES_SWA_donor_alt > 8.5 :
                                    Donor_gain_native_site = 'HIGH'
                                elif MES_SWA_donor_alt < 6.2:
                                    Donor_gain_native_site = 'LOW'
                                else:      
                                    Donor_gain_native_site = 'MODERATE'                             
                            else:
                                Donor_gain_native_site = 'LOW'

                    elif Splicing_var_type == 'Exonic_donor_splice_region':
                        if MES_SWA_donor_ref_comp == '-' or MES_NCSS_downstream_donor == '-' or MES_SWA_donor_alt == '-':
                            Donor_gain_native_site = '-'
                        else:
                            MES_SWA_donor_ref_comp = float(MES_SWA_donor_ref_comp)
                            MES_NCSS_downstream_donor = float(MES_NCSS_downstream_donor)
                            MES_SWA_donor_alt = float(MES_SWA_donor_alt)
                            
                            if MES_SWA_donor_ref_comp != MES_NCSS_downstream_donor :
                            #if MES_SWA_donor_diff < 0:
                                if MES_SWA_donor_alt > 8.5 :
                                    Donor_gain_native_site = 'HIGH'
                                elif MES_SWA_donor_alt < 6.2:
                                    Donor_gain_native_site = 'LOW'
                                else:
                                    Donor_gain_native_site = 'MODERATE'
                            else:
                                Donor_gain_native_site = 'LOW'

                    else:
                        if MES_SWA_donor_ref_comp == '-' or MES_NCSS_downstream_donor == '-' or MES_NCSS_upstream_donor == '-' or MES_SWA_donor_alt == '-':
                            Donor_gain_native_site = '-'
                        else:
                            MES_SWA_donor_ref_comp = float(MES_SWA_donor_ref_comp)
                            MES_NCSS_downstream_donor = float(MES_NCSS_downstream_donor)
                            MES_NCSS_upstream_donor = float(MES_NCSS_upstream_donor)
                            MES_SWA_donor_alt = float(MES_SWA_donor_alt)
                            MES_SWA_donor_diff = float(MES_SWA_donor_diff)   
                            if MES_SWA_donor_diff < 0:
                                if MES_SWA_donor_alt > 8.5 :
                                    Donor_gain_native_site = 'HIGH'
                                elif MES_SWA_donor_alt < 6.2:
                                    Donor_gain_native_site = 'LOW'
                                else:
                                    if MES_SWA_donor_alt > MES_NCSS_downstream_donor and  MES_SWA_donor_alt > MES_NCSS_upstream_donor:
                                        Donor_gain_native_site = 'MODERATE[both]'
                                    elif MES_SWA_donor_alt > MES_NCSS_downstream_donor:
                                        Donor_gain_native_site = 'MODERATE[downstream]'
                                    elif MES_SWA_donor_alt > MES_NCSS_upstream_donor:
                                        Donor_gain_native_site = 'MODERATE[upstream]'
                                    else:
                                        Donor_gain_native_site = 'LOW'
                            else:
                                Donor_gain_native_site = 'LOW'
      
            else:
                Donor_gain_native_site = '-'
                         
            #Annotate donor gain outside native site
            if Splicing_var_type == 'Outside_native':
                if MES_NCSS_downstream_donor == '-' or MES_SWA_donor_diff == '-' or MES_NCSS_upstream_donor == '-' or MES_SWA_donor_alt == '-':
                    Donor_gain = '-'
                else:
                    MES_NCSS_downstream_donor = float(MES_NCSS_downstream_donor)
                    MES_SWA_donor_diff = float(MES_SWA_donor_diff)
                    MES_SWA_donor_alt = float(MES_SWA_donor_alt)
                    MES_NCSS_upstream_donor = float(MES_NCSS_upstream_donor)
    
                    if MES_SWA_donor_diff < 0:
                        if MES_SWA_donor_alt > 8.5 :
                            Donor_gain = 'HIGH'
                        elif MES_SWA_donor_alt < 6.2:
                            Donor_gain = 'LOW'
                        else:
                            if re.search('intron',Consequence): 
                                if MES_SWA_donor_alt > MES_NCSS_upstream_donor:
                                    Donor_gain = 'MODERATE[intronic]'
                                else:
                                    Donor_gain = 'LOW'                                        
                            else:
                                if MES_SWA_donor_alt > MES_NCSS_downstream_donor:
                                    Donor_gain = 'MODERATE'    
                                else:
                                    Donor_gain = 'LOW'
                    else:
                        Donor_gain = 'LOW'
    
            else:
                Donor_gain = '-'    

            outfile.write('\t'.join(row + [Splicing_var_type,Native_loss,Donor_gain_native_site,Donor_gain]) + '\n')
 


	
if __name__ == '__main__':
    main()
