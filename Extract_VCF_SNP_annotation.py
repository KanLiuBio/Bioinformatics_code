#!/usr/bin/python

Usage = """
Extract annotation results from SNPEff result of vcf format.

Usage: -version 1.0 (Python3)
  Extract_VCF_SNP_annotation.py input_vcf_file output_annotation_txt

Kan Liu
liukan.big@gmail.com
"""

from Bio import SeqIO
import sys, re

if len(sys.argv)<3:
    print (Usage)
else:
    with open(sys.argv[2], 'w', newline='\n') as outfile, open(sys.argv[1], 'r', encoding='utf-8') as infile:
        outfile.write('CHROM\tPOS\tREF\tALT\tAllDepth\tInividualDepth\tAnnotation\n')
        for line in infile:
            if re.match("\#", line):
                continue
            arrayLine=re.split('\t',line.strip())
            new_list = list()
            if len(arrayLine) < 4:
                print("No record of this line:" + line.strip())
                continue
            new_list.append(arrayLine[0]) 
            new_list.append(arrayLine[1]) 
            new_list.append(arrayLine[3]) 
            new_list.append(arrayLine[4])
            total_depth = re.match("DP\=\d+", arrayLine[7])
            new_list.append(total_depth[0]) 
            individual_depth = re.search("DP4\=\d+\,\d+\,\d+\,\d+", arrayLine[7])
            new_list.append(individual_depth[0])
            temp1 = [num for num in re.split('\;',arrayLine[7]) if num[0:3] == 'ANN']
            annotaion = re.split('\|',temp1[0])
            new_list.extend(annotaion)
            outfile.write("\t".join(new_list)+"\n")
            

