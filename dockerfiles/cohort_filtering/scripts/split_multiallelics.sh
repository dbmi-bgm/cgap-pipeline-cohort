#!/bin/bash

# variables from command line
input_vcf=$1
reference=$2

# run bcftools
bcftools --version
#bcftools norm -m -any -f $reference -o split_tmp.vcf -O v $input_vcf || exit 1
bcftools norm -m -any -f $reference -o split.vcf.gz --threads 6 -O z $input_vcf || exit 1

# py_script="
# fo = open('split.vcf', 'w')
# with open('split_tmp.vcf') as fi:
#     for line in fi:
#         if line.startswith('#'):
#             fo.write(line)
#         else:
#             line_split = line.rstrip().split('\t')
#             if line_split[4] == '*':
#                 line_split[4] = '-'
#             #end if
#             fo.write('\t'.join(line_split) + '\n')
#         #end if
#     #end for
# #end with

# fo.close()
# "

# python -c "$py_script"

# compress and index output vcf
#bgzip split.vcf || exit 1
tabix -p vcf split.vcf.gz || exit 1

echo ""
echo "== DONE =="
