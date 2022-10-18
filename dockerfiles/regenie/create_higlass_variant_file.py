import click
import json
from granite.lib import vcf_parser

@click.command()
@click.help_option("--help", "-h")
@click.option("-r", "--regenie-output", required=True, type=str, help="Regenie output file")
@click.option("-a", "--annotated-vcf", required=True, type=str, help="Annotated, jointly called VCF")
@click.option("-o", "--out", required=True, type=str, help="the output file name")
def main(regenie_output, annotated_vcf, out):
    """This script takes a variant-based regenie output file and transforms it in a vcf 
       that the Higlass browser can understand

    Example usage: 

    python create_higlass_variant_file.py -r /path/to/out.regenie -a regenie_input_source.vcf -o higlass_variant_tests.vcf

    """

    vcf_obj = vcf_parser.Vcf(annotated_vcf)
    idx_consequence = vcf_obj.header.get_tag_field_idx('CSQ', 'Consequence')
    idx_impact = vcf_obj.header.get_tag_field_idx('CSQ', 'IMPACT')
    idx_cadd = vcf_obj.header.get_tag_field_idx('CSQ', 'CADD_PHRED')

    annotations = {}

    for record in vcf_obj.parse_variants():
        info_split = record.INFO.split(";")
        for field in info_split:
            if "CSQ=" in field:
                csq_split = field.split("|")
                id = record.ID
                consequence = csq_split[idx_consequence]
                impact = csq_split[idx_impact]
                cadd = csq_split[idx_cadd]

                annotations[id] = {
                    "consequence": consequence,
                    "impact": impact,
                    "cadd": cadd
                }


    r_map = {
        "CHROM": 0,
        "GENPOS": 1,
        "ID": 2,
        "ALLELE0": 3,
        "ALLELE1": 4,
        "A1FREQ": 5,
        "INFO": 6,
        "N": 7,
        "TEST": 8,
        "BETA": 9,
        "SE": 10,
        "CHISQ": 11,
        "LOG10P": 12,
        "EXTRA": 13,

    }

    
    with open(regenie_output, 'r') as f_in, open(out, 'w') as f_out:

        f_out.write('##fileformat=VCFv4.3\n')
        f_out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        
        for line in f_in:
            if line.startswith("##") or line.startswith("CHROM"):
                continue

            # Example line
            # CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
            # 1 13613 chr1_13613_T_A T A 0.0514706 1 68 ADD 0.191978 0.695341 0.0762264 0.106528 NA

            result = line.strip().split()
            id = result[r_map["ID"]]
            chr = id.split("_")[0]
            pos = result[r_map["GENPOS"]]
            ref = result[r_map["ALLELE0"]]
            alt = result[r_map["ALLELE1"]]
            alt_allele_freq = result[r_map["A1FREQ"]]
            sample_size = result[r_map["N"]]
            p_value = result[r_map["LOG10P"]]
            info = f"regenie_log10p={p_value};"
            consequence = annotations[id]["consequence"]
            info += f"most_severe_consequence={consequence};"
            impact = annotations[id]["impact"]
            info += f"level_most_severe_consequence={impact};"
            cadd = annotations[id]["cadd"]
            info += f"cadd_scaled={cadd};"
            info += f"alt_allele_freq={alt_allele_freq};"
            info += f"sample_size={sample_size}"
            

            f_out.write(f"{chr}\t{pos}\t{id}\t{ref}\t{alt}\t0\tPASS\t{info}\n")


if __name__ == "__main__":
    main()