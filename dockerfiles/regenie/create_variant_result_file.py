################################################
#   Libraries
################################################

import click
from granite.lib import vcf_parser
import math
from scipy.stats import fisher_exact
from utils import get_worst_consequence, get_worst_transcript


################################################
#   Top level variables
################################################

# AF is reported as 3 siginificant digits
# from Sentieon/GATK so we will use that for
# AF and p-values

# significant digits when calculated below 1
# e.g., AF and p
SIGNIFICANT_DIGITS = 3

#significant digits when calculated above 1
# e.g., OR and log10
ROUND_DIGITS = 4


################################################
#   Top-level variables
################################################

VEP_TAG = 'CSQ'

################################################
#   Functions
################################################


def round_to_significant_digits(v):
    return round(v, SIGNIFICANT_DIGITS - int(math.floor(math.log10(abs(v)))) - 1)

def summarize_genotypes(sample_gts, variant_id):
    ''' 
    Expects a dict of the form
    {
        "sample_1": "0/1",
        "sample_2": "./1",
    }
    and calculates AC, AN and AF across all samples.
    
    '''

    VALID_GENOTYPES = ["./.", "0/0", "1/0", "0/1", "1/1" , "0|0", "1|0", "0|1", "1|1"]

    s_AC = s_AN = 0

    for sample in sample_gts:
        gt = sample_gts[sample]

        if gt not in VALID_GENOTYPES:
            raise Exception(f"Unexpected genotype {gt} found for variant {variant_id}. Did you run bcftools norm multiallelics?")

        # if genotye has not been called, move on
        if gt == "./.":
            continue

        s_AN += 2 # add 2 to total count if there is a genotype called
        s_AC += gt.count('1') # then count the alternative alleles
        
    # with counting done, we need to calculate AF and store all 3 values in the result dict
    s_AF = round_to_significant_digits(s_AC/s_AN) if s_AC > 0 else 0

    return {
        "AC": s_AC,
        "AN": s_AN,
        "AF": s_AF,
    }

def fisher_calculation(proband_alt, proband_ref, gnomAD_alt, gnomAD_ref):
    '''
    This function is called within fisher_exact_gnomAD
    and runs the actual Fisher exact test once the values
    have been generated by the main function.
    '''
    oddsratio, pvalue = fisher_exact([[proband_alt, proband_ref], [gnomAD_alt, gnomAD_ref]], alternative='greater')
    p_value_rounded = round_to_significant_digits(pvalue)
    oddsratio_rounded = round(oddsratio, ROUND_DIGITS)
    minuslog10p = 0 if pvalue == 1 else round(-math.log10(pvalue), ROUND_DIGITS)
    return p_value_rounded, oddsratio_rounded, minuslog10p

def fisher_exact_gnomAD(case_AC, case_AN, gnomAD_AC, gnomAD_AN):
    '''
    This generated the necessary inputs for a  Fisher
    exact test (2 by 2) on the probands (cases)
    versus gnomAD values
    and calls the fisher_calculation function to
    carry out the test. It return a tuple: p_value, oddsratio, minuslog10p
    '''

    #can only run the test if there are values for gnomAD
    if gnomAD_AC != '' and gnomAD_AN != '':
        #store the proband info
        proband_alt = int(case_AC)
        proband_ref = int(case_AN) - proband_alt

        #then access the gnomAD info
        #gnomAD v2 has some entries with multiple values so we need to condition on their absence first
        if "&" not in gnomAD_AC and "&" not in gnomAD_AN:
            gnomAD_alt = int(gnomAD_AC)
            gnomAD_ref = int(gnomAD_AN) - gnomAD_alt
            return fisher_calculation(proband_alt, proband_ref, gnomAD_alt, gnomAD_ref)

        # if the entry does have an "&" we want to select the most rare to compare to
        else:
            # found that some variants with "&" in AC and AN do not have "&" in AF, which means we can't use AF index to pull AC and AN.
            # need to carry out the calculation
            gnomAD_alt_list = [int(x) for x in gnomAD_AC.split("&")]
            gnomAD_AN_list = [int(x) for x in gnomAD_AN.split("&")]

            # if the lists are of different lengths, we can't use a shared index, so we have to ignore variant
            # We want to select the most rare from the list of possibilites. If there is a 0 in the AN list, we
            # can't compute that. Skip variant in that case
            if (len(gnomAD_alt_list) != len(gnomAD_AN_list)) or (0 in gnomAD_AN_list):
                return 'NA', 'NA', 'NA'

            allele_frequencies = [i / j for i, j in zip(gnomAD_alt_list, gnomAD_AN_list)]
            index_min_af = allele_frequencies.index(min(allele_frequencies))

            gnomAD_alt = gnomAD_alt_list[index_min_af]
            gnomAD_ref = gnomAD_AN_list[index_min_af] - gnomAD_alt_list[index_min_af]

            return fisher_calculation(proband_alt, proband_ref, gnomAD_alt, gnomAD_ref)

    return 'NA', 'NA', 'NA'


@click.command()
@click.help_option("--help", "-h")
@click.option("-r", "--regenie-output", required=True, type=str, help="Regenie output file")
@click.option("-a", "--annotated-vcf", required=True, type=str, help="Annotated, jointly called VCF")
@click.option("-c", "--cases", required=True, type=str, help="List of case sample IDs (comma separated)")
@click.option("-o", "--out", required=True, type=str, help="the output file name of the variant level results")
@click.option("-e", "--higlass-vcf", required=True, type=str, help="Higlass VCF file containing the results")
def main(regenie_output, annotated_vcf, cases, out, higlass_vcf):
    """This script takes a variant-based regenie output file and adds Fisher exact test results.
       It also produces a Higlass compatible VCF with some annotations

    Example usage: 

    python create_variant_result_file.py -r /path/to/out.regenie -a regenie_input_source.vcf -o variant_level_results.txt -e higlass_variant_tests.vcf

    """

    vcf_obj = vcf_parser.Vcf(annotated_vcf)
    idx_consequence = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'Consequence')
    idx_canonical = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'CANONICAL')
    idx_impact = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'IMPACT')
    idx_cadd = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'CADD_PHRED')

    # We want gnomAD v2 and v3 allele frequencies etc. to computer Fisher exact test scores
    idx_gnomADg_AC = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'gnomADg_AC')
    idx_gnomADg_AN = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'gnomADg_AN')
    idx_gnomADg_AF = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'gnomADg_AF')
    idx_gnomADe2_AC = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'gnomADe2_AC')
    idx_gnomADe2_AN = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'gnomADe2_AN')
    idx_gnomADe2_AF = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'gnomADe2_AF')

    cohort_sample_ids = vcf_obj.header.IDs_genotypes # This includes cases and controls
    case_sample_ids = cases.split(",")
    control_sample_ids = [id for id in cohort_sample_ids if id not in case_sample_ids]

    # Verify that every case ID is present in the cohort VCF
    if(not set(case_sample_ids).issubset(set(cohort_sample_ids))):
        raise Exception("Not every case ID could be found in the cohort VCF.")
        
    variant_infos = {}

    for record in vcf_obj.parse_variants():
        id = record.ID
        # Retrieve annotations and allele counts
        try: 
            vep_tag_value = record.get_tag_value(VEP_TAG)
            worst_transcript = get_worst_transcript(vep_tag_value, idx_canonical, idx_consequence)
            worst_transcript_ = worst_transcript.split('|')
            worst_consequence = get_worst_consequence(worst_transcript_[idx_consequence])
            impact = worst_transcript_[idx_impact]
            cadd = worst_transcript_[idx_cadd]
            gnomADg_AC = worst_transcript_[idx_gnomADg_AC]
            gnomADg_AN = worst_transcript_[idx_gnomADg_AN]
            gnomADg_AF = worst_transcript_[idx_gnomADg_AF]
            gnomADe2_AC = worst_transcript_[idx_gnomADe2_AC]
            gnomADe2_AN = worst_transcript_[idx_gnomADe2_AN]
            gnomADe2_AF = worst_transcript_[idx_gnomADe2_AF]

            # get the index for genotype (GT) and pull genotypes for all samples
            GT_idx = record.FORMAT.split(":").index("GT")
            case_sample_genotypes = {}
            for sample in case_sample_ids:
                case_sample_genotypes[sample] = record.GENOTYPES[sample].split(":")[GT_idx]

            control_sample_genotypes = {}
            for sample in control_sample_ids:
                control_sample_genotypes[sample] = record.GENOTYPES[sample].split(":")[GT_idx]

            case_sample_gt_summarized = summarize_genotypes(case_sample_genotypes, id)
            control_sample_gt_summarized = summarize_genotypes(control_sample_genotypes, id)

            case_AC = case_sample_gt_summarized["AC"]
            case_AN = case_sample_gt_summarized["AN"]
            control_AC = control_sample_gt_summarized["AC"]
            control_AN = control_sample_gt_summarized["AN"]

            # Perform Fisher calculations for the different control groups
            fisher_p_gnomADg, fisher_or_gnomADg, fisher_ml10p_gnomADg = fisher_exact_gnomAD(case_AC, case_AN, gnomADg_AC, gnomADg_AN)
            fisher_p_gnomADe2, fisher_or_gnomADe2, fisher_ml10p_gnomADe2 = fisher_exact_gnomAD(case_AC, case_AN, gnomADe2_AC, gnomADe2_AN)
            fisher_p_control, fisher_or_control, fisher_ml10p_control = fisher_calculation(case_AC, case_AN-case_AC, control_AC, control_AN-control_AC) 

            variant_infos[id] = {
                "most_severe_consequence": worst_consequence,
                "level_most_severe_consequence": impact,
                "cadd_phred": cadd,
                "case_AC": case_AC,
                "case_AN": case_AN,
                "case_N": int(case_AN/2),
                "case_AF": case_sample_gt_summarized["AF"],
                "control_AC": control_AC,
                "control_AN": control_AN,
                "control_N": int(control_AN/2),
                "control_AF": control_sample_gt_summarized["AF"],
                "gnomADg_AC": gnomADg_AC,
                "gnomADg_AN": gnomADg_AN,
                "gnomADg_AF": gnomADg_AF,
                "gnomADe2_AC": gnomADe2_AC,
                "gnomADe2_AN": gnomADe2_AN,
                "gnomADe2_AF": gnomADe2_AF,
                "fisher_p_gnomADg": fisher_p_gnomADg,
                "fisher_or_gnomADg": fisher_or_gnomADg,
                "fisher_ml10p_gnomADg": fisher_ml10p_gnomADg,
                "fisher_p_gnomADe2": fisher_p_gnomADe2,
                "fisher_or_gnomADe2": fisher_or_gnomADe2,
                "fisher_ml10p_gnomADe2": fisher_ml10p_gnomADe2,
                "fisher_p_control": fisher_p_control,
                "fisher_or_control": fisher_or_control,
                "fisher_ml10p_control": fisher_ml10p_control,

            }
        except Exception: 
            raise ValueError('\nERROR processing variant_infos for variant {0}\n'.format(id))


    # Regenie output file format
    # Example line
    # CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
    # 1 13613 chr1_13613_T_A T A 0.0514706 1 68 ADD 0.191978 0.695341 0.0762264 0.106528 NA
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

    info_list = [
        "case_AC", "case_AN", "case_AF", "control_AC", "control_AN", "control_AF", "gnomADg_AC", "gnomADg_AN", "gnomADg_AF", 
        "gnomADe2_AC", "gnomADe2_AN", "gnomADe2_AF", "most_severe_consequence", "level_most_severe_consequence", "cadd_phred",
        "fisher_or_gnomADg", "fisher_ml10p_gnomADg", "fisher_or_gnomADe2", "fisher_ml10p_gnomADe2", "fisher_or_control", 
        "fisher_ml10p_control", "regenie_ml10p", "regenie_beta", "regenie_chisq", "regenie_se",
    ]

    
    with open(regenie_output, 'r') as f_in, open(out, 'w') as f_out, open(higlass_vcf, 'w') as f_out_hg:
        header = '# CHROM: chromosome\n'
        header += '# GENPOS: position with in the chromosome\n'
        header += '# ID: variant ID\n'
        header += '# ALLELE0: reference allele\n'
        header += '# ALLELE1: alternative allele\n'
        header += '# A1FREQ: frequency of the alternative allele\n'
        header += '# N: sample size (N = CASE_N + CONTROL_N)\n'
        header += '# R_TEST: test performed by Regenie (additive/dominant/recessive)\n'
        header += '# R_BETA: estimated effect sizes (Regenie)\n'
        header += '# R_SE: standard error of the Regenie test\n'
        header += '# R_CHISQ: chi-square test statistics of the Regenie test\n'
        header += '# R_LOG10P: -log10(p) of the Regenie test\n'
        header += '# CASE_AF: case allele frequency\n'
        header += '# CASE_N: number of affected samples\n'
        header += '# CONTROL_AF: control allele frequency\n'
        header += '# CONTROL_N: number of control samples\n'
        header += '# F_LOG10P_GNOMADG: -log10(p) of a Fisher exact test when gnomAD 3 is used as control group\n'
        header += '# F_OR_GNOMADG: odds ratio of a Fisher exact test when gnomAD 3 is used as control group\n'
        header += '# F_LOG10P_GNOMADE2: -log10(p) of a Fisher exact test when gnomAD 2 is used as control group\n'
        header += '# F_OR_GNOMADE2: odds ratio of a Fisher exact test when gnomAD 2 is used as control group\n'
        header += 'CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N R_TEST R_BETA R_SE R_CHISQ R_LOG10P CASE_AF CASE_N CONTROL_AF CONTROL_N F_LOG10P_GNOMADG F_OR_GNOMADG F_LOG10P_GNOMADE2 F_OR_GNOMADE2\n'
        f_out.write(header)

        f_out_hg.write('##fileformat=VCFv4.3\n')
        f_out_hg.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        
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
            freq = result[r_map["A1FREQ"]]
            samples = result[r_map["N"]]
            r_test = result[r_map["TEST"]]
            r_beta = result[r_map["BETA"]]
            r_se = result[r_map["SE"]]
            r_chisq = result[r_map["CHISQ"]]
            r_log10p = result[r_map["LOG10P"]]

            variant_infos[id]["regenie_ml10p"] = r_log10p
            variant_infos[id]["regenie_beta"] = r_beta
            variant_infos[id]["regenie_chisq"] = r_chisq
            variant_infos[id]["regenie_se"] = r_se
            vi = variant_infos[id]
    
            l = f"{chr} {pos} {id} {ref} {alt} {freq} {samples} {r_test} {r_beta} {r_se} {r_chisq} {r_log10p} {vi['case_AF']} {vi['case_N']} {vi['control_AF']} {vi['control_N']} {vi['fisher_ml10p_gnomADg']} {vi['fisher_or_gnomADg']} {vi['fisher_ml10p_gnomADe2']} {vi['fisher_or_gnomADe2']}"
            f_out.write(f"{l}\n")
            

            info = ""
            for field in info_list:
                info+=field+"="
                if vi[field] == '':
                    info+="NA;"
                else:
                    info+=str(vi[field])+";"
            info = info.strip(";")
            
            f_out_hg.write(f"{chr}\t{pos}\t{id}\t{ref}\t{alt}\t0\tPASS\t{info}\n")
  


if __name__ == "__main__":
    main()