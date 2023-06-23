################################################
#   Libraries
################################################

import click
from granite.lib import vcf_parser
from granite.lib.shared_vars import DStags
import math
from scipy.stats import fisher_exact
from utils import get_worst_consequence, get_worst_transcript, clean_dbnsfp, get_maxds
from utils import get_cases, VALID_GENOTYPES

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

dbNSFP_fields = {
    # dbNSFP fields that may be a list
    # and need to be assigned to transcripts
    'Polyphen2_HVAR_pred': 0,
    'Polyphen2_HVAR_score': 0,
    'SIFT_pred': 0,
    'SIFT_score': 0
}


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
@click.option("-s", "--sample-info", required=True, type=str, help="JSON string with sample info")
@click.option("-o", "--out", required=True, type=str, help="the output file name of the variant level results")
@click.option("-f", "--af-threshold-higlass", required=True, type=str, help="Rare variant AF threshold for variants to include in Higlass")
@click.option("-e", "--higlass-vcf", required=True, type=str, help="Higlass VCF file containing the results")
def main(regenie_output, annotated_vcf, sample_info, out, af_threshold_higlass, higlass_vcf):
    """This script takes a variant-based regenie output file and adds Fisher exact test results.
       It also produces a Higlass compatible VCF with some annotations

    Example usage: 

    python create_variant_result_file.py -r /path/to/out.regenie -a regenie_input_source.vcf -o variant_level_results.txt -e higlass_variant_tests.vcf

    """

    vcf_obj = vcf_parser.Vcf(annotated_vcf)

    # Get annotation indices
    idx_transcript = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'Ensembl_transcriptid')
    idx_consequence = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'Consequence')
    idx_canonical = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'CANONICAL')
    idx_impact = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'IMPACT')
    idx_cadd_phred = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'CADD_PHRED') 
    idx_cadd_rs = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'CADD_raw_rankscore')
    idx_enst = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'Feature')

    idx_polyphen_pred = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'Polyphen2_HVAR_pred')
    idx_polyphen_rankscore = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'Polyphen2_HVAR_rankscore')
    idx_polyphen_score = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'Polyphen2_HVAR_score')

    idx_gerp_score = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'GERP++_RS')
    idx_gerp_rankscore = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'GERP++_RS_rankscore')

    idx_sift_rankscore = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'SIFT_converted_rankscore')
    idx_sift_pred = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'SIFT_pred')
    idx_sift_score = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'SIFT_score')

    # idx_spliceai_ds_ag = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'SpliceAI_pred_DS_AG')
    # idx_spliceai_ds_al = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'SpliceAI_pred_DS_AL')
    # idx_spliceai_ds_dg = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'SpliceAI_pred_DS_DG')
    # idx_spliceai_ds_dl = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'SpliceAI_pred_DS_DL')
    # Get SpliceAI ds indexes
    # DStags import from granite.shared_vars
    SpAItag_list, SpAI_idx_list = [], []
    for DStag in DStags:
        tag, idx = vcf_obj.header.check_tag_definition(DStag)
        SpAItag_list.append(tag)
        SpAI_idx_list.append(idx)
    #end for


    # We want gnomAD v2 and v3 allele frequencies etc. to computer Fisher exact test scores
    idx_gnomADg_AC = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'gnomADg_AC')
    idx_gnomADg_AN = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'gnomADg_AN')
    idx_gnomADg_AF = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'gnomADg_AF')
    idx_gnomADe2_AC = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'gnomADe2_AC')
    idx_gnomADe2_AN = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'gnomADe2_AN')
    idx_gnomADe2_AF = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'gnomADe2_AF')

    cohort_sample_ids = vcf_obj.header.IDs_genotypes # This includes cases and controls
    case_sample_ids =  get_cases(sample_info)
    control_sample_ids = [id for id in cohort_sample_ids if id not in case_sample_ids]

    # Verify that every case ID is present in the cohort VCF
    if(not set(case_sample_ids).issubset(set(cohort_sample_ids))):
        raise Exception("Not every case ID could be found in the cohort VCF.")
        
    variant_infos = {}
    variant_ids = []

    # Get VEP indexes
    # Indexes to resolve dbNSFP values by transcript
    for field in dbNSFP_fields:
        dbNSFP_fields[field] = vcf_obj.header.get_tag_field_idx(VEP_TAG, field)
    #end for
    

    for record in vcf_obj.parse_variants():
        id = record.ID
        variant_ids.append(id)
        # Retrieve annotations and allele counts
        try: 
            #vep_tag_value = record.get_tag_value(VEP_TAG)

             # Clean dbNSFP by resolving values by transcript
            VEP_clean = clean_dbnsfp(record, VEP_TAG, dbNSFP_fields, idx_transcript, idx_enst)
            if not VEP_clean: continue

            worst_transcript = get_worst_transcript(VEP_clean, idx_canonical, idx_consequence)
            worst_transcript_ = worst_transcript.split('|')
            worst_consequence = get_worst_consequence(worst_transcript_[idx_consequence])
            impact = worst_transcript_[idx_impact]
            transcript_id = worst_transcript_[idx_enst]

            cadd_phred = worst_transcript_[idx_cadd_phred]
            cadd_raw_rs = worst_transcript_[idx_cadd_rs]

            polyphen_pred = worst_transcript_[idx_polyphen_pred]
            polyphen_rankscore = worst_transcript_[idx_polyphen_rankscore]
            polyphen_score = worst_transcript_[idx_polyphen_score]

            gerp_score = worst_transcript_[idx_gerp_score]
            gerp_rankscore = worst_transcript_[idx_gerp_rankscore]

            sift_rankscore = worst_transcript_[idx_sift_rankscore]
            sift_pred = worst_transcript_[idx_sift_pred]
            sift_score = worst_transcript_[idx_sift_score]

            # spliceai_ds_ag = worst_transcript_[idx_spliceai_ds_ag]
            # spliceai_ds_al = worst_transcript_[idx_spliceai_ds_al]
            # spliceai_ds_dg = worst_transcript_[idx_spliceai_ds_dg]
            # spliceai_ds_dl = worst_transcript_[idx_spliceai_ds_dl]
            # spliceai_score_max = max([spliceai_ds_ag,spliceai_ds_al,spliceai_ds_dg,spliceai_ds_dl])
            # Get max SpliceAI max_ds
            spliceai_score_max = get_maxds(record, SpAItag_list, SpAI_idx_list)
            spliceai_score_max = spliceai_score_max if spliceai_score_max else ''

            
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

            # Include Higlass specific filtering into this logic
            include_for_higlass = True
            if gnomADg_AF and float(gnomADg_AF) > float(af_threshold_higlass):
                include_for_higlass = False

            variant_infos[id] = {
                # We are adding these general infos here as well, so that we don't have to run parse_variants again later
                "chrom": record.CHROM,
                "pos": record.POS,
                "ref": record.REF,
                "alt": record.ALT,

                "transcript": transcript_id,
                "most_severe_consequence": worst_consequence,
                "level_most_severe_consequence": impact,
                "case_AC": case_AC,
                "case_AN": case_AN,
                "case_N": int(case_AN/2),
                "case_AF": case_sample_gt_summarized["AF"],
                "control_AC": control_AC,
                "control_AN": control_AN,
                "control_N": int(control_AN/2),
                "control_AF": control_sample_gt_summarized["AF"],  
                                                                                            
                "cadd_raw_rs": cadd_raw_rs,
                "cadd_phred": cadd_phred,
                "polyphen_pred": polyphen_pred,
                "polyphen_rankscore": polyphen_rankscore,
                "polyphen_score": polyphen_score,
                "gerp_score": gerp_score,
                "gerp_rankscore": gerp_rankscore,
                "sift_rankscore": sift_rankscore,
                "sift_pred": sift_pred,
                "sift_score": sift_score,
                "spliceai_score_max": spliceai_score_max,

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

                # Placeholder for Regenie results
                "regenie_ml10p" : '',
                "regenie_beta" : '',
                "regenie_chisq" : '',
                "regenie_se" : '',
                "regenie_test": '',

                "include_for_higlass": include_for_higlass
            }


            for key in variant_infos[id]:
                if variant_infos[id][key] == '':
                    variant_infos[id][key] = 'NA'

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

    # Add Regenie results to the variant_infos
    with open(regenie_output, 'r') as f_in:
       
        for line in f_in:
            if line.startswith("##") or line.startswith("CHROM"):
                continue

            # Example line
            # CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
            # 1 13613 chr1_13613_T_A T A 0.0514706 1 68 ADD 0.191978 0.695341 0.0762264 0.106528 NA

            result = line.strip().split()
            id = result[r_map["ID"]]
            r_beta = result[r_map["BETA"]]
            r_se = result[r_map["SE"]]
            r_chisq = result[r_map["CHISQ"]]
            r_log10p = result[r_map["LOG10P"]]
            r_test = result[r_map["TEST"]]

            variant_infos[id]["regenie_ml10p"] = r_log10p
            variant_infos[id]["regenie_beta"] = r_beta
            variant_infos[id]["regenie_chisq"] = r_chisq
            variant_infos[id]["regenie_se"] = r_se
            variant_infos[id]["regenie_test"] = r_test

           
    info_list = [
         "transcript", "case_AC", "case_AN", "case_AF", "control_AC", "control_AN", "control_AF", "gnomADg_AC", "gnomADg_AN", "gnomADg_AF", "gnomADe2_AC", "gnomADe2_AN", "gnomADe2_AF", "most_severe_consequence", "level_most_severe_consequence", "cadd_raw_rs", "cadd_phred", "polyphen_pred", "polyphen_rankscore", "polyphen_score", "gerp_score", "gerp_rankscore", "sift_rankscore", "sift_pred", "sift_score", "spliceai_score_max", "fisher_or_gnomADg", "fisher_ml10p_gnomADg", "fisher_or_gnomADe2", "fisher_ml10p_gnomADe2", "fisher_or_control", "fisher_ml10p_control", "regenie_ml10p", "regenie_beta", "regenie_chisq", "regenie_se",
    ]

    with open(out, 'w') as f_out, open(higlass_vcf, 'w') as f_out_hg:
        header = '# CHROM: chromosome\n'
        header += '# GENPOS: position with in the chromosome\n'
        header += '# ID: variant ID\n'
        header += '# ALLELE0: reference allele\n'
        header += '# ALLELE1: alternative allele\n'
        header += '# R_TEST: test performed by Regenie (additive/dominant/recessive)\n'
        header += '# R_BETA: estimated effect sizes (Regenie)\n'
        header += '# R_SE: standard error of the Regenie test\n'
        header += '# R_CHISQ: chi-square test statistics of the Regenie test\n'
        header += '# R_LOG10P: -log10(p) of the Regenie test\n'
        header += '# CASE_AF: case allele frequency\n'
        header += '# CASE_N: number of affected samples\n'
        header += '# CONTROL_AF: control allele frequency\n'
        header += '# CONTROL_N: number of control samples\n'
        header += '# F_LOG10P_CONTROL: -log10(p) of a Fisher exact test with cases vs. control\n'
        header += '# F_OR_CONTROL: odds ratio of a Fisher exact test with cases vs. control\n'
        header += '# F_LOG10P_GNOMADG: -log10(p) of a Fisher exact test when gnomAD 3 is used as control group\n'
        header += '# F_OR_GNOMADG: odds ratio of a Fisher exact test when gnomAD 3 is used as control group\n'
        header += '# F_LOG10P_GNOMADE2: -log10(p) of a Fisher exact test when gnomAD 2 is used as control group\n'
        header += '# F_OR_GNOMADE2: odds ratio of a Fisher exact test when gnomAD 2 is used as control group\n'
        header += '# CADD_RAW_RS: CADD rankscore\n'
        header += '# CADD_PHRED: CADD Phred score\n'
        header += '# POLYPHEN_PRED: PolyPhen 2 prediction\n'
        header += '# POLYPHEN_RANKSCORE: PolyPhen 2 rankscore\n'
        header += '# POLYPHEN_SCORE: PolyPhen 2 score\n'
        header += '# GERP_SCORE: Gerp++ score\n'
        header += '# GERP_RANKSCORE: Gerp++ rankscore\n'
        header += '# SIFT_RANKSCORE: SIFT rankscore\n'
        header += '# SIFT_PRED: SIFT prediction\n'
        header += '# SIFT_SCORE: SIFT score\n'
        header += '# SPLICEAI_MAX_SCORE: SpliceAI predicts whether a variant causes a splice acceptor gain or loss, or a splice donor gain or loss. The score shown here is the max score of these four scores\n'
        header += 'CHROM GENPOS ID ALLELE0 ALLELE1 R_TEST R_BETA R_SE R_CHISQ R_LOG10P CASE_AF CASE_N CONTROL_AF CONTROL_N F_LOG10P_CONTROL F_OR_CONTROL F_LOG10P_GNOMADG F_OR_GNOMADG F_LOG10P_GNOMADE2 F_OR_GNOMADE2 CADD_RAW_RS CADD_PHRED POLYPHEN_PRED POLYPHEN_RANKSCORE POLYPHEN_SCORE GERP_SCORE GERP_RANKSCORE SIFT_RANKSCORE SIFT_PRED SIFT_SCORE SPLICEAI_MAX_SCORE\n'
        f_out.write(header)

        f_out_hg.write('##fileformat=VCFv4.3\n')
        f_out_hg.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

        for id in variant_ids:
            vi = variant_infos[id]

            l = f"{vi['chrom']} {vi['pos']} {id} {vi['ref']} {vi['alt']} {vi['regenie_test']} {vi['regenie_beta']} {vi['regenie_se']} {vi['regenie_chisq']} {vi['regenie_ml10p']} {vi['case_AF']} {vi['case_N']} {vi['control_AF']} {vi['control_N']} {vi['fisher_ml10p_control']} {vi['fisher_or_control']}  {vi['fisher_ml10p_gnomADg']} {vi['fisher_or_gnomADg']} {vi['fisher_ml10p_gnomADe2']} {vi['fisher_or_gnomADe2']} {vi['cadd_raw_rs']} {vi['cadd_phred']} {vi['polyphen_pred']} {vi['polyphen_rankscore']} {vi['polyphen_score']} {vi['gerp_score']} {vi['gerp_rankscore']} {vi['sift_rankscore']} {vi['sift_pred']} {vi['sift_score']} {vi['spliceai_score_max']}"
            f_out.write(f"{l}\n")
            
            info = ""
            for field in info_list:
                if vi[field] == 'NA':
                    continue
                info+=field+"="+str(vi[field])+";"
            info = info.strip(";")
            
            if vi["include_for_higlass"]:
                f_out_hg.write(f"{vi['chrom']}\t{vi['pos']}\t{id}\t{vi['ref']}\t{vi['alt']}\t0\tPASS\t{info}\n")



if __name__ == "__main__":
    main()