from granite.lib.shared_functions import *
import json
import gzip

VALID_GENOTYPES = ["./.", "0/0", "1/0", "0/1", "1/1" , "0|0", "1|0", "0|1", "1|1"]

VEP_ORDER = {
    # HIGH
    'transcript_ablation': 1,
    'splice_acceptor_variant': 2,
    'splice_donor_variant': 3,
    'stop_gained': 4,
    'frameshift_variant': 5,
    'stop_lost': 6,
    'start_lost': 7,
    'transcript_amplification': 8,
    # MODERATE
    'inframe_insertion': 9,
    'inframe_deletion': 10,
    'missense_variant': 11,
    'protein_altering_variant': 12,
    # LOW
    'splice_region_variant': 13,
    'incomplete_terminal_codon_variant': 14,
    'start_retained_variant': 15,
    'stop_retained_variant': 16,
    'synonymous_variant': 17,
    # MODIFIER
    'coding_sequence_variant': 18,
    'mature_miRNA_variant': 19,
    '5_prime_UTR_variant': 20,
    '3_prime_UTR_variant': 21,
    'intron_variant': 22,
    'MODIFIER': 23
}

def get_cases(sample_info):
    sample_info_dec = json.loads(sample_info)
    sample_info_cases = filter(lambda s: s["is_affected"], sample_info_dec)
    return list(map(lambda x: x["sample_id"], sample_info_cases))

def get_worst_consequence(consequence, sep='&'):
    ''' '''
    consequence_tup = []
    for cnsqce in consequence.split(sep):
        try:
            consequence_tup.append((VEP_ORDER[cnsqce], cnsqce))
        except Exception:
            consequence_tup.append((VEP_ORDER['MODIFIER'], cnsqce))
        #end try
    #end for
    return sorted(consequence_tup, key=lambda x_y: x_y[0])[0][1]
#end def


def get_worst_transcript(VEP_val, CNONICL_idx, CONSEQUENCE_idx, sep='&'):
    ''' '''
    # Check transcripts
    worst_trscrpt_tup = []
    trscrpt_list = VEP_val.split(',')
    # Assign worst impact to transcripts
    for trscrpt in trscrpt_list:
        trscrpt_cnsqce = trscrpt.split('|')[CONSEQUENCE_idx]
        worst_cnsqce = get_worst_consequence(trscrpt_cnsqce)
        try: worst_trscrpt_tup.append((VEP_ORDER[worst_cnsqce], trscrpt))
        except Exception:
            worst_trscrpt_tup.append((VEP_ORDER['MODIFIER'], trscrpt))
        #end try
    #end for
    sorted_worst_trscrpt_tup = sorted(worst_trscrpt_tup, key=lambda x_y: x_y[0])
    worst_impact = sorted_worst_trscrpt_tup[0][0]

    # Get worst transcripts and check canonical
    worst_trscrpt_list = []
    for worst_cnsqce, trscrpt in sorted_worst_trscrpt_tup:
        if worst_cnsqce == worst_impact:
            trscrpt_cnonicl = trscrpt.split('|')[CNONICL_idx]
            if trscrpt_cnonicl == 'YES' or \
                trscrpt_cnonicl == '1':
                return trscrpt
            #end if
            worst_trscrpt_list.append(trscrpt)
        else: break
        #end if
    #end for
    return worst_trscrpt_list[0]
#end def

# Taken from cgap-scripts/portal_reformat_vcf
def get_maxds(vnt_obj, SpAItag_list, SpAI_idx_list):
    ''' '''
    # if SpliceAI is within VEP
    # fetching only the first transcript
    # expected the same scores for all transcripts
    SpAI_vals = []
    for i, SpAItag in enumerate(SpAItag_list):
        SpAI_val = get_tag_idx(vnt_obj, SpAItag, SpAI_idx_list[i])
        # if SpliceAI is with VEP and is at the end of Format
        # need to remove , that separate next transcript
        try: SpAI_vals.append(float(SpAI_val.split(',')[0]))
        except Exception:
            return None
        #end try
    #end for
    if SpAI_vals:
        return max(SpAI_vals)
    #end if
    return None
#end def

# Taken from cgap-scripts/portal_reformat_vcf
def clean_dbnsfp(vnt_obj, VEPtag, dbNSFP_fields, dbnsfp_ENST_idx, ENST_idx, sep='&'):
    ''' '''
    # Get VEP
    try: val_get = vnt_obj.get_tag_value(VEPtag)
    except Exception: return None
    #end try
    trscrpt_clean = []
    trscrpt_list = val_get.split(',')
    # Clean transcripts
    for trscrpt in trscrpt_list:
        trscrpt_split = trscrpt.split('|')
        # Get dbnsfp_ENST
        dbnsfp_ENST = trscrpt_split[dbnsfp_ENST_idx].split(sep)
        if len(dbnsfp_ENST) >= 1: # need to assign values by transcripts
            dbnsfp_idx = -1
            trscrpt_ENST = trscrpt_split[ENST_idx]
            # Check index for current transcript in dbNSFP if any
            for i, ENST in enumerate(dbnsfp_ENST):
                if ENST == trscrpt_ENST :
                   dbnsfp_idx = i
                   break
                #end if
            #end for
            for k, v in dbNSFP_fields.items():
                if dbnsfp_idx >= 0:
                    val_ = trscrpt_split[v].split(sep)[dbnsfp_idx]
                    if val_ == '.': val_ = ''
                    #end if
                    trscrpt_split[v] = val_
                else:
                    trscrpt_split[v] = ''
                #end if
            #end for
            trscrpt_clean.append('|'.join(trscrpt_split))
        else:
            trscrpt_clean.append(trscrpt)
        #end if
    #end for
    return ','.join(trscrpt_clean)
#end def

def get_variant_result_file_header():
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
    return header

def get_variant_result_higlass_file_header():
    header = '##fileformat=VCFv4.3\n'
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    return header

def parse_regenie_results(regenie_output):
    # Extract Regenie results
    regenie_results = {}
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
    with gzip.open(regenie_output, 'r') as f_in:
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

            regenie_results[id] = {}
            regenie_results[id]["regenie_ml10p"] = r_log10p
            regenie_results[id]["regenie_beta"] = r_beta
            regenie_results[id]["regenie_chisq"] = r_chisq
            regenie_results[id]["regenie_se"] = r_se
            regenie_results[id]["regenie_test"] = r_test
    return regenie_results