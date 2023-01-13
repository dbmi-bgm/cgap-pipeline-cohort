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