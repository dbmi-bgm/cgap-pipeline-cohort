import click
import json
import os


@click.command()
@click.help_option("--help", "-h")
@click.option("-a", "--annotated-vcf", required=True, type=str, help="Annotated VCF file")
@click.option("-s", "--sample-info", required=True, type=str, help="Encoded JSON with sample information")
def main(annotated_vcf, sample_info):
    """
    Runs Peddy and returns an updated sample info string (JSON encoded)
    """

    PEDDY_FAM_FILE = "tmp.peddy.fam"
    PEDDY_OUT_PREFIX = "tmp.peddy_out"

    sample_info_dec = json.loads(sample_info)

    with open(PEDDY_FAM_FILE, "w") as fam:
        for sample in sample_info_dec:
            sample_id = sample["sample_id"]
            is_case = sample["is_affected"]
            phenotype_value = 2 if is_case else 1
            fam.write(f'{sample_id}\t{sample_id}\t0\t0\t0\t{phenotype_value}\n')


    os.system(f'peddy -p 4 --sites hg38 --prefix {PEDDY_OUT_PREFIX} {annotated_vcf} {PEDDY_FAM_FILE} || exit 1')

    # parse Peddy output
    # The Peddy output file has the following header
    # #family_id	sample_id	paternal_id	maternal_id	sex	phenotype	duplicates	het_call_rate	het_ratio	het_mean_depth	het_idr_baf	ancestry-prediction	PC1	PC2	PC3	sex_het_ratio	sex_fixed
    sample_id_idx = 1
    ancestry_pred_idx = 11

    ancestry_dict = {}
    with open(PEDDY_OUT_PREFIX + ".peddy.ped", "r") as p_out:
        for line in p_out:
            if line.startswith("#"):
                continue
            line_ = line.split("\t")
            sample_id = line_[sample_id_idx]
            ancestry_pred = line_[ancestry_pred_idx]
            ancestry_dict[sample_id] = ancestry_pred

    # Add to ancestry to sample_info
    for sample in sample_info_dec:
        sample_id = sample["sample_id"]
        sample["ancestry"] = ancestry_dict[sample_id]

    # This is the "return" value that will be assigned to sample_info
    print(json.dumps(sample_info_dec))

if __name__ == "__main__":
    main()
