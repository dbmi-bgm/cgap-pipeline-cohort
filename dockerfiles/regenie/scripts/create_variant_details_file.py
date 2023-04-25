import click, json
from utils import VALID_GENOTYPES
from granite.lib import vcf_parser

@click.command()
@click.help_option("--help", "-h")
@click.option("-a", "--annotated-vcf", required=True, type=str, help="Jointly called, annotated and filtered VCF")
@click.option("-s", "--sample-info", required=True, type=str, help="Encoded JSON with sample information")
@click.option("-o", "--output", required=True, type=str, help="File name of phenotype file")
def main(annotated_vcf, sample_info, output):
    """
    This script takes the annotated, filtered VCF and sample
    information and produces a VCF files that contains the variants together with the sample info.

    """

    vcf_obj = vcf_parser.Vcf(annotated_vcf)

    sample_info_dec = json.loads(sample_info)
    sample_info_dict = {}
    for sample in sample_info_dec:
        sample_id = sample["sample_id"]
        sample_info_dict[sample_id] = {
            "linkto_id": sample["linkto_id"],
            "is_affected": bool(sample["is_affected"]),
            "tissue_type": sample["tissue_type"] or "",
            "contact": sample["contact"] or "",
        }

    variant_samples = {}
    with open(output, 'w') as f_out:
        f_out.write('##fileformat=VCFv4.3\n')
        f_out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for record in vcf_obj.parse_variants():
            id = record.ID
            samples = record.IDs_genotypes
            GT_idx = record.FORMAT.split(":").index("GT")
            variant_samples[id] = []
            ##samples is comma separated list with SAMPLE_ID:LINKTO_ID:IS_AFFECTED:TISSUE_TYPE:CONTACT
            info="samples="
            for sample in samples:
                gt = record.GENOTYPES[sample].split(":")[GT_idx]
                if gt in VALID_GENOTYPES and gt not in ["./.", "0/0", "0|0"]:
                    variant_samples[id].append(sample)
                    linkto_id = sample_info_dict[sample]["linkto_id"]
                    is_affected = sample_info_dict[sample]["is_affected"]
                    tissue_type = sample_info_dict[sample]["tissue_type"]
                    contact = sample_info_dict[sample]["contact"]
                    info += f"{sample}:{linkto_id}:{is_affected}:{tissue_type}:{contact},"

            f_out.write(f"{record.CHROM}\t{record.POS}\t{record.ID}\t{record.REF}\t{record.ALT}\t0\tPASS\t{info}\n")

    
if __name__ == "__main__":
    main()