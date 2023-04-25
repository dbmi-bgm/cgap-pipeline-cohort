import click
from granite.lib import vcf_parser
from utils import get_worst_consequence, get_worst_transcript

@click.command()
@click.help_option("--help", "-h")
@click.option("-a", "--annotated-vcf", required=True, type=str, help="VEP annotated VCF, filteres and with IDs")
@click.option("-c", "--high-cadd-threshold", required=True, type=float, help="High CADD threshold")
def main(annotated_vcf, high_cadd_threshold):
    """This script takes an annotated VCF file as input and created the annotations and mask files needed by regenie

    Example usage: 

    python create_phenotype.py -a /path/to/annotated_vcf.vcf

    """

    VEP_TAG = 'CSQ'

    vcf_obj = vcf_parser.Vcf(annotated_vcf)
    idx_gene = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'Gene')
    idx_consequence = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'Consequence')
    idx_canonical = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'CANONICAL')
    idx_cadd_phred = vcf_obj.header.get_tag_field_idx(VEP_TAG, 'CADD_PHRED')


    """
    In order to create the set list file, we need to aggregate the data in the right way
    set_list_data = {
        "<GENE>": {
            "chr": <CHR>,
            "pos": <POS>, # According to docs, the "physical position of the gene". We are using the position of the first encountered variant
            "variants": [<VAR_ID_1>,<VAR_ID_2>,...]
        }
    }
    """
    set_list_data = {}


    """
    Variants are assigned the following categories:
    missense
    high_cadd
    rare
    They are combined into approriate masks below
    """
    with open("regenie_input.annotation", "w") as output_file:
        all_categories = []
        for record in vcf_obj.parse_variants():
            id = record.ID
            vep_tag_value = record.get_tag_value(VEP_TAG)
            worst_transcript = get_worst_transcript(vep_tag_value, idx_canonical, idx_consequence)
            worst_transcript_ = worst_transcript.split('|')
            worst_consequence = get_worst_consequence(worst_transcript_[idx_consequence])
            gene_symbol = worst_transcript_[idx_gene]
            is_missense = worst_consequence == "missense_variant"
            is_nonsense = worst_consequence == "stop_gained"
            is_essential_splice = (worst_consequence == "splice_acceptor_variant") or (worst_consequence == "splice_donor_variant")
            cadd_phred = float(worst_transcript_[idx_cadd_phred]) if worst_transcript_[idx_cadd_phred] else False
            is_high_cadd = cadd_phred >= high_cadd_threshold
            af = float(record.get_tag_value("AF"))


            categories = []
            if is_missense:
                categories.append("missense")
            if is_high_cadd:
                categories.append("high_cadd")
            if is_nonsense:
                categories.append("nonsense")
            if is_essential_splice:
                categories.append("essential_splice")


            if len(categories) > 0:
                category = "_".join(categories)
                all_categories.append(category)
                output_file.write(f'{id} {gene_symbol} {category}\n')
            else: 
                # The None category is not present in the mask file and is
                # therefore ignored downstream
                output_file.write(f'{id} {gene_symbol} None\n')

            if gene_symbol not in set_list_data:
                set_list_data[gene_symbol] = {
                    "chr": record.CHROM,
                    "pos": str(record.POS),
                    "variants": [id]
                }
            else:
                set_list_data[gene_symbol]["variants"].append(id)


    # Create the set list file from  set_list_data
    with open("regenie_input.set_list", "w") as output_file:
        for gene in set_list_data.keys():
            data = set_list_data[gene]
            variants_str = ','.join(data["variants"])
            output_file.write(f'{gene} {data["chr"]} {data["pos"]} {variants_str}\n')

    # remove duplicate values
    all_categories = list(set(all_categories))
    all_categories.sort()

    # Create the mask file - currently hardcoded
    with open("regenie_input.masks", "w") as output_file:
        missense_categories = [cat for cat in all_categories if "missense" in cat]
        output_file.write(f'mask_missense {",".join(missense_categories)}\n')
        high_cadd_categories = [cat for cat in all_categories if "high_cadd" in cat]
        output_file.write(f'mask_cadd {",".join(high_cadd_categories)}\n')
        missense_cadd_categories = [cat for cat in all_categories if "missense" in cat and "high_cadd" in cat]
        output_file.write(f'mask_missense_cadd {",".join(missense_cadd_categories)}\n')
        nonsense_splice_categories = [cat for cat in all_categories if "nonsense" in cat or "essential_splice" in cat]
        output_file.write(f'mask_nonsense_splice {",".join(nonsense_splice_categories)}\n')


if __name__ == "__main__":
    main()