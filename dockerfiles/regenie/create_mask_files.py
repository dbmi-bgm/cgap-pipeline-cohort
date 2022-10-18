import click
from granite.lib import vcf_parser

@click.command()
@click.help_option("--help", "-h")
@click.option("-a", "--annotated-vcf", required=True, type=str, help="VEP annotated VCF, filteres and with IDs")
def main(annotated_vcf):
    """This script takes an annotated VCF file as input and created the annotations and mask files needed by regenie

    Example usage: 

    python create_phenotype.py -a /path/to/annotated_vcf.vcf

    """

    

    vcf_obj = vcf_parser.Vcf(annotated_vcf)
    idx_gene = vcf_obj.header.get_tag_field_idx('CSQ', 'Gene')
    idx_consequence = vcf_obj.header.get_tag_field_idx('CSQ', 'Consequence')
    #idx_cadd_raw = vcf_obj.header.get_tag_field_idx('CSQ', 'CADD_RAW')
    idx_cadd_phred = vcf_obj.header.get_tag_field_idx('CSQ', 'CADD_PHRED')

    HIGH_CADD_TRHESHOLD = 12.0

    """
    In order to create the set list file, we need to aggregate the data in the right way
    set_list_data = {
        "<GENE>": {
            "chr": <CHR>,
            "pos": <POS>, # According to docs, the "physical positin of the gene". We are using the position of the first encountered variant
            "variants": [<VAR_ID_1>,<VAR_ID_2>,...]
        }
    }
    """
    set_list_data = {}


    """
    Variants are grouped into the following categories - currently hardcoded:
    missense_no_cadd
    no_missense_no_cadd
    missense_low_cadd
    no_missense_low_cadd
    missense_high_cadd
    no_missense_high_cadd
    """
    with open("regenie_input.annotation", "w") as output_file:
        for record in vcf_obj.parse_variants():
            info_split = record.INFO.split(";")
            for field in info_split:
                if "CSQ=" in field:
                    csq_split = field.split("|")
                    id = record.ID
                    gene_symbol = csq_split[idx_gene]
                    consequence = csq_split[idx_consequence]
                    is_missense = consequence == "missense_variant"
                    cadd_phred = float(csq_split[idx_cadd_phred]) if csq_split[idx_cadd_phred] else False
                    category = "no_missense_low_cadd"
                    if is_missense and cadd_phred == False:
                        category = "missense_no_cadd"
                    if not is_missense and cadd_phred == False:
                        category = "no_missense_no_cadd"
                    elif is_missense and cadd_phred >= HIGH_CADD_TRHESHOLD:
                        category = "missense_high_cadd"
                    elif not is_missense and cadd_phred >= HIGH_CADD_TRHESHOLD:
                        category = "no_missense_high_cadd"
                    elif is_missense and cadd_phred < HIGH_CADD_TRHESHOLD:
                        category = "missense_low_cadd"
                    #print(record.ID, gene_symbol, category)
                    output_file.write(f'{id} {gene_symbol} {category}\n')

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


    # Create the mask file - currently hardcoded
    with open("regenie_input.masks", "w") as output_file:
        output_file.write('mask_missense missense_no_cadd,missense_low_cadd,missense_high_cadd\n')
        output_file.write('mask_cadd missense_high_cadd,no_missense_high_cadd\n')
        output_file.write('mask_missense_plus_cadd missense_no_cadd,missense_low_cadd,missense_high_cadd,no_missense_high_cadd\n')


    

if __name__ == "__main__":
    main()