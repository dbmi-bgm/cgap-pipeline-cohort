import click, json, os
from utils import VALID_GENOTYPES
from granite.lib import vcf_parser

CHUNK_SIZE = 1000000
CHUNK_PREFIX = "variant_details_chunk"

@click.command()
@click.help_option("--help", "-h")
@click.option("-a", "--annotated-vcf", required=True, type=str, help="Jointly called, annotated and filtered VCF")
@click.option("-s", "--sample-info", required=True, type=str, help="Encoded JSON with sample information")
@click.option("-o", "--output", required=True, type=str, help="File name of details file")
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

    num_variants = 0
    chunk = 0
    chunk_files = []
    f_out = None

    for record in vcf_obj.parse_variants():

        if num_variants % CHUNK_SIZE == 0:
            compress_and_close_chunk(chunk-1, f_out)
            chunk_file = f"{CHUNK_PREFIX}_{chunk}.vcf"
            f_out = open(chunk_file, "w")
            chunk_files.append(f"{chunk_file}.gz")
            vcf_obj.write_header(f_out)
            chunk += 1
        
        num_variants += 1

        samples = record.IDs_genotypes
        GT_idx = record.FORMAT.split(":").index("GT")
        ##samples is comma separated list with SAMPLE_ID:LINKTO_ID:IS_AFFECTED:TISSUE_TYPE:CONTACT
        info="samples="
        for sample in samples:
            gt = record.GENOTYPES[sample].split(":")[GT_idx]
            if gt in VALID_GENOTYPES and gt not in ["./.", "0/0", "0|0"]:
                linkto_id = sample_info_dict[sample]["linkto_id"]
                is_affected = sample_info_dict[sample]["is_affected"]
                tissue_type = sample_info_dict[sample]["tissue_type"]
                contact = sample_info_dict[sample]["contact"]
                info += f"{sample}:{linkto_id}:{is_affected}:{tissue_type}:{contact},"


        f_out.write(f"{record.CHROM}\t{record.POS}\t{record.ID}\t{record.REF}\t{record.ALT}\t0\tPASS\t{info}\n")
    
    compress_and_close_chunk(chunk-1, f_out)
    # Files need to be in the correct order to produce a sorted vcf. This is required for tabix to work.
    cmd_result = os.system(f"bcftools concat {' '.join(chunk_files)} -o {output} --threads 6 -O z  || exit 1") 
    if cmd_result > 0: # exit status shows error
        raise Exception(f"bcftools concat command failed.")
    cmd_result = os.system(f"tabix -p vcf -f {output} || exit 1")
    if cmd_result > 0: # exit status shows error
        raise Exception(f"tabix command failed.")
    os.system(f"rm -f {CHUNK_PREFIX}*")

def compress_and_close_chunk(chunk:int, file_handle):
    if chunk < 0 or not file_handle or file_handle.closed:
        return
    file_handle.close()
    file_name = f"{CHUNK_PREFIX}_{chunk}.vcf"
    os.system(f"bgzip --threads 6 -c {file_name} > {file_name}.gz || exit 1")
    os.system(f"rm -f {file_name}")
    
if __name__ == "__main__":
    main()