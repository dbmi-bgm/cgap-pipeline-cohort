################################################
#   Libraries
################################################

import click
from granite.lib import vcf_parser
import os


################################################
#   Top level variables
################################################

# Thresholds
FS_SNP = 60.0
FS_INDEL = 200.0
INBREEDING = -0.8
MQ_RANK_SUM = -12.5
QD = 2.0
READ_POS_RANK_SUM_SNP = -8.0
READ_POS_RANK_SUM_INDEL = -20.0
SOR_SNP = 3.0
SOR_INDEL = 10.0

CHUNK_SIZE = 1000000
CHUNK_PREFIX = "gatk_filter_chunk"


################################################
#   Functions
################################################


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-a",
    "--annotated-vcf",
    required=True,
    type=str,
    help="Annotated, jointly called VCF (gzipped)",
)
@click.option(
    "-o",
    "--out",
    required=True,
    type=str,
    help="the output file name of the gzipped VCF after filtering",
)
def main(annotated_vcf, out):
    """This script applies GATK best practice filter. It will be applied to CHUNK_SIZE variants at a time.
    The intermediate filtered VCF files are gzipped and merged in the end. This prevents the creation of an
    uncompressed VCF containing all variants.

    GATK best practices filters
    - Include SNPs: QD > 2.0, FS < 60, MQRankSum > -12.5, ReadPosRankSum > -8.0, SOR <= 3
    - Include Indels: QD > 2.0, FS < 200, ReadPosRankSum > -20.0, SOR < 10.0
    - inbreeding coefficient > -0.8 (vcftools --het > produces .het files which can be checked — optional)


    Example usage:

    python apply_gatk_filter.py -a annotated_vcf.vcf.gz -o annotated_vcf_filtered.vcf.gz

    """

    vcf_obj = vcf_parser.Vcf(annotated_vcf)

    num_variants = 0
    num_excluded = 0
    num_missing_tags = 0


    print(f"SNP filtering: FS < {FS_SNP}, InbreedingCoeff > {INBREEDING}, MQRankSum > {MQ_RANK_SUM}, QD > {QD}, ReadPosRankSum > {READ_POS_RANK_SUM_SNP}, SOR <= {SOR_SNP}")
    print(f"INDEL filtering: FS < {FS_INDEL}, QD > {QD}, ReadPosRankSum > {READ_POS_RANK_SUM_INDEL}, SOR <= {SOR_INDEL}")

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

        try:
            id = record.ID
            ref = record.REF
            alt = record.ALT

            is_indel = len(ref) > 1 or len(alt) > 1

            v_fs = float(record.get_tag_value("FS"))
            v_inbreeding = float(record.get_tag_value("InbreedingCoeff"))
            v_mqrs = float(record.get_tag_value("MQRankSum"))
            v_qd = float(record.get_tag_value("QD"))
            v_rprs = float(record.get_tag_value("ReadPosRankSum"))
            v_sor = float(record.get_tag_value("SOR"))

            # Replace * with -. VEP does not work for record.ALT=="*"
            record.ALT = alt.replace("*", "-")
            record.ID = id.replace("*", "-")

            if is_indel and (
                v_fs < FS_INDEL
                and v_qd > QD
                and v_rprs > READ_POS_RANK_SUM_INDEL
                and v_sor <= SOR_INDEL
            ):
                vcf_obj.write_variant(f_out, record)
            elif (not is_indel) and (
                v_fs < FS_SNP
                and v_inbreeding > INBREEDING
                and v_mqrs > MQ_RANK_SUM
                and v_qd > QD
                and v_rprs > READ_POS_RANK_SUM_SNP
                and v_sor <= SOR_SNP
            ):
                vcf_obj.write_variant(f_out, record)
            else:
                #print(id, is_indel, v_fs, v_inbreeding, v_mqrs, v_qd, v_rprs, v_sor)
                num_excluded += 1
            
        except ValueError: # This is thrown by Granite if the tag is not there
            num_missing_tags += 1
            num_excluded += 1
            pass
        except Exception:
            raise Exception(
                "\nERROR applying GATK filter for variant {0}\n".format(id)
            )

    compress_and_close_chunk(chunk-1, f_out)
    # Files need to be in the correct order to produce a sorted vcf. This is required for tabix to work.
    cmd_result = os.system(f"bcftools concat {' '.join(chunk_files)} -o {out} --threads 6 -O z  || exit 1") 
    if cmd_result > 0: # exit status shows error
        raise Exception(f"bcftools concat command failed.")
    cmd_result = os.system(f"tabix -p vcf -f {out} || exit 1")
    if cmd_result > 0: # exit status shows error
        raise Exception(f"tabix command failed.")
    os.system(f"rm -f {CHUNK_PREFIX}*")

    # Merge chunks and then remove them


    print("GATK best practice filtering done.")
    print(f"Original number of variants: {num_variants}")
    print(f"Variants excluded: {num_excluded}. {num_missing_tags} of those had missing tags.")
    print(f"New number of variants: {num_variants-num_excluded}")



def compress_and_close_chunk(chunk:int, file_handle):
    if chunk < 0 or not file_handle or file_handle.closed:
        return
    file_handle.close()
    file_name = f"{CHUNK_PREFIX}_{chunk}.vcf"
    os.system(f"bgzip --threads 6 -c {file_name} > {file_name}.gz || exit 1")
    os.system(f"rm -f {file_name}")



if __name__ == "__main__":
    main()
