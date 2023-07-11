import click
import csv, gzip

@click.command()
@click.help_option("--help", "-h")
@click.option("-r", "--regenie-output", required=True, type=str, help="Regenie output file")
@click.option("-s", "--snp-list", required=True, type=str, help="Mask SNP list file")
@click.option("-g", "--gene-info", required=True, type=str, help="Gene inserts file from portal")
@click.option("-a", "--aaf-bin", required=True, type=str, help="AAF bin to extract")
@click.option("-o", "--out", required=True, type=str, help="the output file name")
def main(regenie_output, snp_list, gene_info, aaf_bin, out):
    """This script takes a gene-based regenie output file and transforms it in a vcf 
       that the Higlass browser can understand

    Example usage: 

    python create_higlass_gene_file.py -r /path/to/out.regenie  -g gene_inserts_from_portal.tsv -o higlass_gene_tests.vcf

    """

    gene_mapping = {}

    with open(gene_info, 'r') as f:
        tsv_file = csv.reader(f, delimiter="\t")
     
        # header: ens_id	chr	start	end	strand	gene
        for info in tsv_file:
            if info[0] == "ens_id":
                continue

            gene_id = info[0]
            gene_mapping[gene_id] = {
                "chrom": info[1],
                "start": info[2],
                "end": info[3],
                "symbol": info[5]
            }

    #print(gene_mapping)

    # contains information about which mask contains which SNPs
    mask_snp_list = {}
    with gzip.open(snp_list, 'r') as f:
        for line in f:

            # Example line
            # ENSG00000164002.mask_cadd.0.01	chr1_40515148_C_T,chr1_40515335_CTG_C

            line_ = line.strip().split("\t")
            mask_id = line_[0]
            snps = line_[1]

            if aaf_bin != "1" and (not mask_id.endswith(aaf_bin)):
                continue

            mask_id_ = mask_id.split(".")
            gene_id = mask_id_[0]
            mask = mask_id_[1].upper()

            if gene_id not in mask_snp_list:
                mask_snp_list[gene_id] = ""
            mask_snp_list[gene_id] += f"{mask}_SNPS={snps};"



    r_map = {
        "CHROM": 0,
        "GENPOS": 1,
        "ID": 2,
        "ALLELE0": 3,
        "ALLELE1": 4,
        "A1FREQ": 5,
        "N": 6,
        "TEST": 7,
        "BETA": 8,
        "SE": 9,
        "CHISQ": 10,
        "LOG10P": 11,
        "EXTRA": 12,

    }

    genes_without_annotations = []

    # Formatted regnie results:
    # regenie_results = {
    #     "ENSG00000177465": {
    #         "mask_cadd": {
    #             "BURDEN": 0.9,
    #             "SKAT": 0.8,
    #             ...
    #         },
    #         ...
    #     },
    #     ...
    # }
    regenie_results = {}

    with gzip.open(regenie_output, 'r') as f_in:
        
        for line in f_in:
            if line.startswith("##") or line.startswith("CHROM"):
                continue

            # Example line
            # CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
            # 1 169859186 ENSG00000000457.mask_cadd.all ref mask_cadd.all NA 67 ADD-ACATO NA NA 2.73635 1.00838 NA

            result = line.strip().split()

            regenie_id = result[r_map["ID"]]
            regenie_id_arr = regenie_id.split(".")
            gene_id = regenie_id_arr[0]
            mask = regenie_id_arr[1]

            if gene_id in genes_without_annotations:
                continue

            if gene_id not in gene_mapping:
                print(f"WARNING: {gene_id} not found in annotations. Skipping.")
                genes_without_annotations.append(gene_id)
                continue
            
            p_value = result[r_map["LOG10P"]]
            regenie_test = result[r_map["TEST"]]
            test_used = ""
            if regenie_test.startswith("ADD-"):
                test_used = regenie_test.split("-")[1]
            elif regenie_test == "ADD":
                test_used = "BURDEN"
            else:
                print(f"WARNING: {regenie_test} not recognized")
                continue

            # Select the right lines based on the test used.
            # if test_used == "BURDEN":
            #     if aaf_bin == "1" and (not regenie_id.endswith("all")):
            #         continue
            #     if aaf_bin != "1" and (not regenie_id.endswith(aaf_bin)):
            #         continue
            # else: # all other tests only consider all variants in the mask
            #     if not regenie_id.endswith("all"):
            #         continue
            
            # Select the right lines based on the test used.
            if aaf_bin == "1" and (not regenie_id.endswith("all")):
                continue
            if aaf_bin != "1" and (not regenie_id.endswith(aaf_bin)):
                continue

            if gene_id not in regenie_results:
                regenie_results[gene_id] = {}

            if mask not in regenie_results[gene_id]:
                regenie_results[gene_id][mask] = {}

            regenie_results[gene_id][mask][test_used] = p_value
               
    with open(out, 'w') as f_out:
        f_out.write('##fileformat=VCFv4.3\n')
        f_out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

        for gene_id in regenie_results.keys():
            gene_annotation = gene_mapping[gene_id]

            chr = "chr" + gene_annotation["chrom"]
            gene_start = gene_annotation["start"]
            gene_end = gene_annotation["end"]
            gene_symbol = gene_annotation["symbol"]

            info = f"END={gene_end};SYMBOL={gene_symbol};"
            for mask in regenie_results[gene_id].keys():
                for test in regenie_results[gene_id][mask].keys():
                    mask_u = mask.upper()
                    p_value = regenie_results[gene_id][mask][test]
                    info += f"{mask_u}_{test}={p_value};"
            
            if gene_id in mask_snp_list:
                info += mask_snp_list[gene_id]

            f_out.write(f"{chr}\t{gene_start}\t{gene_id}\t.\t.\t0\tPASS\t{info}\n")

        

if __name__ == "__main__":
    main()