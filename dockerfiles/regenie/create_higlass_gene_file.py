import click
import json

@click.command()
@click.help_option("--help", "-h")
@click.option("-r", "--regenie-output", required=True, type=str, help="Regenie output file")
@click.option("-g", "--gene-info", required=True, type=str, help="Gene inserts file from portal")
@click.option("-a", "--aaf-bin", required=True, type=str, help="AAF bin to extract")
@click.option("-o", "--out", required=True, type=str, help="the output file name")
def main(regenie_output, gene_info, aaf_bin, out):
    """This script takes a gene-based regenie output file and transforms it in a vcf 
       that the Higlass browser can understand

    Example usage: 

    python create_higlass_gene_file.py -r /path/to/out.regenie  -g gene_inserts_from_portal.json -o higlass_gene_tests.vcf

    """

    gene_mapping = {}

    with open(gene_info, 'r') as f:
        data = json.load(f)
        for genes in data:
            gene_id = genes["ensgid"]
            gene_mapping[gene_id] = {
                "chrom": genes["chrom"],
                "start": genes["spos"],
                "end": genes["epos"],
                "symbol": genes["gene_symbol"]
            }
    #print(gene_mapping)

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

    with open(regenie_output, 'r') as f_in:
        
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
            if test_used == "BURDEN":
                if aaf_bin == "1" and (not regenie_id.endswith("all")):
                    continue
                if aaf_bin != "1" and (not regenie_id.endswith(aaf_bin)):
                    continue
            else: # all other tests only consider all variants in the mask
                if not regenie_id.endswith("all"):
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

            info = f"END={gene_end};NAME={gene_symbol};"
            for mask in regenie_results[gene_id].keys():
                for test in regenie_results[gene_id][mask].keys():
                    mask_u = mask.upper()
                    p_value = regenie_results[gene_id][mask][test]
                    info += f"{mask_u}_{test}={p_value};"

            f_out.write(f"{chr}\t{gene_start}\t{gene_id}\t.\t.\t0\tPASS\t{info}\n")

        

if __name__ == "__main__":
    main()