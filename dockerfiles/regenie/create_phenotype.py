import click

@click.command()
@click.help_option("--help", "-h")
@click.option("-s", "--sample-file", required=True, type=str, help="BGEN sample file")
@click.option("-c", "--cases", required=True, type=str, help="Comma separaten list of affected sample IDs")
@click.option("-o", "--output", required=True, type=str, help="File name of phenotype file")
def main(sample_file, cases, output):
    """This script takes a BGEN sample file and produces a phenotype file for use in regenie. 

    Example usage: 

    python create_phenotype.py -s /path/to/sample/MSA.sample -c ID1,ID2,ID3,ID4 -o MSA.phenotype

    python create_phenotype.py -s /path/to/sample/MSA.sample -o MSA.phenotype -c MSA_singleton_99_40-WES,MSA_singleton_782532_01_C_1_S8-WES,MSA_singleton_551913_P1_1_EXO_551913_01-WES,MSA_singleton_541376_01_C_1_S7-WES,MSA_singleton_300683_01_C_1_S6-WES,MSA_singleton_22699_MSA_A_1_S23-WES,MSA_singleton_22692_MSA_A_2_S24-WES,MSA_singleton_22655_MSA_A_2_S21-WES,MSA_singleton_22571_MSA_A_1_S20-WES,MSA_singleton_22561_MSA_A_1_S19-WES,MSA_singleton_22549_MSA_A_2_S18-WES,MSA_singleton_22527_MSA_A_1_S17-WES,MSA_singleton_22514_MSA_A_1_S16-WES,MSA_singleton_22494_MSA_A_1_S15-WES,MSA_singleton_22365_MSA_A_1_S14-WES,MSA_singleton_198677_01_C_1_S5-WES,MSA_singleton_171099_01_MSA_EXO-WES,MSA_singleton_14-49_A_2_S13-WES,MSA_singleton_126939_01_C_1_S2-WES,MSA_singleton_12_18-WES,MSA_singleton_11_46-WES,MSA_singleton_100584_01_C_2_S1-WES,MSA_singleton_07_36-WES,MSA_singleton_07_03-WES,MSA_singleton_04_56-WES,MSA_singleton_04_51-WES,MSA_singleton_03_55-WES,MSA_994562_P1_1_EXO_994562_01-WES,MSA_994562_M1_2_EXO_994562_02-WES,MSA_994562_F1_1_EXO_994562_03-WES,MSA_985648_S1_EXO_S1_MSA-WES,MSA_985648_P1_EXO_Proband_MSA-WES,MSA_985648_M_EXO_M_MSA-WES,MSA_828071_S3_2_EXO_828071_06-WES,MSA_828071_S1_1_EXO_828071_04-WES,MSA_828071_P1_2_EXO_828071_01-WES
    """

    case_list = cases.split(",")

    with open(sample_file) as input_file, open(output, "w") as output_file:
        output_file.write("FID IID Y1\n")
        
        next(input_file) # Skip header line
        next(input_file) # Skip 0 0 line
        for line in input_file:
            line_split = line.split()
            sample_id = line_split[1]
            case_control = "1" if sample_id in case_list else "0"
            output_file.write(line_split[0]+" "+sample_id+" "+case_control+"\n")


if __name__ == "__main__":
    main()