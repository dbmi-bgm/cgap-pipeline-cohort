import click
from utils import get_cases

@click.command()
@click.help_option("--help", "-h")
@click.option("-s", "--sample-file", required=True, type=str, help="BGEN sample file")
@click.option("-c", "--sample-info", required=True, type=str, help="Encoded JSON with sample information")
@click.option("-o", "--output", required=True, type=str, help="File name of phenotype file")
def main(sample_file, sample_info, output):
    """This script takes a BGEN sample file and produces a phenotype file for use in regenie. 
    """

    case_list = get_cases(sample_info)

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