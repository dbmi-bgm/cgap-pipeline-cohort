import click
import json

@click.command()
@click.help_option("--help", "-h")
@click.option("-s", "--sample-info", required=True, type=str, help="Encoded JSON with sample information")
@click.option("-o", "--output", required=True, type=str, help="Popmap output")
def main(sample_info, output):
    """
    This script takes the same information and creates the popmap file that the HWE filter required
    """

    sample_info_dec = json.loads(sample_info)
    sample_info_dict = {}
    for sample in sample_info_dec:
        sample_id = sample["sample_id"]
        sample_info_dict[sample_id] = sample["ancestry"]


    # Create the set list file from  set_list_data
    with open(output, "w") as output_file:
        for sample_id in sample_info_dict.keys():
            ancestry = sample_info_dict[sample_id]
            output_file.write(f'{sample_id}\t{ancestry}\n')

if __name__ == "__main__":
    main()
