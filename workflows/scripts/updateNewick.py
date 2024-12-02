from ete3 import Tree
import csv

def read_csv_mapping(csv_file):
    """
    Reads a CSV file mapping sequence numbers to actual names.
    Expects two columns: 'number' and 'name'.
    :param csv_file: Path to the CSV file.
    :return: Dictionary mapping sequence numbers to names.
    """
    mapping = {}
    with open(csv_file, mode='r') as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) != 2:
                continue
            number, name = row
            mapping[number] = name
    return mapping

def replace_names_in_newick(newick_file, mapping, output_file):
    """
    Replaces sequence numbers with actual names in a Newick file.
    :param newick_file: Path to the Newick file.
    :param mapping: Dictionary mapping numbers to names.
    :param output_file: Path to save the modified Newick file.
    """
    # Load the Newick tree
    tree = Tree(newick_file, format=1)

    # Replace node names based on the mapping
    for leaf in tree:
        if leaf.name in mapping:
            leaf.name = mapping[leaf.name]

    # Write the modified tree to a new file
    tree.write(outfile=output_file, format=1)
    print(f"Modified Newick tree saved to: {output_file}")

if __name__ == "__main__":
    import argparse

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Replace sequence numbers in a Newick file with actual names.")
    parser.add_argument("newick_file", help="Path to the input Newick file.")
    parser.add_argument("csv_file", help="Path to the CSV file containing number-to-name mapping.")
    parser.add_argument("output_file", help="Path to save the modified Newick file.")

    args = parser.parse_args()

    # Read the mapping and process the Newick file
    mapping = read_csv_mapping(args.csv_file)
    replace_names_in_newick(args.newick_file, mapping, args.output_file)
