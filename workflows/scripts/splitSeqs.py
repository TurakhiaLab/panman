from Bio import SeqIO
import os


def splitFasta(input_file, output_dir):

    seqNameMap=""
    ff = open("temp_dir", "w")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    count=1


    for record in SeqIO.parse(input_file, "fasta"):
        filename = f"{count}.fasta"
        seqNameMap = str(count) + "," + record.id
        ff.write(seqNameMap)
        ff.write("\n")
        output_path = os.path.join(output_dir, filename)

        with open(output_path, "w") as output_file:
            SeqIO.write(record, output_file, "fasta")
        count+=1
        print(f"Saved: {output_path}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Split sequences in a FASTA file into individual files.")
    parser.add_argument("inp", help="Path to the input FASTA file.")
    parser.add_argument("outDir", help="Directory to save the split FASTA files.")

    args = parser.parse_args()

    splitFasta(args.inp, args.outDir)

