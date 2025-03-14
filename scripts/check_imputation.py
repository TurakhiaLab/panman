import argparse  # For parsing command-line arguments

# IUPAC ambiguity codes and their possible resolutions
IUPAC_CODES = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G", "R"}, "Y": {"C", "T", "Y"}, "S": {"G", "C", "S"},
    "W": {"A", "T", "W"}, "K": {"G", "T", "K"}, "M": {"A", "C", "M"},
    "B": {"C", "G", "T", "B", "Y", "S", "K"},
    "D": {"A", "G", "T", "D", "R", "W", "K"},
    "H": {"A", "C", "T", "H", "Y", "W", "M"},
    "V": {"A", "C", "G", "V", "R", "S", "M"},
    "N": {"A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"}
}

def read_fasta(filepath: str) -> dict[str, str]:
    """Read a FASTA file and return a dictionary of sequences indexed by sequence ID."""
    sequences = {}
    with open(filepath, "r") as f:
        seq_id = None
        seq_list = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id is not None:
                    sequences[seq_id] = "".join(seq_list)
                seq_id = line[1:]  # Remove '>' from header
                seq_list = []
            else:
                seq_list.append(line)
        if seq_id is not None:
            sequences[seq_id] = "".join(seq_list)
    return sequences


def compare_sequences(original: dict[str, str], imputed_filepath: str):
    """Compare original and imputed sequences line by line to reduce memory usage."""
    illegal_change_ids = []
    current_seq_id = None
    imputed_seq = []
    
    with open(imputed_filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq_id is not None:
                    orig_seq = original.get(current_seq_id, None)
                    imputed_seq = "".join(imputed_seq)
                    if orig_seq and len(orig_seq) == len(imputed_seq):
                        is_illegal = [i not in IUPAC_CODES.get(o, {o}) for o, i in zip(orig_seq, imputed_seq)]
                        if any(is_illegal):
                            illegal_change_ids.append(current_seq_id)
                            print(f"Error: {current_seq_id}")
                            print("".join(orig_seq[i] for i in range(len(orig_seq)) if is_illegal[i]))
                            print("".join(imputed_seq[i] for i in range(len(orig_seq)) if is_illegal[i]))
                            return
                    elif orig_seq:
                        print(f"Error: Sequence length mismatch for {current_seq_id}")
                        return
                current_seq_id = line[1:]  # Remove '>' from header
                imputed_seq = []
            else:
                imputed_seq.append(line)
        
        if current_seq_id is not None:
            orig_seq = original.get(current_seq_id, None)
            imputed_seq = "".join(imputed_seq)
            if orig_seq and len(orig_seq) == len(imputed_seq):
                if any(i not in IUPAC_CODES.get(o, {o}) for o, i in zip(orig_seq, imputed_seq)):
                    illegal_change_ids.append(current_seq_id)
            elif orig_seq:
                print(f"Error: Sequence length mismatch for {current_seq_id}")
                return
    
    if illegal_change_ids:
        print("Sequences with illegal changes:", ", ".join(illegal_change_ids))
    else:
        print("No illegal changes detected.")


def main():
    parser = argparse.ArgumentParser(description="Check imputation results between original and imputed FASTA files.")
    parser.add_argument("--original", required=True, help="Path to the original FASTA file.")
    parser.add_argument("--imputed", required=True, help="Path to the imputed FASTA file.")
    args = parser.parse_args()

    original_sequences = read_fasta(args.original)
    compare_sequences(original_sequences, args.imputed)


if __name__ == "__main__":
    main()
