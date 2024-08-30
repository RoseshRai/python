from collections import Counter


class CholeraePathogenicBacterium:
    def __init__(self, fasta_file: str):
        self.fasta_file = fasta_file
        self.nucleotide_counts: Counter[str] = Counter()

    def count_nucleotides(self) -> dict[str, int]:
        """Count the number of each nucleotide (A, T, C, G) in the FASTA file."""
        with open(self.fasta_file, "r") as file:
            for line in file:
                if not line.startswith(">"):
                    self.nucleotide_counts.update(line.strip())
        return dict(self.nucleotide_counts)


if __name__ == "__main__":
    print(
        "Nucleotide counts:",
        CholeraePathogenicBacterium("../resources/vibrio_cholerae_pathogenic_bacterium.fasta").count_nucleotides(),
    )
