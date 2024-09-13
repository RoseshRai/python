import itertools
import math
from collections import Counter

from module_2 import Tasks

NINE_MER_PROBABILITY: float = 0.25  # It could be A, C, T, G = 1/4


class Module3:
    """
    What is the expected number of occurrences of a 9-mer in 500 random DNA strings, each of length 1000?
    Assume that the sequences are formed by selecting each nucleotide (A, C, G, T) with the same probability (0.25).

    Each DNA sequence is 1000 bases long, and a 9-mer is a sequence of 9 bases. The 9-mer can start at several positions within the DNA sequence, but it needs to have space for all 9 bases. If you start at position 1, the 9-mer will use positions 1 to 9.
    If you start at position 2, it will use positions 2 to 10, and so on.
    1. 1000−9+1=992
    A 9-mer can start at 992 different positions in a sequence of 1000 letters. This is because we count the first position it starts on meaning at the first 9-mer here would be 1-9 second would be 2-10 (as you can see were count the 2 as the first position)
    2. 500 sequences, each with 992 places for a 9-mer, gives:
        500×992=496000

    3.Each letter has a 1 in 4 chance (0.25) of being A, C, G, or T
    For 9 letters, the probability is
    0.25^9=0.0000038147

    4. Multiply the number of spots (496,000) by the probability of the 9-mer appearing:
    496000 × 0.0000038147 = 1.89209
    """

    @staticmethod
    def calculate_9_mer_occurrences(sequence_length: int, number_of_sequences: int) -> float:
        num_of_possible_9_mers: int = sequence_length - 9 + 1
        total_9_mers_in_all_sequences: int = number_of_sequences * num_of_possible_9_mers
        probability_of_specific_9_mer: float = NINE_MER_PROBABILITY**9

        return total_9_mers_in_all_sequences * probability_of_specific_9_mer

    """
    Implement the below Pseudo code:
        MotifEnumeration(Dna, k, d)
        Patterns ← an empty set
        for each k-mer Pattern in the first string in Dna
            for each k-mer Pattern’ differing from Pattern by at most d mismatches
                if Pattern' appears in each string from Dna with at most d mismatches
                    add Pattern' to Patterns
        remove duplicates from Patterns
        return Patterns 
    """

    @staticmethod
    def motif_enumeration(dna_list: list, motif_len: int, allowed_matches: int, module_2: Tasks) -> set:
        patterns = set()
        first_dna = dna_list[0]

        for i in range(len(first_dna) - motif_len + 1):
            pattern = first_dna[i : i + motif_len]
            neighborhood = module_2.neighbors(pattern, allowed_matches)

            for neighbor in neighborhood:
                if all(
                    module_2.approximate_pattern_count(dna_str, neighbor, allowed_matches) > 0 for dna_str in dna_list
                ):
                    patterns.add(neighbor)

        return patterns

    @staticmethod
    def maximum_possible_motif_score(num_motifs: int, motif_length: int) -> int:
        # In the worst case, at each position, the most frequent nucleotide appears 3 times,
        # and 7 motifs have different nucleotides.
        most_frequent_count = num_motifs // 4 + 1  # 3 times, as explained
        mismatches_per_position = num_motifs - most_frequent_count  # 7 mismatches per position

        # The total score is the sum of mismatches across all positions
        total_score = mismatches_per_position * motif_length

        return total_score

    """
        This function calculates the total entropy of a motif matrix.
    
        :param motif_matrix: List of strings representing the motif matrix
        :return: Total entropy of the motif matrix
    """

    @staticmethod
    def calculate_motif_entropy(motif_matrix: list[str]) -> float:
        motif_matrix = [seq.upper() for seq in motif_matrix]

        num_rows = len(motif_matrix)
        num_cols = len(motif_matrix[0])

        total_entropy = 0.0

        for col in range(num_cols):
            column_nucleotides = [motif_matrix[row][col] for row in range(num_rows)]

            counts = Counter(column_nucleotides)
            total_count = sum(counts.values())

            entropy = 0.0
            for nucleotide in "ACGT":
                frequency = counts.get(nucleotide, 0) / total_count
                if frequency > 0:
                    entropy -= frequency * math.log2(frequency)
            total_entropy += entropy

        return total_entropy

    """
        MedianString(Dna, k)
            distance ← ∞
            for each k-mer Pattern from AA…AA to TT…TT
                if distance > d(Pattern, Dna)
                     distance ← d(Pattern, Dna)
                     Median ← Pattern
            return Median
    """

    @staticmethod
    def generate_kmers(k: int) -> list:
        """Generates all possible k-mers of length k from the nucleotides A, C, G, T."""
        return ["".join(kmer) for kmer in itertools.product("ACGT", repeat=k)]

    @staticmethod
    # Distance between a k-mer Pattern and a set of DNA sequences
    def d(pattern: str, dna: list) -> int:
        """Calculates the sum of minimum Hamming distances between 'pattern' and each DNA sequence in the list."""
        total_distance = 0
        k = len(pattern)

        for sequence in dna:
            min_dist = 999999
            # Slide the pattern along the sequence and calculate the minimum Hamming distance
            for i in range(len(sequence) - k + 1):
                k_mer = sequence[i : i + k]
                min_dist = min(min_dist, Tasks.hamming_distance(pattern, k_mer))
            total_distance += min_dist

        return total_distance

    def median_string(self, dna: list, k: int) -> str:
        """Finds the k-mer (median string) that minimizes the distance between it and a given set of DNA sequences."""
        best_distance = float("inf")
        median = ""

        # Iterate over all possible k-mers
        for pattern in self.generate_kmers(k):
            current_distance = self.d(pattern, dna)
            if current_distance < best_distance:
                best_distance = current_distance
                median = pattern

        return median

    """
        Computes the probability of a k-mer given a profile matrix.
    
        :param k_mer: The k-mer sequence
        :param profile: A dictionary representing the profile matrix with keys 'A', 'C', 'G', 'T'
        :return: The probability of the k-mer given the profile matrix
    """

    @staticmethod
    def compute_probability(k_mer: str, profile: dict) -> float:
        probability = 1.0  # Start with a probability of 1

        for i, nucleotide in enumerate(k_mer):
            probability *= profile[nucleotide][
                i
            ]  # Multiply by the profile probability for the nucleotide at position i

        return probability

    """
        Finds the Profile-most probable k-mer in a given text.
    
        :param text: A DNA string
        :param k: The length of the k-mer
        :param profile: A dictionary representing the profile matrix with keys 'A', 'C', 'G', 'T'
        :return: The Profile-most probable k-mer in text
    """

    def profile_most_probable_kmer(self, text: str, k: int, profile: dict) -> str:
        max_probability = -1.0  # Initialize to a very low value
        most_probable_kmer = text[:k]  # Initialize to the first k-mer as the default

        # Iterate over all k-mers in the text
        for i in range(len(text) - k + 1):
            k_mer = text[i : i + k]
            probability = self.compute_probability(k_mer, profile)  # Compute the probability of the k-mer

            # Update if this k-mer has a higher probability
            if probability > max_probability:
                max_probability = probability
                most_probable_kmer = k_mer

        return most_probable_kmer

    def score_motifs(self, motifs: list[str]) -> int:
        consensus = self.consensus_motif(motifs)
        score = 0
        for motif in motifs:
            score += Tasks.hamming_distance(motif, consensus)
        return score

    # Function to find the consensus sequence for a motif matrix
    @staticmethod
    def consensus_motif(motifs: list[str]) -> str:
        k = len(motifs[0])
        consensus = ""
        for j in range(k):
            counts = {"A": 0, "C": 0, "G": 0, "T": 0}
            for motif in motifs:
                counts[motif[j]] += 1
            consensus += max(counts, key=lambda x: counts[x])
        return consensus

    @staticmethod
    def consensus_from_profile(profile: dict) -> str:
        consensus = ""
        for i in range(len(profile["A"])):
            max_nucleotide = max("ACGT", key=lambda nucleotide: profile[nucleotide][i])
            consensus += max_nucleotide
        return consensus

    # Function to construct a profile matrix from a list of motifs
    @staticmethod
    def build_profile(motifs: list[str]) -> dict:
        k = len(motifs[0])
        profile = {"A": [0.0] * k, "C": [0.0] * k, "G": [0.0] * k, "T": [0.0] * k}

        for j in range(k):
            col_count = {"A": 0, "C": 0, "G": 0, "T": 0}
            for motif in motifs:
                col_count[motif[j]] += 1
            for nucleotide in "ACGT":
                profile[nucleotide][j] = col_count[nucleotide] / len(motifs)

        return profile

    @staticmethod
    def build_profile_with_pseudocounts(motifs: list[str]) -> dict:
        k = len(motifs[0])
        profile = {"A": [0.0] * k, "C": [0.0] * k, "G": [0.0] * k, "T": [0.0] * k}
        t = len(motifs)  # number of motifs

        # Applying Laplace's Rule of Succession: Add 1 to each count, divide by t+4
        for j in range(k):
            col_count = {"A": 1, "C": 1, "G": 1, "T": 1}  # Initialize with pseudocounts (1)
            for motif in motifs:
                col_count[motif[j]] += 1
            for nucleotide in "ACGT":
                profile[nucleotide][j] = col_count[nucleotide] / (t + 4)  # Divide by total + 4

        return profile

    def greedy_motif_search(self, dna: list[str], k: int, t: int) -> list[str]:
        # Initialize BestMotifs as the first k-mers in each string
        best_motifs = [seq[:k] for seq in dna]

        # Iterate over each k-mer in the first string
        for i in range(len(dna[0]) - k + 1):
            motifs = [dna[0][i : i + k]]  # Start with the k-mer from the first string
            # Build motifs iteratively
            for j in range(1, t):
                # profile = self.build_profile(motifs)
                profile = self.build_profile_with_pseudocounts(motifs)
                next_motif = self.profile_most_probable_kmer(dna[j], k, profile)
                motifs.append(next_motif)

            # Update BestMotifs if the current set has a better score
            if self.score_motifs(motifs) < self.score_motifs(best_motifs):
                best_motifs = motifs

        return best_motifs

    """
        Calculates the total distance between a pattern and a list of DNA strings.
    
        :param pattern: The k-mer pattern to compare
        :param dna: A list of DNA strings
        :return: The total distance between the pattern and all DNA strings
    """

    @staticmethod
    def distance_between_pattern_and_strings(pattern: str, dna: list[str]) -> int:
        k = len(pattern)
        total_distance = 0

        for text in dna:
            min_hamming_distance = 99999  # Initialize Hamming distance to infinity

            # Iterate through all k-mers in the current string (text)
            for i in range(len(text) - k + 1):
                k_mer = text[i : i + k]
                current_hamming_distance = Tasks.hamming_distance(pattern, k_mer)

                # Find the minimum Hamming distance between the pattern and k-mer in the string
                if current_hamming_distance < min_hamming_distance:
                    min_hamming_distance = current_hamming_distance

            # Add the minimum Hamming distance for this string to the total distance
            total_distance += min_hamming_distance

        return total_distance

    """
        MedianString(Dna, k)
            distance ← ∞
            Patterns ← AllStrings(k)
            for i ← 0 to |Patterns|
                Pattern ← Patterns[i]
                if distance > DistanceBetweenPatternAndStrings(Pattern, Dna)
                    distance ← DistanceBetweenPatternAndStrings(Pattern, Dna)
                    Median ← Pattern
            return Median
    """


if __name__ == "__main__":
    module_three = Module3()
    # print(module_three.calculate_9_mer_occurrences(1000, 500))
    # print(" ".join(map(str, module_three.motif_enumeration(["CGAAACTTGAAGCCATTCTCTTTGG", "TATGCATAAGGCCGGTCGAATTTCG", "CAAAATAACCAGTATTCCGGTTTGG", "CCCGTCCATACGTTTTTTTGCTACC", "ATAATAACTTTTTGGATACTCAATA", "TACATTGGTATTTAGAGTGTTTCAA"], 5, 1, Tasks()))))

    # motif_matrix_input = [
    #     "TCGGGGgTTTtt",
    #     "cCGGtGAcTTaC",
    #     "AcGGGGATTTtC",
    #     "TtGGGGAcTTtt",
    #     "aaGGGGAcTTCC",
    #     "TtGGGGAcTTCC",
    #     "TCGGGGATTcat",
    #     "TCGGGGATTcCt",
    #     "TaGGGGAacTaC",
    #     "TCGGGtATaaCC"
    # ]
    #
    # print(module_three.calculate_motif_entropy(motif_matrix_input))

    # print(
    #     "".join(
    #         map(
    #             str,
    #             module_three.median_string(
    #                 [
    #                     "ACAGAGGAAACTAGGGCTAATCGATACGACTAAGACCCCACC",
    #                     "GCTTGTGATGCAATAACGGGGCACTCATGCAGGAAGGAAACT",
    #                     "GAAGCTGCCCATACTGCTGTTCAAGCGATACAGTTGATCAGC",
    #                     "TTTGAGGGGGGTCAATATTAAATCAGAACTATTCGCGAAACT",
    #                     "AGCTTGGATGTATAGTGCAATCGTCATAAGCGTATAGAACCT",
    #                     "CCAAGGCGGACCGTCAGGCCGCTCGAAGCTCTGCTGTTCTAA",
    #                     "GAGTGAGCATCAGCCCGGCATCGAATCTATGAAGCTCGTCAC",
    #                     "TTCAGTTCTCAGGGGCATGAACCTCAATCAGTTCTTCTGTGA" "GAAACTGGTAGCCCGACACGCGGAGTCGAGCCTAGATTTTGA",
    #                     "ATGTTTCTAATAGGGGGCTCTCCACAACTCATGGGGAACCT",
    #                 ],
    #                 6,
    #             ),
    #         )
    #     )
    # )
    #
    # profile_input = {
    #     'A': [0.2, 0.2, 0, 0, 0, 0, 0.9, 0.1, 0.1, 0.1, 0.3, 0],
    #     'C': [0.1, 0.6, 0, 0, 0, 0, 0, 0.4, 0.1, 0.2, 0.4, 0.6],
    #     'G': [0, 0, 1, 1, 0.9, 0.9, 0.1, 0, 0, 0, 0, 0],
    #     'T': [0.7, 0.2, 0, 0, 0.1, 0.1, 0, 0.5, 0.8, 0.7, 0.3, 0.4]
    # }
    # sequence_input = "TCGTGGATTTCC"
    #
    # print(module_three.compute_probability(sequence_input, profile_input))

    # profile_input = {
    #     'A': [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
    #     'C': [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
    #     'G': [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
    #     'T': [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]
    # }
    # sequence_input = "GAGCTA"
    #
    # print(module_three.compute_probability(sequence_input, profile_input))

    # profile_input = {
    #     "A": [0.132, 0.276, 0.316, 0.224, 0.211, 0.316, 0.184, 0.25, 0.303, 0.237, 0.184, 0.303, 0.224],
    #     "C": [0.329, 0.237, 0.211, 0.316, 0.316, 0.224, 0.224, 0.211, 0.25, 0.276, 0.276, 0.316, 0.25],
    #     "G": [0.303, 0.184, 0.237, 0.303, 0.263, 0.289, 0.197, 0.25, 0.289, 0.263, 0.25, 0.145, 0.237],
    #     "T": [0.237, 0.303, 0.237, 0.158, 0.211, 0.171, 0.395, 0.289, 0.158, 0.224, 0.289, 0.237, 0.289],
    # }
    # print(
    #     module_three.profile_most_probable_kmer(
    #         "CGAGGCAGGTTGGAAGATTATTTTCCCGTGGACTCCGCCGAAACCAGTGGAGCCCTTGGAAGTCCAATTATGTCGATGAAACCTCTGGGGGCCCGAGCAATTTCTCAGAACTAATGCGAGGAGTGCAAGTGCGATGTAGCCTCACAAATTAATTAGGACGCTCCATGGCCTTAGTATTCGACTACGTGATATATGGAAATCCCCAAATACTAGAGCCCCCTACACACACGCAGTTATTGCCGCGCCTTTGTCAAGGTAACTAGTCGGCATCATATCTCCATCATGCCTTTAGTCGGCCTGGGGCGACATCATGGTGCAGTCATTTTCGGCAAGTGCGAAACCGAGAATATCGTCGTAGTAATGGCTGTGTACTAGCCAGTGAGGGTAGCGGTACACGGAAGTTGTGTCGAGGCGTCACTAGCGTGCCTCTCCTCACACCCAACCTAGTGGTACAATTTGGGAGTGGCCTACATCGTAACTGGTGGTTTGACAATACATCGACGCTCGGTGTTTTTCGGCTGCTAACCGGGCCTGTCGTGTACGGGATTGTGAGTCGAGCCCTAAAATGACTATAGACTGAAGGCGCACTTTACATACCGGGAGAACTGTACGGTATGCAGTCCAGGAAATGAATACAATAAGAGGCGGATTTGTTTGCTAAGCGGTGTTGGTCAAAGCAAGTGATTCGAGGTTGTGCTCCTATGATGTAGTGGTGACCTGAGCACCGCGGAACCGAGTTCTTTGCCTAAAGTGTCCAATAAACGAACCAGAAAATGTTCATGTACGGCATTAAAACCCGGCTCGGCAGCTGCGACCGAATAAGGCGATGCCAAGCTCATAATTCGGTTGTCCCAGTCGGCGCTATTCTGTTCCGTTAGTTGACCAAGGGGTTGTAAAGTAGTCGGTGGCCTGTCATCCACCTCAGTGCGATAGATTCGACGAACCAGCAATCACGCTAGCCACTCTTACAGAAGGCCTGATAAGAATTGTAGCATCGTAA",
    #         13,
    #         profile_input,
    #     )
    # )

    # print(
    #     " ".join(
    #         map(
    #             str,
    #             module_three.greedy_motif_search(
    #                 [
    #                     "GGATACGATATGAACGTATACTCGACCTCTCGATAACGATTAAGTTACACTCATTACGTCCGCTTATACGAATAAGATTTCTCTTAATCCATGAATCGATACATCCCCACATATGACACAGTCTCGTCCGTGCAGCGCTTCATCAAGGTCGCTCTT",
    #                     "GTGCTAAGTGGTAAAGCCATGAATATTTAGCAAGTTAGGAGAGCGATATGGACGATTCTCCACTCCTGACCCCCGGTGCATACTCTAGGCATTGCAGTAACTTGACGGGGGGTGATGTTCGATCGGAGCGGGTAGGAGTCCTTACCCGAGTAAATT",
    #                     "GAAACATGACCATCTGACTTTTCTACTATGGGATATTTCACAAAGAGCGACGGCATGGCCATTTTATGGCTCATAACGTACAATTAGCCCATGAACCCTTTCGCCAGCAGACAGAACTGGCTCGCCGCCGTTGGAAACACGGACCTCCGCACGGAA",
    #                     "ATATTCCACACACCGATGGGTCTCGAGGAATACGCGCTAGATACATCCGGTTCCAAACATACTGTTCACGATGAAAGCGTCGCCTACACCATGAACTTAGGTGACACTCATAGCACCGTTTCATCCATCCACGGACAAACGAGGCGAGTAGGTGGG",
    #                     "GACGATAAGCTGAAACTTCCCGTGGACCTATTGGGAAATCAAGTACGATAGACGTCACGCTAATACTACAGGTATACACAGATTGTGATCAGCAACATGGTGTTTGTCAAGGCCATGAATCTTATTTGTCTCCCCGTCACCAGTGCGACGCATACT",
    #                     "ACATAGAAGAGTATAGTGTCTCAATTGTGCGCTTTTGGCGACCCCAGGCAGCATTCGAATTATCGGGGATAGTATTGAAAGGTGGGAAGCAAAACTAAGTGAACCTATGTCACTGGGTCAGAAACGACGTCTCAGGCCATGAATTACGGTTTGGTC",
    #                     "AAGGGAACCGCCAGGCAGCTGCTTCGCAGCCTCTCTCTATTCAGAGTTCCAACATCCGAACACGTCACATCCTGTTGGGGGAGGCAGACCATGAATGCTCAGGATGGCTCGGAGACTAGTCGTTTTCCCGGGGAACGTCGTTTAGCGCCAGTATTT",
    #                     "GCGCGGCAGAGGAAATCCATGAAATGCCGCCGTATGGGGGTAGGAGACATGAGTTTGTTATAAGCCTCACCCCCGATTGCTCCCGCGGTGGTTGACAGGCGGACCGAAGTCACGTTCAGCCAAGGCCATTGAGCCGCGACTTTGCCGGGGCAGGGC",
    #                     "TGTAGATAAGCTAATGCCATGAACTGTTCCTTTGGTTTCGCCTCACGTTAGCCGAGTTGACCGTGTACTCTGGGCTATAATAAGGTACCCTTGTCCTGGGGTTACTTCAACTTTATCAGCATATTATTCCCTACATACAGTCCTAAATAAAGGGTT",
    #                     "TCAGGTCCGGCAATTGGCAACGTGGATCCAAACATAATCAACTCCAGATAGTCCATGAATGTACCCGCACCAGGCCCGCACGGAAAAAGGGGTTGTAACAGCCTACAACTCTTGGATTACCCCATTTTCCGGAATCCGGGGCTTAGAACCGGTCCG",
    #                     "CCGATGCGAACATAGAGCCTACTCTAGGGCCCATCATAAAGTCACCGGAACCCCTCCTCAACTCAAGTGCATAGCCCCCCTCAGGATGCCATGAACTGCAATTACGTCCTACGCTAGCTGCTACCCACCGCATCGTGGGATATCGGTGGCAGTGAT",
    #                     "GAAAGGTCCAAAACCGAAGAATTATCTCATTCTACTACAATGACATGAAAACCCATGAAGATCCCACCCAGAGGTTGGCGACGACTCACCTCTGTGTTCTGATACTTGAACCGCTCCACTTCGAGACTACTCTTCCAAAACGCCACGATTCATAAC",
    #                     "ACGGTGCCAGTTCAAAGAGATCAAAGTATTGCGTACCGGGTTCCACCGTTACTCGACATTCAGGTCGTAGGTTCCGTCGTACTTGGGACAAGTGTGAATGAGATCCAGAGCGCAAGCTTCTATACCATGAACGAATAGAGTCGGGAACAGGTCGGG",
    #                     "GGGGAGTGGGGCCAGATGATCAACGGAGAAAAGGCGAAACCCATGAAGATCTAGGGGGCGATCGCCTCGTCGTGGGACACTTTTGGCTGCGGCAATCTAGTTGCGTAGAGAAGAGGATCTGGAGTAGTATGAATGGACGAATGTAGGTACCATTAG",
    #                     "CGGGGAATTACCCTTTGCTCGAGTCCCACGCGTTGGACGGTGAGGGTATTGGCTTAAGTCCGTAATCATGCAGAACCCATGAATGAGCGTTGTCACCACGGCCTATGTTTTGACTTCCTCAGGTGTGATGCTGGGAGCAATGTGTACTATTTGAAT",
    #                     "GAGTCCATGAAAGTTTCGCATATCGGCCTTACTGATGTAATAGTGCGACTACGTAGCTACGGCCGACGAAGCATACTCATCGGGCCGCGTTAATACCTGCTCGGCGTTGGCAGGCAGCCGGGCCACATTAACTTCATCCAGAACTGCTAATACAGG",
    #                     "TATTACGCAATCGAGCCTACGGGAATCACTGCCGGAACTTTGGATAGCCCCGCGCGCATAGATCATCAACCCTGTGATAGGCACAAAAGTCGAGCGCCATGGGGTACCCCGCGAAACTAGGAAGCCATGAAACCCCGAGGGTATTAAGGTGTAGGT",
    #                     "CTTCGAAGATGAGTGAACATCTGTCTCCGACTACCTCGAAGGATCCGACGCAGCGTCGCTCGGGCGAAGATGGAGGGGCACCGCTATTCCATGAACGCATAATTGTTCTGGTGTCGAGGACCTTCCCTACCGCAATCACGGCCTTGACTCCATACC",
    #                     "AACTCCATGAATCCACTATTCTCCCTGTCTATCAGCCGATCTTTCTGGTACTGACTTCTAGTTAACCCGCATAAGTAATCCACCCTCTAGTCCAGCGCGCAAGTTGCACCTCAACGGGTGATCCCGACATTGACGCGATCGTGTATTCATATCTAG",
    #                     "TTTCCTTGGGCTCTTGATTAAGCAAAACGCCGGACCCTGATGTCAGTAGAAAACACCTGAAAAACAGATGCCGAAAAAGCCAAAGTTTGCAATTAGGAACGTCATGGAAACTACGTGCGGGATCCCATGAAAAGCAACAGTGTGGAGATTTTATGT",
    #                     "CTTCTTTGGCCACCGAGCTAAATTGTCAGAGTGTACCCGCCCCTTTAATATTCCATGAAGAAACCCTGACCATAAGACCGATGTCTTGTGGCCCTTGGGCGATAGGAACCAATGTCGATCCTCCAACCCTGTTTTTACTGTCTATTTAGCGGTATT",
    #                     "CCGGCACCAGTCAGCGAGTGCTGGTAAGTTGCCCATCACCCGAGAAGTAGTCTTTTGCGATTGCCTGACTAATTCTTGATGCGGGAGGCCATGAATAATATGTTAATGACTTTAACGCTCATTGTATTGTTATCTAGGGTCTGGTTTTCATCGCCC",
    #                     "ATTGAAAGGACTCTTGCTAGCATGACTACACAGAAGCTGATCCTGCACAAGCAGTGTTACGGTCCGGGGCTCCTAGGGGCCGTGCCCGCATCCGCCATTTGAATATATCGACATCAATATTTCGAGAGGCGCCCTCGAAGTCTACATCCCATGAAT",
    #                     "CTCGTAGATGGTGTTATCTAGTTTCTATTGCCCGCCCGATTTCGTGATTGCCAACTATTGACTGTTACTCAGAACTCCATGAACCTCACCTTTGCCAGTGGGTACTCCGGACTCGCGAGCAACTATGTGGCGATGACGTCCGGTAGTCTTAAGAAC",
    #                     "CCGAGAAGGTGTAATTCGCGAGGACGCGCCTTATAAGAGCCCATGTTAAGGTTCTCTCGCTTGTCGGATATGCCTAAACCATTCCGTGATCTTTGCCAATACGAACATTGTAAGATCTGGTCGGTGATGATTCGCCTGGGAACGCAAGCCATGAAG",
    #                 ],
    #                 12,
    #                 25,
    #             ),
    #         )
    #     )
    # )

    # print(
    #     " ".join(
    #         map(
    #             str,
    #             module_three.greedy_motif_search(
    #                 [
    #                     "ATTGGCTCTCTTGCTGTCTAGACGCCGGTAGTAACGGACGCGCTACTTCACGACACAGGCTGTCAGAGGTCTAGAGTGACCAGCCTAGGCCCCTCTGTGGTCGTATACATGGATCAGCTGCAAACACAGGGCGGTCACAGCTTTATCTGCCCATCC",
    #                     "ACTGACTCGACGCGTAATTTCCCTTCTTTATAGCTGGAGTAGATCCTAACAAACAACATCTGTCGAGCAGGGTAACCCCGGTGCCGTTGTTACTTCTGCACAACCTGAAAAAACGAGGGGATAGTCCCAACCAGGGCGATCATCGTGCTATGGGTC",
    #                     "TCAGTAGTTTCGAGTAAAGCCTCGCCTGGCCGGACGTGTGTGACCAGCACTCTAACTAGCTCCAATAAACATAGCGAGATGGTTGTCCGGCAAGCGGATGCTTAAACGTCAGTGGGTATTCTTCCGTCGTCCCGAGATAAAACGGTCCTAAGGCAA",
    #                     "CTGATAAAGTTCCTAGTGCTTAAAGCTGTCTTGACGCAACACAGGTCGTTCCACGTTTCAAGATGGGGGTATCCATTACTCGGTATTGCCCTTGGTAACTACAGGCGGATAGACCCTACAGTAGGAGGAAATCCGGGACCTCCACTTAATTCTGCA",
    #                     "TGAGGTCTGCTTTGACTCATTATATATCCATGAAGCCCTGCCGCGACGCAATTTGTAAGAATCTGTGTTATGATCTTCCGCTTTATGGTCGGGGAGGTAATATCCATACGGGCCGGCAAAGATCGCACCGTCGGTATAGCCCACTCGCCATCATGA",
    #                     "GGCCCTAGATTGAGCAGATTGAATTTTGCACTTAGCTGCGTTACTGAATTCAACCATAGGCGATGAGCGTCCACGTGCAGCATGTCTAAGGTTAATTCTAAGCTTAGCGAAAATCTTTATGCTGCCTTGACGATGAGTGTTCTCTCGTTTGTCAAA",
    #                     "CATGTGTTTATAAAATTTCAATCGGGCTATCAATCCACTGACTAGACGTTGCCCTCTACTTCACTCAGATCAGGAATTCCGAGGCGAGACTAAGAATGCGCAGATGTTGTCGCTTTAAGCTAGGACTCGCTCGGTAGCGCAATGCAGGCTCGTTGT",
    #                     "GAGGAGCTTCGAGGGCTCTGAAGATGGTAGTGGGCAGGTATGATAAATCCTGGCTCGACGTGCAGCAAGGGAGGTGCGATGACACTACGCATGTCCTGTCGTTCGCCCTAGGCGCGATTGCGTTGAGGTTCGAGCCGTTCTCCACCAACCTTTATG",
    #                     "TAGGACAAAGGGCAGTTACGATCCCAAGCGGCTTATGTGCAACTTACTCACTCATGCCATCGTACATAGCCTTGACCCCATGCCACTGACAAGACGGACTAAGAACTCTCCGTGTGCTAGAAATGCCACGTAGACTGCGTAAAATGATCGGGAGGC",
    #                     "ACTGGCCCGACGGGAAGCAGTTTATCTCTAAGTCGTTCACTACTATAGAATTTGTTGGGCTGCGTCGGATCGGAGGAAGATACTTAGATTCATAATACTCGTCCCGCCGAGAAAAAGGAGCAAGTCAACTCACCCAGAACAAGTCTTTGCGTGGAG",
    #                     "CCATAAAGGGCCAGTCAAAGACCATTTCAAAATACAGCTGTCTGGACGCACCAAGGAGTGAGTTCCCTTGCCGTGCTAGCGATGATCTGTTCCAAGATTACATAGCGGCCTATTATTAGGTCAGAGAGAACAGTCGCGTTTAGTGGTGGATTGAAA",
    #                     "GAATAGACCCTATTCAAATCCACGTTAGCACGAACCAGCAGGCGAATGTACGAATATGATTCAGGCCCCAATCACTACACCCGGAACAGACGGGGGCTATGTGTTCATGAGCCGTGAGCAATTCGTGTCGCGGGCAAACCCGGTCCTGACTCGACG",
    #                     "AACCCAAGCCACATCTAGTGCAGTTACTGCCTGTGTAGCGGTCTGCACTCCATGGGCTACAACGTATCATCGGCTGCCCAGAGATGGAGAGCAGGTAGTTATGCAACTTCTGTCGCGACGGCCCCGTTTGCCCCCCTGCAAACAACACGTAGGGAA",
    #                     "CCCCCAGGCTACCCTGACGGGACGAACCGTAGGACCAAGCGGCACAACAATATTCCAAGTAGACATAGGATGACACGGGAGTTACTTCATCGTGTACTAATTTGAAGGATCACCTGTCCTCTCACGGGACGGCAATAATGCCTGTCGAGAACAAGT",
    #                     "TACCAATAGAACGCCGACACCTGCTAAAATCCTCGTAAGCCGTAACTAGCTGCCTAGACGTCCTAGGCAGCAATTTATGAAGTAAATCATTCACAATCTGATTCTAGTAGGAACCGAACACCTGAGCCGCGGGTAACGCGAGTTGACACACCGGCG",
    #                     "CCCTCGTAGTCAAAGTACCTGCGAATCACTAAGACGTGTAGCATCAGCTCAGGGAGTATATCACCATAATGTGATCATCCCATGTGCCTTCCCTTTACTTTACTAAGTCGTTCTCTCTCGCTTTGCGTGAACCCTGCCGTGACGACGAACGAGCTG",
    #                     "AAACCTCTTTACGAGCGCCTATCGATCCAGCTGTATTCGTTATGTGGGAGTAATTTATCTGCTGGCCTGACGCGTCCTTCATAGAGACTTAGTACACCCAGGCTGCTAGTCTCGCGCCTTGCCTCTCAGCAGACGTGGCTACGGAGTCCAGCTCAA",
    #                     "CCTTCCTGAAGCCGAACGTTTCAGGCGAGTGACCCGTTATTCCATGACCTGTACTGGTATGGCCCTTCATTTCGCTGGTAGATCTTTCCACATTAACGCAAGTCGTTGTGCTTCGAGTCTAACACCGCACACGCTGGCCTGACGGACAGCTCGTAA",
    #                     "CGAGTCGCCAAGAGCCTGTGATTTGCGCACTTCCCCGCAAAAGACTGGCCTGACCCGACGGGCAGGCGGGATTTGACGCTGGTCCGCAGAAAGCTCAAGCGCTTAGACTTTCAATGCTAACAATATCCCCGATAGCTCGTGTCAGAGCCCCGGTAG",
    #                     "ATCGGAGGTCATGCCCGGTGACTTTCCGTGCGGATGAGTGGCTTATTAACAAGGTAACGAGCAGCCTGCACTTCTGCCGCGACGTTGTCGTTAGCAACTAACAAGGGATCGTACAAGAGACTCAGATGGAAATTTATACAGCTAATTATCGTATGA",
    #                     "CCCACTAAGAGGCAAGTTGGGTAACCTGATCGCTATACGGCAGGTTTATGCGATTTGACGTCTGACAAGACGCGAGGAATAGACCCAGCTAGTCCTTTGCGTCTCAATGATGACATGCCCCTCGTCGAGAGTAGAAGTGGTCCGTACACTGCAGTT",
    #                     "CCTGTCACGACGTTGAGCGAATATATTAATGTTTACGTGAATCTTGGGTCCAATCCAGTACCATTTCCCTGACGAATATGATCGTGCGTCCCGCCTAGTGCGTCCAGTATAAACTAAGAGCGATGACGTCAGGGCGCGACCGGGTCGTGAAATTGC",
    #                     "AGTAGGTACTCTAACAGCGCCCCAGTGTGATCCTTTTCGAAGGAAATTAACGCGCGGGACGTAGGAGTATGTATCTGAACGGCAGATTATATTATACGGACTAATTTAAGGAATGAATATAACACGGTAACGGCTGTCAAGACGTGGCGCCTCCAA",
    #                     "CTTTTGGGATCAATCTTCCACGAAGTCCGGCCCTTATTTCTGTTCACATCTACCCTCCAGCAACGTCTCGCCTTCCGGTTGGATGAATCGAAAGTCATATTGCAAGGAGTCAAGCGTGTGACTGACTGGACGGCTAACCGGGCTGGAGCTAGAGCC",
    #                     "ATTGCTTTTAAGAGAGTCACACGCTTGATAGAGGTTACTGGCACGACGGGACGTAACTCTCTCCTTCCTTAGCATTGCAACTACTGCATTCGCTCAAACACATTAGGTCGAGTTCGCCCAGTAGGCATTTGGGGCGGCTTTTCAAATCCCTGATCC",
    #                 ],
    #                 12,
    #                 25,
    #             ),
    #         )
    #     )
    # )

    # print(
    #     module_three.distance_between_pattern_and_strings(
    #         "CAATG",
    #         [
    #             "TCCATTTTTTGAGACGTGGTTCGAGGGAGGAAGGCTAGGACATAATTCCGATCAGGACATAGAAACGACACGTGAGTCTCACTCACAGTT",
    #             "ACACGCCTTACTTTTAGAATCTTCAGCGTACATTTTGCCACTATCGCACTCTCTGGAGCCACGTACGATAAAAAGCCGTACCTAGTGTCT",
    #             "AGGAGAAACCTGAGGCCTGCAGCCGTGACGATTTCAATTCTTCGCACACACGCTCGGCATACTGTAGACGGACGAACACTTTAGGAAGCT",
    #             "GCGCTATCGTCGAAAGCCTCGAACACGATCCGTTTTTTACATATATGTGGTCGGCGGGATAAGAGAACTACTATTGAGTTATGAAATCGG",
    #             "CCCTCATATAGCCTTATTATGTTAGCGAATGCTCCCATGAATAAAGAAACGGGGTGCGTAATTCATCAACTACTATTCTCCCGAGTTCGC",
    #             "TACTCCCCCTTTTTTCGTGCCAACACCGAAGGAAGCGCCAACGGTAACGACGACCATTCCTATACCCGGAGTCACGTTCCAAGAAGACAC",
    #             "GCAGATGTGGCAATGATCGGCTCTTGCTAGATCTGTCCAAATGAGTTTGCACCACGTTAACTTGGACCTATCCACGTGGGAGCACTAGCG",
    #             "CGACGGGCGTCCACCACGTCTCCACAAGTTAACTCCATACTCGCTCGGCCTTGGTTGATGTTCTATCGGGCCACGGGATCTTCAGGGGTA",
    #             "GTAAAGGCATCTGTAAGCCTTGGCGCTCATATTACCGAGCGCTCATTGATCAAGCCCGCTTAATCGCGGCGTCTGCGTAAGGAGCCTGAC",
    #             "CGCATCATGTCGCGTAAAAGCCCATTGCCATCAGCCCTTACATTGGTGCAGCAACTTGGTACACAGCAGATCAGCTGACAGCAAGTCTTC",
    #             "TGAATGCACCGAACAAGTGGGATGCTTGTCAATATCACAACATCTTGTCCCCCCAGCCGGGTAGGATGCGCCTATTGCCCCGTTGATTCG",
    #             "CCCGCGCCAATAGCTCGTTAATACCTAGCAGGGTCGAACAGCCGGTGCTCCGGGTATCGGTGCACATCAACATGAGCCCCCCAAACATTT",
    #             "GCTATGTTTGCGCTGTAATGTTCGTGATATGTCTGTCGCTCTTCTATGCGGGGGAATTATTTTCCATACCAGTAGTTATTTTTGCAGCCT",
    #             "TGCCCACACGGGATAGGTGGATTCTTTTTCCCCTTGACTGATTCTATCTGGTGACGCTCTACATGCCTGTTGCCTCAACGGGATCGTGCA",
    #             "TCGAATTGGAGATATACTCCTTAGACCCTTCATTTATATCGGAAGATATTGGGGTAACGAGCCTTGGTCTTGATGCGTAGCGAATGGACC",
    #             "TCGTAAATATCCCCTGGTTAAATTCCACCGTAGAAACTGGAATTTTCTGAGCACCATCGTCAGATACTAGGTGAGTGTGTCACAGAGCAT",
    #             "CTGCTCGTAGCAGGTCGGGTCTCGCTCCGCGGGCGTTTGCCATACCCATTTTAAGGAAGCTCGAGTTCCATCACGCCAGACGCCTAGCTC",
    #             "GAAACAACGGGTCTTGCCGGATCGGAGGTCCCCTAGTACTCTACAATGGACCCACTTGACCACGCGTTTTCACAGTCATAATCACGGGAA",
    #             "CTTGATAACAAGGGCATCCTTAAGATACACATCTCGTACTCCCTGGTTTTTACTGCGCATGGTCTACACTGTTGTTTTATTTCAGCAGAT",
    #             "CTTGCGAACCCTCACCTATCCGCTGTTACTATTGTCAGCCGTTCCCCCGCATCCCAGACCTTGTTGCAAAAGTCGGCCTCATCAGGCATA",
    #             "TGTTTAGAAATTCACTGTGTACTACGATATGCTGAGCCGGTTGACGCGACAACGATCCATAGGTTATCCAGGAGTGCCCTGAACCGGCCA",
    #             "TAATAAGTGTCCATATTACGACATCCTCGGATGCCTCAATGTACGACGGCCCACATCTCCGGTACTGTAACCCTGTACAGCCGTAACTAG",
    #             "GGACCCCCGGATGAACGCTATTACACTTTGCACGGACCGGCCTTTTGGTGAGGCCCACCATAGTCCATTAATGCACGGGCGGCCGATATA",
    #             "TTATACGATGAGTTGCTCCAACCCCAACTTAATGGATACTGATGCTTATGCCATGGAATGTTCACTCTTACCTTTGAGACTTATCCGGGG",
    #             "TTAACAATAGAGGTGTTGGACCGGCAATACCCAGTGACGACGACTCTCGTACGCCGCACGTCGTGGAGGACAGAGATGGCTATATAGGCC",
    #             "CTGACGGCAACATTTGGATCTGTCGTCATCGTGAAAGGACGCCGATCGAATCCCAATCAGCCAACAGCCAGTGCGCAATTTCACACTAGC",
    #             "TATCCTTCTAGGCTCGATCCGAGAAAGGGATATGTACACACCGTAGTTTTTCACCATATCGGAAGACCCAGCCCATAGTCAGCCGCAACG",
    #             "CCGTCGTGAGGGTAGTGGATGTCTCTACCCGCCTTTGCATCTTCTCACTACTAGGCCAGCAACTTTTCGCGCTTAACGCGTTAGCTTAGA",
    #             "TGGTAGATGACATGCGCGTATCAAACGTTGGGGTGTGACCGCAATCGTCTACCGGGAATGCCCGGACCAAACATTAAGTGACTCCTGGCT",
    #             "TGGGGGTAGCAGTGCGTTCGAGCCCGATGCAGTGTGTAACTAGCATGATACGAAACTGACGCCCGCAAAGCTGACGAGTGGAAGTATGAC",
    #             "TGGCTGTAAGACTTCGTCGAGACTGTTGTACGCTTTACAAGGCAGCTGACACTCCACAAGATTGAGGCTGGAAAGTCGCGCTTGAGGAGG",
    #             "TAGGACCGGTGAAGAGGGCTTTTGTCCGGTGATAGGGTGTAAAACATTTCTGTAGATCCTTTGTTGTTTGCGCCTTGGTATCCAATCAAG",
    #             "GGTTGTTGGACCAAGATATATGGTTGAAGTTTGTTTCGTTTCTACAAAGTATCGACCAGCCCCAGACAGAGTTAAAATTCAGCGGACGCA",
    #         ],
    #     )
    # )

    # print(
    #     "".join(
    #         map(
    #             str,
    #             module_three.median_string(
    #                 [
    #                     "CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC",
    #                     "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC",
    #                     "GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"
    #                 ],
    #                 7,
    #             ),
    #         )
    #     )
    # )

    profile = {
        "A": [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
        "C": [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
        "G": [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
        "T": [0.3, 0.1, 0.0, 0.4, 0.5, 0.0],
    }
    print(module_three.consensus_from_profile(profile))

    motif_matrix = [
        "CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC",
        "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC",
        "GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG",
    ]

    # 7-mer candidates to check
    candidates = ["GGTTACT", "AATCCTA", "AACGCTG", "ATAACGG", "CGTGTAA", "GAACCAC"]

    for candidate in candidates:
        distance = module_three.distance_between_pattern_and_strings(candidate, motif_matrix)
        print(f"Distance for {candidate}: {distance}")
