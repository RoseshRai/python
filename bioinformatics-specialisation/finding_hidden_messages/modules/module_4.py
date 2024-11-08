import random


class Module4:
    @staticmethod
    def random_motifs(dna: list[str], k: int, t: int) -> list[str]:
        """Randomly select k-mers from each DNA sequence."""
        motifs = []
        for i in range(t):
            start = random.randint(0, len(dna[i]) - k)
            motifs.append(dna[i][start : start + k])
        return motifs

    @staticmethod
    def random_motifs_v2(dna: list[str], k: int) -> list[str]:
        """Randomly select k-mers from each DNA sequence."""
        motifs = []
        for seq in dna:
            start = random.randint(0, len(seq) - k)
            motifs.append(seq[start : start + k])
        return motifs

    @staticmethod
    def profile_with_pseudocounts(motifs: list[str], k: int) -> dict[str, list[float]]:
        """Construct a profile matrix with pseudocounts."""
        t = len(motifs)
        profile = {"A": [1.0] * k, "C": [1.0] * k, "G": [1.0] * k, "T": [1.0] * k}

        for motif in motifs:
            for j in range(k):
                profile[motif[j]][j] += 1

        for nucleotide in profile:
            for j in range(k):
                profile[nucleotide][j] /= t + 4

        return profile

    @staticmethod
    def profile_with_pseudocounts_v2(motifs: list[str]) -> dict[str, list[float]]:
        """Construct a profile matrix with pseudocounts."""
        k = len(motifs[0])
        profile = {nucleotide: [1.0] * k for nucleotide in "ACGT"}

        for motif in motifs:
            for index, nucleotide in enumerate(motif):
                profile[nucleotide][index] += 1.0

        # Normalize the counts to probabilities
        for index in range(k):
            total = sum(profile[nucleotide][index] for nucleotide in "ACGT")
            for nucleotide in "ACGT":
                profile[nucleotide][index] /= total

        return profile

    @staticmethod
    def profile_most_probable_kmer(text: str, k: int, profile: dict[str, list[float]]) -> str:
        """Find the most probable k-mer in text based on the profile matrix."""
        max_prob = -1.0
        most_prob_kmer = text[:k]

        for i in range(len(text) - k + 1):
            kmer = text[i : i + k]
            prob = 1.0
            for j in range(k):
                prob *= profile[kmer[j]][j]
            if prob > max_prob:
                max_prob = prob
                most_prob_kmer = kmer

        return most_prob_kmer

    def generate_motifs(self, profile: dict[str, list[float]], dna: list[str], k: int) -> list[str]:
        """Generate motifs from the profile matrix for each DNA sequence."""
        motifs = []
        for seq in dna:
            motifs.append(self.profile_most_probable_kmer(seq, k, profile))
        return motifs

    @staticmethod
    def score(motifs: list[str]) -> int:
        """Compute the score of the motifs based on the consensus sequence."""
        consensus = ""
        k = len(motifs[0])
        for j in range(k):
            count = {"A": 0, "C": 0, "G": 0, "T": 0}
            for motif in motifs:
                count[motif[j]] += 1
            consensus += max(count, key=lambda nucleotide: count[nucleotide])

        score = 0
        for motif in motifs:
            for i in range(k):
                if motif[i] != consensus[i]:
                    score += 1
        return score

    def randomized_motif_search(self, dna: list[str], k: int, t: int) -> list[str]:
        """Perform randomized motif search algorithm."""
        motifs = self.random_motifs(dna, k, t)
        best_motifs = motifs
        while True:
            profile = self.profile_with_pseudocounts(motifs, k)
            motifs = self.generate_motifs(profile, dna, k)
            if self.score(motifs) < self.score(best_motifs):
                best_motifs = motifs
            else:
                return best_motifs

    def run_randomized_motif_search(self, dna: list[str], k: int, t: int, iterations: int = 1000) -> list[str]:
        """Run randomized motif search multiple times to find the best motifs."""
        best_motifs = self.randomized_motif_search(dna, k, t)
        best_score = self.score(best_motifs)

        for _ in range(iterations - 1):
            motifs = self.randomized_motif_search(dna, k, t)
            current_score = self.score(motifs)
            if current_score < best_score:
                best_motifs = motifs
                best_score = current_score

        return best_motifs

    @staticmethod
    def compute_probability() -> float:
        """Compute the probability of capturing at least one implanted 15-mer."""
        # Parameters
        total_positions = 586  # Total number of possible 15-mers
        total_trials = 10  # Number of randomly selected 15-mers
        successful_positions = 585  # Number of positions where the implanted 15-mer does not appear

        # Probability of not capturing the implanted 15-mer in one trial
        p_not_capturing_one = successful_positions / total_positions

        # Probability of not capturing the implanted 15-mer in all trials
        p_not_capturing_all = p_not_capturing_one**total_trials

        # Probability of capturing at least one implanted 15-mer
        p_capturing_at_least_one = 1 - p_not_capturing_all

        return p_capturing_at_least_one

    @staticmethod
    def probability_at_least_one() -> float:
        """Compute the probability of capturing at least one implanted 15-mer."""
        total_positions = 586
        p_not_capture_single = (total_positions - 1) / total_positions

        # Probability of capturing at least one implanted 15-mer
        return 1 - p_not_capture_single**10

    @staticmethod
    def probability_exactly_one() -> float:
        """Compute the probability of capturing exactly one implanted 15-mer."""
        total_positions = 586
        p_capture_single = 1 / total_positions
        p_not_capture_single = (total_positions - 1) / total_positions

        # Probability of capturing the 15-mer on exactly one specific string and not on others
        p_exactly_one = 10 * p_capture_single * (p_not_capture_single**9)

        return p_exactly_one

    def compute_probability_at_least_two(self) -> float:
        """Compute the probability of capturing at least two implanted 15-mers."""
        p_at_least_one = self.probability_at_least_one()
        p_exactly_one = self.probability_exactly_one()

        # Probability of capturing at least two implanted 15-mers
        return p_at_least_one - p_exactly_one

    @staticmethod
    def profile_randomly_generated_kmer(text: str, k: int, profile: dict[str, list[float]]) -> str:
        """Randomly select a k-mer from text based on the profile probabilities."""
        probabilities = []
        kmers = []
        for i in range(len(text) - k + 1):
            kmer = text[i : i + k]
            prob = 1.0
            for j in range(k):
                prob *= profile[kmer[j]][j]
            probabilities.append(prob)
            kmers.append(kmer)

        total_prob = sum(probabilities)
        if total_prob == 0:
            probabilities = [1 / len(probabilities)] * len(probabilities)
        else:
            probabilities = [p / total_prob for p in probabilities]

        kmer = random.choices(kmers, weights=probabilities, k=1)[0]
        return kmer

    def gibbs_sampler(self, dna: list[str], k: int, t: int, n: int) -> list[str]:
        """Perform Gibbs sampling algorithm."""
        random.seed()  # You can set a specific seed here if you want reproducible results
        motifs = self.random_motifs_v2(dna, k)
        best_motifs = motifs[:]
        best_score = Module4.score(best_motifs)

        for _ in range(n):
            i = random.randint(0, t - 1)  # Randomly choose an index
            motifs_excluding_i = [motifs[j] for j in range(t) if j != i]
            profile = self.profile_with_pseudocounts_v2(motifs_excluding_i)
            motifs[i] = self.profile_randomly_generated_kmer(dna[i], k, profile)
            current_score = self.score(motifs)
            if current_score < best_score:
                best_motifs = motifs[:]
                best_score = current_score
        return best_motifs

    def gibbs_sampler_v2(self, dna: list[str], k: int, t: int, n: int) -> list[str]:
        """Perform Gibbs sampling algorithm."""
        random.seed(0)  # Set a fixed random seed for reproducibility
        motifs = self.random_motifs_v2(dna, k)
        best_motifs = motifs[:]
        best_score = Module4.score(best_motifs)

        for _ in range(n):
            i = random.randint(0, t - 1)  # Randomly choose an index
            motifs_excluding_i = [motifs[j] for j in range(t) if j != i]
            profile = self.profile_with_pseudocounts_v2(motifs_excluding_i)
            motifs[i] = self.profile_randomly_generated_kmer(dna[i], k, profile)
            current_score = self.score(motifs)
            if current_score < best_score:
                best_motifs = motifs[:]
                best_score = current_score
        return best_motifs


if __name__ == "__main__":
    module_four = Module4()
    # k = 15
    # t = 20
    # dna_input = [
    #     "TGGGGGTTATCTCAGGGGACTCAGCGCACCTGGGCGTACAATCTGAAAGACGTGGTCATAGCCATGATGTCCATCCCTTATACACACTGGCAATGCCTACCGGCGATTAGTACGCTAAGAGGGGGACACACCAGGCGCACTGAGCGGTGACGATTGCTAATTATGGGGGTTATCTCAG",
    #     "GGGACTCAGCGCACCTGGGCGTACAATCTGAGAAGTTTGCCCTGGTAAGACGTGGTCATAGCCATGATGTCCATCCCTTATACACACTGGCAATGCCTACCGGCGATTAGTACGCTAAGAGGGGGACACACCAGGCGCACTGAGCGGTGACGATTGCTAATTATGGGGGTTATCTCAG",
    #     "CCCCTTCTCTCGAACTGTCGGCAAATGTCTACATGTTCGTTCATAAAACCCACAGGACGGGGATATACCACTAACATGCCGCAATTATGGTAGCAGCGCAAGCCCCTGACAAAAAACGGCCTGCACCTTTTAAAAGTGGACAATAAAGTCCTTGGTCCCGAACAGTGCCCTCAGGGGC",
    #     "ACTCAGGGTGCCGGGGGCCGCCCCTACGCGCAATTGTAATCCCCTTGTTTGAATGCGGGGTCCTAGTACATGAACAGTGCCCCCCTGCTGTTTCCAAGGTTCTATTTAACAAATCGGCACCCATCCACCTCAGCGTATACCTACTTTGGTAACATCTAAGGAGCCACAACTTGTCCGT",
    #     "AGGCCTAATAACTGTGCGAACAGTGCCACCGTACGCACTAATCGCTAGCAGTTGGCCATCGACGAGCCCCCAACCGCAAGTACACCGTGGATGCCGACGCCCTCGGGAGCCTCGAGGGTGAAGGCTCTCGAGCTTGGACACCGTTAAGGTTCTTGGCATGTCCATGCTTTTGGCTCAG",
    #     "TAGGCCGATATCGAATAGTATATATGCACCCCCCTGAAAAAAATACCCCTCGGGGCTCAGATGGCTAAAAGGGACGAGCCTCATGGGAATGCTGCCCTGGTAGGTCTGTATCGACAAGGCAGAGTAAAAATTGTGTACCAGCAAACGACTCTGTAATCGCCATCTAATGAGGCTTCGG",
    #     "AAAAAGTAGTCCTTGCGGCCACATGAGGTACCTTGCCGTTGACGGCTTCAGAATTCGTGACGTTAAGCGATTGTGATTCCCAGGTACATCATAGTTGATACGAAACACAAGACTAACATGGCTCGTGACATTGAACCTCCCAGTGCCCTGGTGAGAGGTTCAGCAGCCATCCCGCAAT",
    #     "GGCCGGCGAGTGCCCTGGTGTCAGGGCCTCGACCTACAATCATGGAAACGCGTAGTTAGTAAAATGACGGGTAGAGAGAGGGGATATATCATGCTATTGAACCTGTCGGCAGGCTACTGACGACGGCAAACAATAGGGTATGTCTGATAACCCCAGAGTGGCCGAATGAAACCGCAGT",
    #     "TCTATCTCAAGTCAACTATTGAAGGGAACTACGCCCTGGTCAATTTGGTTACCGACTTGTACCACCGATTTCTCAACGTGCCCCCCGCGACGCCGATCAGCGTCAATTATCGAAGTGGGGTTCACGTTTTTACATACTGTAAACCTCCGGACGGCCCACGCTCCTGTCCGGGAGCCAC",
    #     "CTAACTACCATTAAACTCGGATAGTTACTCATACTAGATAATCGTTTAAGAGCTGAGTGAAGCGGACGCCCGTGGGACGCCGCCGGTGCCGGCTCTCTGTGAAGCCGGAACGCTGTCGAACAGTGTTGTGGTGCAACGCCGTGTGGAACAGTTCCGGCGGAAGCAGGTTGACGGGATA",
    #     "AGCCTTGCGGATAGCGGTGTTCGCCAGTATGGATGCACGCGTAGTTGAACGGCTAGTGCCCATATCGTAAGATGAACACCGGAGGCTTGGTTCACGGAATCAACGGAGCCGGATAAAAATTTACTGTAGAACAGTGCGAGGGTTTTGCTTTTTTAGGCATGGGTCGAAGTGGGAATAA",
    #     "ACCGTACTACACTACACATAGGTTGGGCTAGATTCAAATGGAGTGCCCGCCACGGCTCTTGAGGCCGGTGTTTACATGAGCTGTTCAGCGCAGGATATATTTGAACAGGAACCTGGTTCATTGTAGTGGAAATATGACGGCTCTATCCAGACCGCCCTATTTGGTGCGCTACGAACTT",
    #     "AAGCATCCATTCTCACCAAACAAATGCTTATCATGCTAGGAAGTTCGGAAGGGGGCAGAGAGGAACCTGGCCCTGGTCCAGAGTATCGCGGGCTCAAACGCACTTGGGTAGCCGTGGTGATACTCACTCGTTTCCAGAGCAGCGCAAACATTTTGGCACAGTTTATGATAGTTAGCCA",
    #     "TAACAGTGCCCTGCCCAAAAGCTAGTCGCCTACCGTTGTTATATGGCCCAGTTCGAGATTGTCTGTGCAAACATATGGTATCGGCGCACGTTAATCATTCCCGAGGAGCTGACTGAGTTGTGAGAATGCCCTCGACAAAGCACGGCCCGGAACCTCTATAGGCAAGAGTAGATCACCT",
    #     "ATTGGGTCTAGCTTCATGGGATCGTTTGGCGCCGGTCCTGTCCTGATTAGGCTATCAGACGAAAGTCAGGGCTCTCAGCGCTATAATATGCGTCCGGGAACAGTTATCTGGTGAGTGATCGTCTTTTGGTACCGCCCCCGCTCGTATTGAGTCTTATGATGTGGGATCAGGGATCTAT",
    #     "CCCGAGATCCGGTTTAAGTGCTGATATAACGCATATGAGAGTAGTGTGAGAACATGCCCCTGGTGTTGAACACCGCATAATCGCAGTGTGATCCGTCATATGGTTTATCAAAACCGCCGTCGGTACTGCACAAGGTCGTCCGAGAACCCTAGATCTGATCAGCGCCCGCAGCTTCGAT",
    #     "TGCCTCGCGAACAGCAGCCTGGTCCGTTTTTGCGCGGCTGCGGGCTCCAACAACCAAAACATGGCAGATAGTTCTTACTACGAGCCTTATTAGTTTGTCCAGGGACGACAATTATTATATTAAGGTAGAGATTTTTAACGTGAGCTAGTAGTCACGTGACTAGCGTGGGAAGACCATC",
    #     "GGTGGGGCTAGTGGGGACCCGAACGGGTTGTGCGTACGTCAATGCAACTACAAAGTCTCAGAGACGTGCCCTGGTTTTGTCTGTTGAGGTGACGTAGCCAGCTGGCTACCGCTCAGCCACTACTTACGTCTTTAGGCCTTCCGGCCACGTGACAACAGAATGGCTGGATGTGGTGGCT",
    #     "CCCGCGGTTGAGGCCTGATACATGCTGCAAAGCCGGTTTCTTGTTGCTCAGAGTCTATCCCGCCATTGAGGCACCACAGTGCCCTGGCTGGGGAGTACCTCTCGCTGCATCGTGTTGGCGTCTAGAGTACCAGACGCTTTCGGCCTATACTTTCAACTAACGAAACTGACGATTAATG",
    #     "TTATGTACACGAAATGAAAAAATTGTCCGCCAGCACCAGGCTCCAAGACACCTTATATGCCAGATGGGACTTCAACCTATAATTTCCCCTTTATTCCCTTTAATGAGAGCCTTGTTGGCCAAGAGCGGGACTCTTGTTTAGACGGTCTGCGTAGAACACCACCCTGGTCTCGTGAAGA",
    # ]
    #
    # best_motifs_out = module_four.run_randomized_motif_search(dna_input, k, t, 1000)
    # for motif_iter in best_motifs_out:
    #     print(motif_iter)
    #
    # # print(module_four.compute_probability())
    # print(module_four.compute_probability_at_least_two())

    # input_dna = [
    #     "TTCTTACGAAGGGAGTTAAAATTGAAGCCGTGCTCAAGATGAGACAATCAGTCAATAACCCGATAGTCGTACCCGCTAAAAAGGGAGACAATCTCGTCAAATCGGCCTTACACCTTGCTGATGACGTATTGTGTACGGGATTTTAATAGGCCTAACAAGTCGCTTGGTAGGACTCCTGGATCCCGGAGTAAACCTTTCGTTGGGCTGTCTTGCTAAGCGCGGCTGAAGTTGTCATTCCGAGCTCACTGGCCGTCGCTGAGGGGTAATAACTCGTGTAACGGAGCGCCTCGCATAAGTAGGTATTTCTTACGAAGGGAG",
    #     "TTAAAATTGAAGCCGTGCTCAAGATGAGACAATCAGTCAATAACCCGATAGTCGTACCCGCTAAAAAGGGAGACAATCTCGTCAAATCGGCCTTACACCTTGCTGATGACGTATTGTGTACGGGATTTTAATAGGCCTAACAAGTCGCTATGTTCTGCTAAGTATGGTAGGACTCCTGGATCCCGGAGTAAACCTTTCGTTGGGCTGTCTTGCTAAGCGCGGCTGAAGTTGTCATTCCGAGCTCACTGGCCGTCGCTGAGGGGTAATAACTCGTGTAACGGAGCGCCTCGCATAAGTAGGTATTTCTTACGAAGGGAG",
    #     "TTACTCCGTTTATGCGAAAATAAGTAGAATCACGTAGAAGTTTAGGCCCCTTCACACCCAGTGAGTTCGAGAACTGTTATGGACAAAGTGACACTTTACACTTTCTTATTTTAGCCACTCGAGCCAGTAGAGCCTAGGACCCCTTGAATTAAGCAATAGGACTCAAATTAGGCAAGAAAAACGCGGTCAGACGGTTAAGGACGCTTAGTCACTGTTACTGACTACCCACCGTTGCCAGCATTGACGTAGTCTCTGCAGCATGGTTTTGAAACGGGATGCTGTGTCCTGTTCTACCGTCAGTGATGAACACTCGCCTTT",
    #     "CTGCAAGTCACAGTTAGGCCGGATAGGTCGAATGGTGTTCTAAGCGTTAGTTGCAGCCGCGGTAAGTGGAATTGCCCTGATTCACAGGCCACCCATCCGACTGGAATAATATGCTAAGTATTTCACAACGAGGGACCGTTGTGTTGCTTGACCAATTTTGCGTTAATTCCTCAATGAGCTTCAGGCGCGACCAACCTACGACTTGTCAATTGCTACTACCGCGACGATCTGATAAATAGGATTGTGATTGGACACCTCGCTGTGCACATGATATGAGTGAACCATTGACAGATCACAAGCAGGTAGCTCCATTGCCTA",
    #     "GGCCTGACTAGGTAAATCTTCGACCCTCGGGGGCGAGATAGTAATCGAAAAAGCAACGCATTACACGAAGGGAGTAGAGCCGTGAACTGGTTAGAGCATCCTGATGTAGCCGACTTACCCGGAACAATGCGATGCATTGTAGTCCACACTGTAAGGGAACCTACAACGACTCAGACGTCCAATTGCCTTATGCGGAGTTGCCGTGTGCCTGCGCTCGCTTCCTCCACTGTGTCTTCAGCTAGTTCCCCACATACGAGTACGCATAAAGAGTCAATTCGGTTATCTAATTGATCCATATGGCTGAGGTCACCGCAGCTC",
    #     "GACTTGTTAGAAATCATTTATAAACTAACTACTTTAGCGCCGTGGCGGAAGGCGGAATTCCTATTGTCCCGCATGGGCACCGTGTCGACCTCCACTGCATGTCATGTGCGAGCCGACCCAGGTTAGCAGACGGTAGCGCCTACTTATTCCATATCTGTTATAGTAGTTTCGGCCCCCCTCTGTTAACCAAAATGGACGAGGCTTGACCTATTCGTTATCAATATGACGCATCCTACGGGGCTGGTCGTTTAAGCGTGAAATGCGATCGGAATATCCCTCCGGAACTAGCTCGAGATGCGATGCTTTATAACAAGACCT",
    #     "ACTCTCCTGAGGATTAGAATGGACGAGGACCTTGCGCTCGCACACAAGAGGGGTCTGTTAGGATTGGGTTCCTTAGAAATCTCCCGATATTTAAAGGTGCGGGGTATCGCTTTCCCACATGCATGCCTCATCTGACGGTCTGAACGCACACTAGAACGAGTGCTCGTTGTGAGTTTCGCAGTCCGCCACCTTCTCGATCCCATAAACACCCCCCTGAAGTTGCGAACGAAGTACTCGATGGGAGTATCGGCTGCACTCTCTCGTTCTAGATATGCGATGCTAATGCCGTACATCAGCTGAATGCTATGTACGGGCATA",
    #     "TATGTATTGATAGGCTAGAGTTGCTGAACGTAACTGTAGATTAAGAAGTTGGCGTTACGGATCCTTCCAGCACGGTGGTACCCAAGTCGGAACGATTATAACTCTCAAGAGTTATCGACAGGCTAACCCGGTTTTAATGGCTCGGGGGCGGAGGCTGTGTGAGCATAAAGAATGCGCATCTAAGTATCGTTAGCTCACCAGGATAGTAGAAAATTATACAGCGTCCCTAGGGTCGAGCCATAAATATGATGAGATCTACCCTTAACATCATAGACGAGAAGTGACAGATAGTAAGGTCGCGCGCAGGTCCCGATGCTA",
    #     "ACCGTGGGTTGTTACCGTAGCGGATGCTAAATTGCAGCTCCAGGTGCACAAAGCATAGGTTAGATGTAGTGTACTTATGCCGCGCTAAGTACCCAGACAATAGGCTTAGTCGGTCGCCTTGTCACTGGTTTAGCTCAAAACGTTTTCTTCAGTTCCATTATTCGCGCGGGCATGCCGTGTGGCTGAGTCCTGGCGCATACCTAGAGATTCCACGTAAATAACGTCCCCCCTATAAGCCTGCGCATGACAAAAGATTCCCTCACTTACGTGCTAAGCGGTTTTTAATTATTAAGGGACTAGATTAGGGCCCACATACTA",
    #     "GGAAGCGACTGAGGCGTTGAGGCGTCTTAGGACTACAACACCATTAATCTAAGAGGCCCTCAATGCGATTACAAGTAGAGTTCCTCTCGTTACGTCCTGCTGTGGCTCAAGTGGAAGTGGTCTTGTCGAAACGATCGACCGGTGAGGGGTAAGCACGTGGGTTGCTCGAGAGAATCGGTAACACTTCCGTTATTTCATTACCGTGAAAGGTACCGCTTTCTCACGGAGGCGGCATGAAGGCGAATCTCTTCTAGGCGGGAATTCACGGGATCGAGTTCTTGGTTCACTAAGATTCCGGACTACATACTACAAAAAGCA",
    #     "CTTTACATTCTACACAAACATATCGGAGTCGCAACAGAGGGCCTAACCTATTTCAGAGTATGTCATAGAAGAATCAAGTTAGATAAACTCTCGCTACGGGGTATCCAGATGCTTCTACTCTGACCCAATGCGATGCTATAGACGGCACCCATGCGTCTCGCACAATCGTCGAAGGAGCTCTGTGTCCCACTACTCTAGCCGGAATTGGGTATCCTGCCTAATGTGATAAGCTTGGGAACCCCCCGACCAGACCCCGCACGTGACCGCTAAAAAACAGAGGGTCGGAGAATAGAGAGTCAAGCTGCCCACCACCAGCCG",
    #     "TTTTGTCCGCTGGATTAAGTGATTGAAAACCGAGCGCACCACGCCGTCCGACCTCATAATCATCGGAAGGATGGTGCCGTTGCTACCCGGCGCGATGCTAAGTTTGTGCGCACGGTCCTCCCGAGTCAACCGCAAGGAAACACATCTTAGGCTATCCTGCGACCGTTCAACACGCTATACCCACACGCGCGGGTGCTTGGCAACTCGGGGCTCTTCATCAATTGTCCGCCTGGGAGGCACATCAGCTAATAACTCTGTTTGATATGCAAAGTAGGTTGGATTTGAGTGCACCTTTATAAGTTGTTCGGTTCCGGTCCA",
    #     "CGACCCCTAAATGGTGCGCCACGGCTCGTCATTCCGCTCATCGATTGTTATTGAAATCCATTTGTGGTTAACGCCGATGACAGTAATGGACACCGCGTTTCAATCATTGTTTTACCGTTAGTACTGACGTCTCAACAGGAAAGCCGGGAGCCCAATCTCCCGGATACGCTCCGTTTTCCGAGAGGGAAGATGCGATGTCGAGTATCACAGACCGGTACGCTGACGTTACTAAGAGAGGAGACCGAACTGCACCTATGCCCAATGCTGGCGATGTGGCATGTGGGGTAGGAGCGCACACCTTAAAGGCAACATCGCTTC",
    #     "GTTGGAACCCCGAGGCATCCACGGCTCGGGTAACGTCATCCTAATAACCGTGGCATCAGCACGTACCGAACGGCGTAAGGGATTGAAAGCGGATAGCACTATGCTGGGGATACATCGCCATCCCAAGAAGGACGTGCGGCTGCTTCCGACGTTGGGGGCCGTTGTTCGCAAGGCAGCCCCGCGCGGAATAGGATTGGCTCACCTAGCGGAACATGCGGCACTAAGTAGATGAGTTGTGTACATTGCCGCTAGGATCCACGGAGGGACGGCGATGGAAAACTTAATGGGTAAAGCTATATACAACGGCGGTATCTATCG",
    #     "CGATGTACTCAGAAGGGAGAGGATGCTAAGTACTTGTCGGAGCGAGGTGCTGGGCCCGTATCATGCCAGCGTCACTTGCGCGTGTGAGGGGTACTGACCGGACAGACACATGCCGTTATAAAAAATAAGAATACATGCGATGTGCCTGCGACGCGGCCTCGGGGTGCGAGCCTCCCGGTAGGCCCAAAGACGACGCGGGGTACGTCTAGATTCGGGCGAGCGCGTTTGAGGCTGGTGCGCCAGAACGCCATGCCACCGACATCCGAGGCTGGTAAGAGCATATTTACGGCAGTCGTCGTGCGCTTTACCTTACCCTCG",
    #     "TCGATCTGTATGTTTGCAGCCACTTCTGACGCCTCTCCAGAATTATCCCCATGAAAGCCAAGCAAAATAATTCCGCGTTGATGATGACCCGGTGCGTCTAGAGAGTGTGCGAACCCACTCTGTGTGGATACTTGGTTTATAGAACGTAGGATTGCAGTCCCATGCGAACTTAAGTAACTTCCTGCCGAATCAGATGCTATCCTTACTGTCTACCCTAGAGACTGTGAAAGCGTAAGTAAGCACCAGGTGAGCCACCTACCTTCTGGGCAGGGACCATCTCGCGGTCAGGCACCAGGCATCTGGATCCAGAACACTCCC",
    #     "ATGTCGTTGACATATCCACTCAAAAGCGCGGACTTTAGAATCCGTCTGGCCCCGAGAGGTGTGTTACTGAGCCTCCCGACGATAAGCGCCCCCGTTATTACGCGCTTATTCTCGCGCAGATTCCCCTCCGGCTAGGTAGACACGCAGACCTTGGGGGGTGCAGCCGAAAACCGATGCGCTGAGGGGCAAGGAAATAGGGCTCCGAGGCAGAGAACCTTTGTGTTTATTCCATAGCAAGCGTACCATCTGCGATGCTAAGAGCACTCACTAGAGCCAACGTCATCGCCATACGAGCTATCCTTACGTCCACGCATTTCC",
    #     "AATACAGATCCATCATGACAATGCGCCCTAGTGGGCAATCCCTACTACCAGCCAAAGCCGATCTGGGAAGGCGACGGGTCGAGTTGCGTACTCTCTCCCTACAGATCCGCGACTAGTCGTCGATTAATATAGGGTGCGGAGGAGAGGTCAGCCGGTGGGTCTAGTCTCCTGTGTTCAGCTTAGTCGACGTCACCCCGGGTCCTTTCCCTCCCATCCTGTACGCATCATGCCTGGCTAAGTAGAACATCAGTTGATAATACAGACAAGAGTGTTATTTGGAATCGTGTTTTAAACTGTCCGTCACGGCAACTAATCACC",
    #     "CGTCGCGAATGTTGGAATTGCTAATTTAACTACCAGGACTAGTGGTGGGGCCTTGACACGATGCATGTCCTGCTAAGTATGCCCTGAGGAGAGGTGAGCCAGCGAGCTTGTGCCGTTACTTAATGGTCTTAAAGATTCATTAGATCCAAGTTAATACTTGCAGTGGCTCACAATCCTCTCTGATCTTGGCCTTAGTACCGAGCATGAATTCGTGGGAAGAGCTGTACCATTTGCTGTAAAAGTTGTGCTTTGACGATTTCCTTTTGCGAATCCAACGTAAGGTGGACGGTTATTCCATGGTGTTATTGGCGCGAATTG",
    #     "GGGCACTACGGAGTGCGAGGCTGTCGGCTATGGTTCGGGTTCGCGACGCCGGCAGGCCGGGTAATGATTAAAACCTGAGTAGAGAGCCGTTGGCATCAATGGGGTACGTGTGTTCTGAAAGCAGCTTAACACCTACCCATATGACGCCTATCGATTTCATAATCACCTTCCCGCAATTACGCACTCGGCAAGATGTTGGTGCAGTTAAATGTTTACGCAGTATTCAGTTCTTTGACTCTCGATGCTAAGTAGTAAATGCTACTCGTTGGCGCTATAACAGATTGATTGCACCTAACTAGGATGTGAGAAGGCGCAGTG"
    # ]
    # k = 15
    # t = 20
    # n = 200  # Number of iterations
    # runs = 20
    #
    # best_motifs = module_four.gibbs_sampler_v2(input_dna, k, t, n)
    # for motif in best_motifs:
    #     print(motif)
