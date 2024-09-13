import re
from collections import defaultdict

""" Please note that this is my personal attempt at decoding the message and doesn't reflect the actual content of the bioinformatics course. 
    While the course includes a section on encoded messages and deciphering techniques, this specific example is purely for exploration and is unrelated to the biological focus of the course.
    You can find the official solution here: https://en.wikipedia.org/wiki/The_Gold-Bug. """


class GoldBug:
    def __init__(self, file_path: str):
        self.file_path = file_path
        self.encoded_message: str = ""
        self.patterns: dict = defaultdict(int)

    def find_patterns(self, min_length: int = 2, max_length: int = 4):
        with open(self.file_path, "r") as file:
            self.encoded_message = file.read()
            non_alpha_message = re.sub(r"[A-Z]", "", self.encoded_message)

            for length in range(min_length, max_length + 1):
                for i in range(len(non_alpha_message) - length + 1):
                    pattern = non_alpha_message[i : i + length]
                    self.patterns[pattern] += 1

            self.patterns = {pattern: count for pattern, count in self.patterns.items() if count > 1}
            # for pattern, count in sorted(
            #     self.patterns.items(), key=lambda item: item[1], reverse=True
            # ):
            #     print(f"Pattern: {pattern}, Count: {count}")

    def replace_patterns(self, pattern_mapping: dict[str, str]) -> str:
        decoded_message: str = self.encoded_message
        for pattern, replacement in pattern_mapping.items():
            decoded_message = decoded_message.replace(pattern, replacement)
        return decoded_message


if __name__ == "__main__":
    gold_bug = GoldBug(file_path="../resources/the_gold_bug.txt")
    gold_bug.find_patterns()
    char_mappings: dict[str, str] = {
        "6": "A",
        "·": "B",
        "5": "C",
        ")": "D",
        "‡": "F",
        "1": "G",
        "(": "I",
        "9": "J",
        "0": "K",
        "†": "L",
        "3": "M",
        "2": "N",
        ":": "O",
        "?": "P",
        "—": "Q",
        ";": "R",
    }
    print(gold_bug.replace_patterns(char_mappings))
