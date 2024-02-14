import RNA
import pandas as pd


def features(gene_sequence, start_codon_index, stop_codon_index):
    start_codon = "AUG"
    stop_codon = ["UAA", "UAG", "UGA"]

    if gene_sequence[start_codon_index-1:start_codon_index + 2] == start_codon and \
            gene_sequence[stop_codon_index-1:stop_codon_index + 2] in stop_codon:
        cds_sequence = gene_sequence[start_codon_index + 2:stop_codon_index-1]
    else:
        print("Invalid start or stop codon index.")
        return None

    gene_length = len(cds_sequence)//3
    if len(cds_sequence) % 3 != 0:
        print("Please make sure coding sequence is in triplets")
    else:
        gene_length = gene_length+2
        five_prime_utr = gene_sequence[:start_codon_index-1]
        length_of_5prime_utr = len(five_prime_utr)
        folding_energy_80 = calculate_folding_energy_80(
            gene_sequence, start_codon_index, length_of_5prime_utr)
        folding_energy_70 = calculate_folding_energy_70(gene_sequence)
        koz_score = kozak_score(gene_sequence, start_codon_index)
        in_frame_AUG = calculate_in_frame_AUG(five_prime_utr)

        encoding = {"A": 1, "U": 2, "G": 3, "C": 4}

        N1 = encoding.get(gene_sequence[start_codon_index - 6])
        N4 = encoding.get(gene_sequence[start_codon_index - 3])

        gene_features = {
            "gene_length": gene_length,
            "folding_energy_70": folding_energy_70,
            "folding_energy_80": folding_energy_80,
            "length_of_5prime_utr": length_of_5prime_utr,
            "kozak_score": koz_score,
            "N1": N1,
            "N4": N4,
            "in_frame AUG": in_frame_AUG
        }
    return gene_features


def calculate_folding_energy_70(gene_sequence):
    sequence_70 = gene_sequence[:70]
    (ss, mfe) = RNA.fold(sequence_70)
    return float("{:.2f}".format(mfe))


def calculate_folding_energy_80(gene_sequence, start_codon_index, length_of_5prime_utr):
    if length_of_5prime_utr <= 40:
        sequence_80 = gene_sequence[:start_codon_index + 42]
    else:
        sequence_80 = gene_sequence[start_codon_index -
                                    41:start_codon_index + 42]
    (ss, mfe) = RNA.fold(sequence_80)
    return float("{:.2f}".format(mfe))


def kozak_score(gene_sequence, start_codon_index):
    koz = gene_sequence[(start_codon_index-1)-6:(start_codon_index-1)] + \
        gene_sequence[(start_codon_index-1)+3:(start_codon_index-1)+6]
    score = 0

    if len(koz) < 9:
        return 0

    if koz[0] == "A" or koz[0] == "U":
        score += 1
    if koz[1] == "A":
        score += 1
    if koz[2] == "A" or koz[2] == "C":
        score += 1
    if koz[3] == "A":
        score += 1
    if koz[4] == "A" or koz[2] == "C":
        score += 1
    if koz[5] == "A":
        score += 1
    if koz[6] == "U":
        score += 1
    if koz[7] == "C":
        score += 1
    if koz[8] == "U" or koz[8] == "C":
        score += 1

    return score


def calculate_in_frame_AUG(five_prime_utr):
    num_in_frame_aug = 0
    for i in range(0, len(five_prime_utr), 3):
        if five_prime_utr[i:i + 3] == "AUG":
            num_in_frame_aug += 1
    return num_in_frame_aug


def main():

    gene_sequence = input("Enter the gene sequence: ")
    five_prime_utr = gene_sequence[:17]
    start_codon_index = len(five_prime_utr)+1
    stop_codon_index = len(gene_sequence[:-2])
    gene_features = features(
        gene_sequence, start_codon_index, stop_codon_index)
    df = pd.DataFrame([gene_features], columns=["gene_length", "folding_energy_70",
                      "folding_energy_80", "length_of_5prime_utr", "kozak_score", "N1", "N4", "in_frame AUG"])
    print(df)


if __name__ == "__main__":
    main()
