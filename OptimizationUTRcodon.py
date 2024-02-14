import numpy as np
import pandas as pd
import random
import features
import InitiationRate
import warnings
warnings.filterwarnings("ignore")


def OptimizationUTRcodon(gene_sequence, five_prime_utr, start_codon_index, stop_codon_index, target_initiation_rate, iterations, kbt_values:list):
    seq_final = gene_sequence
    tir = target_initiation_rate
    utr = five_prime_utr
    sequence = seq_final[len(utr) + 1:]

    
    # Sequence to be used for UTR and codon mutation 
    if len(utr) < 70:
        sequence_70 = seq_final[:start_codon_index] + seq_final[start_codon_index + 2:start_codon_index + (71 - len(utr))]
    else:
        sequence_70 = seq_final[start_codon_index - 70:start_codon_index]

    len(sequence_70)

    if len(utr) >= 40:
        sequence_80 = seq_final[start_codon_index - 40:start_codon_index] + seq_final[start_codon_index + 2:start_codon_index + 42]
        net_seq = utr + seq_final[start_codon_index + 2:start_codon_index + 42]

    elif len(utr) < 40:
        if len(utr) <= 30:
            sequence_80 = utr+seq_final[start_codon_index + 2:start_codon_index + 42]
            net_seq = sequence_70

    elif len(utr) > 30 and len(utr) < 40:
        sequence_80 = utr + seq_final[start_codon_index + 2:start_codon_index + 42]
        net_seq = utr + seq_final[start_codon_index + 2:start_codon_index + 42]

    index_utr = len(utr)
    utr_mut = net_seq[:index_utr]
    cdsmut = net_seq[index_utr + 1:]
    rest_cds_seq = sequence[len(cdsmut) + 2:]
    if len(cdsmut) % 3 == 0:
        cds_mut = cdsmut
    elif len(cdsmut) % 3 == 1:
        cds_mut = cdsmut[:-1]
    elif len(cdsmut) % 3 == 2:
        cds_mut = cdsmut[:-2]

    seq = utr_mut + cds_mut
    
    
    # Dictionary of non synonymous codons
    syn_codons = {
        "GCU": ["GCC", "GCA", "GCG"],
        "GCC": ["GCU", "GCA", "GCG"],
        "GCA": ["GCC", "GCU", "GCG"],
        "GCG": ["GCC", "GCA", "GCU"],
        "CGU": ["CGC", "CGA", "CGG", "AGA", "AGG"],
        "CGC": ["CGU", "CGA", "CGG", "AGA", "AGG"],
        "CGA": ["CGU", "CGC", "CGG", "AGA", "AGG"],
        "CGG": ["CGU", "CGC", "CGA", "AGA", "AGG"],
        "AGA": ["CGU", "CGC", "CGA", "CGG", "AGG"],
        "AGG": ["CGU", "CGC", "CGA", "CGG", "AGA"],
        "AAU": ["AAC"],
        "AAC": ["AAU"],
        "GAU": ["GAC"],
        "GAC": ["GAU"],
        "UGU": ["UGC"],
        "UGC": ["UGU"],
        "GAA": ["GAG"],
        "GAG": ["GAA"],
        "CAA": ["CAG"],
        "CAG": ["CAA"],
        "GGU": ["GGC", "GGA", "GGG"],
        "GGC": ["GGU", "GGA", "GGG"],
        "GGA": ["GGU", "GGC", "GGG"],
        "GGG": ["GGU", "GGC", "GGA"],
        "CAU": ["CAC"],
        "CAC": ["CAU"],
        "AUU": ["AUC", "AUA"],
        "AUC": ["AUU", "AUA"],
        "AUA": ["AUU", "AUC"],
        "UUA": ["UUG", "CUU", "CUC", "CUA", "CUG"],
        "UUG": ["UUA", "CUU", "CUC", "CUA", "CUG"],
        "CUU": ["UUA", "UUG", "CUC", "CUA", "CUG"],
        "CUC": ["UUA", "UUG", "CUU", "CUA", "CUG"],
        "CUA": ["UUA", "UUG", "CUU", "CUC", "CUG"],
        "CUG": ["UUA", "UUG", "CUU", "CUC", "CUA"],
        "AAA": ["AAG"],
        "AAG": ["AAA"],
        "AUG": ["AUG"],
        "UUU": ["UUC"],
        "UUC": ["UUU"],
        "CCU": ["CCC", "CCA", "CCG"],
        "CCC": ["CCU", "CCA", "CCG"],
        "CCA": ["CCU", "CCC", "CCG"],
        "CCG": ["CCU", "CCC", "CCA"],
        "UCU": ["UCC", "UCA", "UCG", "AGU", "AGC"],
        "UCC": ["UCU", "UCA", "UCG", "AGU", "AGC"],
        "UCA": ["UCU", "UCC", "UCG", "AGU", "AGC"],
        "UCG": ["UCU", "UCC", "UCA", "AGU", "AGC"],
        "AGU": ["UCU", "UCC", "UCA", "UCG", "AGC"],
        "AGC": ["UCU", "UCC", "UCA", "UCG", "AGU"],
        "ACU": ["ACC", "ACA", "ACG"],
        "ACC": ["ACU", "ACA", "ACG"],
        "ACA": ["ACU", "ACC", "ACG"],
        "ACG": ["ACU", "ACC", "ACA"],
        "UGG": ["UCG"],
        "UAU": ["UAC"],
        "UAC": ["UAU"],
        "GUA": ["GUC", "GUU", "GUG"],
        "GUC": ["GUA", "GUU", "GUG"],
        "GUU": ["GUA", "GUC", "GUG"],
        "GUG": ["GUA", "GUC", "GUU"],
        "UGA": ["UGA"],
        "UAA": ["UAG","UGA"],
        "UAG": ["UAA","UGA"],
        "UGA": ["UAA","UAG"]
    }
    

    # Initialize parameters
    s = list(seq)
    s_old = s
    nuc = ["A", "G", "C", "U"]
    met = ["A", "U", "G"]
    list_utr = list(utr_mut)
    list_cds = list(cds_mut)
    list_mut_cod = [cds_mut[i:i + 3] for i in range(0, len(cds_mut), 3)]
    data = []
    I_old = round(InitiationRate.InitiationRate(seq_final, start_codon_index, stop_codon_index)[1], 4)
    I = I_old
    # kbt_values = [1e-01, 1e-02, 1e-03, 1e-03, 1e-03, 1e-04, 1e-04, 1e-04, 1e-05, 1e-05]
    index = 0
    s = list(seq)
    s_old = s
    
    # Optimization loop
    for i in range(iterations):
        if i % (iterations // 10) == 0 and index < len(kbt_values):
                kbt = kbt_values[index]
                index += 1

        r1 = random.randint(0, (len(utr_mut) + int(len(cds_mut) / 3)))
        r2 = random.randint(0, 3)
        if r1 < len(utr_mut):
            list_utr[r1] = nuc[r2]
        elif r1 > len(utr_mut):
            r3 = r1 - len(utr_mut)
            syn_cod = syn_codons[list_mut_cod[r3 - 1]]
            r4 = random.randint(0, len(syn_cod) - 1)
            list_mut_cod[r3 - 1] = syn_cod[r4]

        s = list_utr + met + list_mut_cod + list(rest_cds_seq)

        gene = ''.join(s)

        I_new = round(InitiationRate.InitiationRate(gene, start_codon_index, stop_codon_index)[1], 4)

        E1 = round(abs(tir - I), 4)
        E2 = round(abs(tir - I_new), 4)

        delE = E2 - E1

        if delE >= 0:
            p = np.exp(-(delE / kbt))
            r3 = random.uniform(0, 1)
            if r3 < p:
                I_old = I
                I = I_new
                s_old = s
                s = gene
            else:
                I = I
                I_old = I_old
                s = s
                s_old = s_old
        else:
            I_old = I
            I = I_new
            s_old = s
            s = gene
        
        data.append([tir, I, gene])

    df = pd.DataFrame(data, columns=['tir','I','gene'])   
    df = df.tail(1)
    return df

    
def main():

    gene_sequence = input("Enter the gene sequence: ")
    five_prime_utr = gene_sequence[:17]
    start_codon_index = len(five_prime_utr)+1
    stop_codon_index = len(gene_sequence[:-2])
    target_initiation_rate = 0.11
    iterations = 10000
    opt = OptimizationUTRcodon(gene_sequence, five_prime_utr, start_codon_index, stop_codon_index, target_initiation_rate, iterations)
    print(opt)

if __name__ == "__main__":
    main()

