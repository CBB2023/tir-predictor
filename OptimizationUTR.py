import numpy as np
import pandas as pd
import random
import InitiationRate
import warnings
warnings.filterwarnings("ignore")


def OptimizationUTR(gene_sequence, start_codon_index, stop_codon_index, target_initiation_rate, iterations, kbt_values: list):
    seq = gene_sequence
    tir = target_initiation_rate

    # Initialize parameters
    s = list(seq)
    s_old = s
    nuc = ["A", "G", "C", "U"]
    ul = len(seq[:start_codon_index - 2])
    data = []
    I_old = round(InitiationRate.InitiationRate(seq, start_codon_index, stop_codon_index)[1], 4)
    I = I_old
    # kbt_values = [1e-01, 1e-02, 1e-03, 1e-03, 1e-03, 1e-04, 1e-04, 1e-04, 1e-05, 1e-05]
    index = 0

    # Optimization loop
    for i in range(iterations):
        if i % (iterations // 10) == 0 and index < len(kbt_values):
                kbt = kbt_values[index]
                index += 1
        r1 = random.randint(0, ul-1)
        r2 = random.randint(0, 3)

        s = list(s_old)
        s[r1] = nuc[r2]  # mutated string
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
    print(df)
    return df 

def main():

    gene_sequence = input("Enter the gene sequence: ")
    five_prime_utr = gene_sequence[:17]
    start_codon_index = len(five_prime_utr)+1
    stop_codon_index = len(gene_sequence[:-2])
    target_initiation_rate = 0.08
    iterations = 10000
    opt = OptimizationUTR(gene_sequence, start_codon_index, stop_codon_index, target_initiation_rate, iterations)
    print(opt)

# if __name__ == "__main__":
#     main()

