import numpy as np 
import pandas as pd
import features
import pickle
import warnings
warnings.filterwarnings("ignore")


def evaluate_model(model, X_test):
    # Perform the prediction using the provided ML model
    y_pred = model.predict(X_test)
    return y_pred

def InitiationRate(gene_sequence, start_codon_index, stop_codon_index):
    gene_features = features.features(gene_sequence, int(start_codon_index), int(stop_codon_index))
    # Create a dataset from the calculated features
    dataset = pd.DataFrame([gene_features], columns=["gene_length", "folding_energy_70", "folding_energy_80", "length_of_5prime_utr", "kozak_score", "N1", "N4", "in_frame AUG"])

    # Load the ML model
    rf_model = pickle.load(open("tir_rf_model.pkl", "rb"))
    # Perform initiation rate prediction
    initiation_rate = evaluate_model(rf_model, dataset)
    return gene_features, round(initiation_rate[0], 4)

def main():
    """
    gene_sequence = "CACCAGGUUUUUGGCUUUUUAGAUUUUAUCCCCUUCCAGCAUGAGGAUUGGCACGGAUGCUAACGUGAUAAUCUGGCUGUAG"
    start_codon_index = 41
    stop_codon_index = 80
    """
    #gene_sequence = input("Enter the gene sequence: ")
    #start_codon_index = int(input("Enter the start codon index: "))
    #stop_codon_index = int(input("Enter the stop codon index: "))
    
    gene_sequence = input("Enter the gene sequence: ")
    five_prime_utr = gene_sequence[:17]
    start_codon_index = len(five_prime_utr)+1
    stop_codon_index = len(gene_sequence[:-2])
    


    rate = InitiationRate(gene_sequence, start_codon_index, stop_codon_index)
    print(rate)

if __name__ == "__main__":
    main()
