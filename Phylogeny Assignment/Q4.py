import os
import numpy as np
import parameters as params

def juke_cantor(seq1, seq2):
    # initiating mismatch counter
    mismatch = 0

    # iterating over each pair
    for (x,y) in zip(seq1,seq2):
        if x != y:
            mismatch += 1

    # calculating the prbability of each base
    beta_val = 0.75
    prob = float(mismatch) / float(max(len(seq1),len(seq2)))
    
    #calculating distance
    distance = -beta_val * np.log(1-prob/beta_val)
    
    return distance


def kimura(seq1, seq2):
    # initating beta, gamma, delta counters
    beta = 0
    gamma = 0
    delta= 0

    # iterating over all pairs
    for (x,y) in zip(seq1,seq2):
        if x+y in params.trans:
            beta += 1
        elif x+y in params.transv1:
            gamma += 1
        elif x+y in params.transv2:
            delta+= 1

    # calculating the beta, gamma and delta values
    beta = float(beta)/float(max(len(seq1),len(seq2)))
    gamma = float(gamma)/float(max(len(seq1),len(seq2)))
    delta = float(delta)/float(max(len(seq1),len(seq2)))

    # calculating the distances
    kimura3_dist = -0.25 * np.log((1 - 2*beta - 2*gamma) * (1 - 2*beta - 2*delta) * (1 - 2*gamma - 2*delta))
    kimura2_dist = -0.5 * np.log((1 - 2*beta - (gamma + delta))*np.sqrt(1 - 2*(gamma + delta)))

    return kimura2_dist, kimura3_dist


# Main loop
if __name__ == "__main__":

    seq1_path = input("\nEnter path to first sequence:\n")
    seq2_path = input("\nEnter path to second sequence:\n")

    if os.path.isfile(seq1_path) and os.path.isfile(seq2_path):
        with open(seq1_path, 'r') as file:
            seq1 = file.read().replace('\n', '')
        with open(seq2_path, 'r') as file:
            seq2 = file.read().replace('\n', '')

        print("\nJukes-Cantor distance", juke_cantor(seq1, seq2), sep='\n', end='\n')
        
        kimura2_dist, kimura3_dist = kimura(seq1, seq2)
        print("\nKimura2 distance", kimura2_dist, "\nKimura3 distance", kimura3_dist,sep='\n', end='\n')
    
    else:
        print("Invalid input path")