import numpy as np

# seq1 = input('Seq 1:') or "GTATCCAACGGTTGTGTGAGTAAAATTCTGGGCAGGTATTACGAGACTGGCTCCATCAGA" # Pax6 Mouse gene
# seq2 = input('Seq 2:') or "GTGTCCAACGGTTGTGTCAGTAAAATCCTGGGCAGATACTATGAAACAGGATCCATCAGA" # Shark eye control
seq1 = "GTATCCAACGGTTGTGTGAGTAAAATTCTGGGCAGGTATTACGAGACTGGCTCCATCAGA" # Pax6 Mouse gene
seq2 = "GTGTCCAACGGTTGTGTCAGTAAAATCCTGGGCAGATACTATGAAACAGGATCCATCAGA" # Shark eye control

mismatch = 0
for pair in zip(seq1,seq2):
    if pair[0] != pair[1]:
        mismatch += 1
b = 0.75
p = float(mismatch) / len(seq1)
jc_d = -b * np.log(1-p/b)
print('Jukes-Cantor Distance:',jc_d)

transitions = [ "AG", "GA", "CT", "TC"] #beta,P
transversions1 = [ "AT", "TA", "GC", "CG"] #gamma,Q
transversions2 = [ "AC", "CA", "GT", "TG" ] #delta,R

P = 0
Q = 0
R = 0
for (x,y) in zip(seq1,seq2):
    if x+y in transitions:
        P += 1
    elif x+y in transversions1:
        Q += 1
    elif x+y in transversions2:
        R += 1

P /= len(seq1)
Q /= len(seq1)
R /= len(seq1)

k3_d = -0.25 * np.log((1 - 2*P - 2*Q) * (1 - 2*P - 2*R) * (1 - 2*Q - 2*R))

Q1 = Q + R
k2_d = -0.5 * np.log((1 - 2*P - Q1)*np.sqrt(1 - 2*Q1))

print('Kimura-2 Distance:', k2_d)
print('Kimura-3 Distance:', k3_d)