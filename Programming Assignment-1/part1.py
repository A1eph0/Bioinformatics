import parameters # importing necessary macros

def reverse_complement_finder(strand):
# Gives reverse complement of given DNA sequence

    flip_strand = strand[::-1]
    output = ""

    for nucleotide in flip_strand:
        if nucleotide == "a":
            output +="t";
        elif nucleotide == "t":
            output += "a";
        elif nucleotide == "c":
            output += "g";
        elif nucleotide == "g":
            output += "c";
        else:
            print("Invalid input string\n")
            exit()

    return output


def rna_finder(strand):
# Gives the RNA transcribed by the coding DNA (forward DNA) 

    output = strand.replace("t", "u")

    return output


def amino_finder(strand, in_rna = ""):
# Gives amino acid sequence translated by the given genome
    
    rna = rna_finder(strand) + in_rna
    output = ""

    index = rna.find("aug")
    if index == -1:
        return output
    
    
    output += "Met"
    index += 5

    while index < len(rna):
        read = rna[index-2:index+1]

        if read == "uuu" or read == "uuc":
            output += "-Phe"
        elif read == "uua" or read == "uug" or read == "cuu" or read == "cuc" or read == "cua" or read == "cug":
            output += "-Leu"
        elif read == "auu" or read == "auc" or read == "aua":
            output += "-Ile"
        elif read == "aug":
            output += "-Met"
        elif read == "guu" or read == "guc" or read == "gua" or read == "gug":
            output += "-Val"
        elif read == "ucu" or read == "ucc" or read == "uca" or read == "ucg":
            output += "-Ser"
        elif read == "ccu" or read == "ccc" or read == "cca" or read == "ccg":
            output += "-Pro"
        elif read == "acu" or read == "acc" or read == "aca" or read == "acg":
            output += "-Thr"
        elif read == "gcu" or read == "gcc" or read == "gca" or read == "gcg":
            output += "-Ala"
        elif read == "uau" or read == "uac":
            output += "-Tyr"
        elif read == "uag" or read == "uaa" or read == "uga":
            break
        elif read == "cau" or read == "cac":
            output += "-His"
        elif read == "caa" or read == "cag":
            output += "-Gln"
        elif read == "aau" or read == "aac":
            output += "-Asn"
        elif read == "aaa" or read == "aag":
            output += "-Lys"
        elif read == "gau" or read == "gac":
            output += "-Asp"
        elif read == "gaa" or read == "gag":
            output += "-Glu"
        elif read == "ugu" or read == "ugc":
            output += "-Cys"
        elif read == "ugg":
            output += "-Trp"
        elif read == "cgu" or read == "cgc" or read == "cga" or read == "cgg" or read == "aga" or read == "agg":
            output += "-Arg"
        elif read == "agu" or read == "agc":
            output += "-Ser"   
        elif read == "ggu" or read == "ggc" or read == "gga" or read == "ggg":
            output += "-Gly"
        
        else:
            print("error:", read)
            exit()
        
        index += 3
    
    
    if index+1 < len(rna):
        remains = rna[index+1:]
        output+="\n" + amino_finder("", in_rna=remains)

    return output

def trunc_amino_finder(strand, in_rna = ""):
# Gives amino acid sequence translated by the given genome

    rna = rna_finder(strand) + in_rna
    output = ""

    index = rna.find("aug")
    if index == -1:
        return output
    
    
    output += "Met"
    start_occurence = index + 2 # occurence of start codon
    threshold = parameters.truncation_threshold # threshold for Met occurence (set in parameters.py)
    index += 5

    while index < len(rna):
        read = rna[index-2:index+1]

        if read == "uuu" or read == "uuc":
            output += "-Phe"
        elif read == "uua" or read == "uug" or read == "cuu" or read == "cuc" or read == "cua" or read == "cug":
            output += "-Leu"
        elif read == "auu" or read == "auc" or read == "aua":
            output += "-Ile"
        elif read == "aug":
            if index <= start_occurence + threshold*3: # accounting for truncation if within threshold
                output += " (truncated) \nMet"
                start_occurence = index # reassigning occurence of start codon
            else:
                output += "-Met"
        elif read == "guu" or read == "guc" or read == "gua" or read == "gug":
            output += "-Val"
        elif read == "ucu" or read == "ucc" or read == "uca" or read == "ucg":
            output += "-Ser"
        elif read == "ccu" or read == "ccc" or read == "cca" or read == "ccg":
            output += "-Pro"
        elif read == "acu" or read == "acc" or read == "aca" or read == "acg":
            output += "-Thr"
        elif read == "gcu" or read == "gcc" or read == "gca" or read == "gcg":
            output += "-Ala"
        elif read == "uau" or read == "uac":
            output += "-Tyr"
        elif read == "uag" or read == "uaa" or read == "uga":
            break
        elif read == "cau" or read == "cac":
            output += "-His"
        elif read == "caa" or read == "cag":
            output += "-Gln"
        elif read == "aau" or read == "aac":
            output += "-Asn"
        elif read == "aaa" or read == "aag":
            output += "-Lys"
        elif read == "gau" or read == "gac":
            output += "-Asp"
        elif read == "gaa" or read == "gag":
            output += "-Glu"
        elif read == "ugu" or read == "ugc":
            output += "-Cys"
        elif read == "ugg":
            output += "-Trp"
        elif read == "cgu" or read == "cgc" or read == "cga" or read == "cgg" or read == "aga" or read == "agg":
            output += "-Arg"
        elif read == "agu" or read == "agc":
            output += "-Ser"   
        elif read == "ggu" or read == "ggc" or read == "gga" or read == "ggg":
            output += "-Gly"
        
        else:
            print("error:", read)
            exit()
        
        array.append(output)
        index += 3

    output += " (functional)"

    if index+1 < len(rna):
        remains = rna[index+1:]
        output+="\n" + trunc_amino_finder("", in_rna=remains)

    return output


dna_strand_unedited = input("Enter DNA strand in 5' -> 3' order:\n")
dna_strand = dna_strand_unedited.lower()

print("\nReverse strand for the DNA strand is:",reverse_complement_finder(dna_strand), sep='\n', end='\n')
print("\nRNA transcribed by the DNA strand is:",rna_finder(dna_strand), sep='\n', end='\n')
print("\nAmino acid (assuming no truncation) synthesised by the DNA strand is:",amino_finder(dna_strand), sep='\n', end='\n')
print("\nAmino acids (assuming truncation) synthesised by the DNA strand is:",trunc_amino_finder(dna_strand), sep='\n', end='\n')


