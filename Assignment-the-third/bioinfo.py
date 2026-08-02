#!/usr/bin/env python

# Author: Hannah Balmer hbalmer@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.4"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set('ATGCNatcgn')
RNA_bases = set('AUGCNaucgn')

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score as an integer'''
    return ord(letter) - 33


def qual_score(phred_score: str) -> float:
    '''Takes a Phred+33 score quality line as a string character and outputs the average quality score as a float value.'''
    total_P = 0
    for P in phred_score:
        total_P += convert_phred(P)
    
    return total_P / len(phred_score)


def validate_base_seq(seq,RNAflag=False):
    '''This function takes a sequence string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    seq = seq.upper()
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)


def gc_content(DNA: str) -> float:
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1. Case insensitive.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used a DNA sequence?"
    
    DNA = DNA.upper()
    Gs = DNA.count("G")
    Cs = DNA.count("C")
    return (Gs+Cs)/len(DNA)

def calc_median(lst: list) -> float:
    """Calculates the median of a sorted list, outputting it as a float. calc_median() function takes a list as its only argument."""
    if len(lst)%2 == 1: #If the length of the list is an odd number
        middle = int((len(lst) - 1) / 2)
        median = lst[middle]
    else: #If the length of the list is an even number
        middle = int(len(lst) / 2)
        median = (lst[middle-1] + lst[middle]) / 2

    return float(median)

def oneline_fasta(fafile: str):
    '''Takes a FASTA file and outputs a new FASTA file with all sequence lines stripped of newlines so that they are only one line.
 	Function will output new FASTA file with .formatted suffix.'''
    with open(fafile, "r") as fa, open(f"{fafile}.formatted", "w") as newfa:
        for i, line in enumerate(fa):
            if line.startswith(">") == True and i == 0:
                newfa.write(f"{line}")
            elif line.startswith(">") == True and i != 0:
                newfa.write(f"\n{line}")
            else:
                newfa.write(f"{line.strip()}")


if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    assert qual_score("RTA$C$$") == 25.0, "wrong average quality score for RTA$C$$"
    assert qual_score(">=;;") == 27.250, "wrong average quality score for >=;;"
    assert qual_score("+))012F") == 15.857142857142857, "wrong average quality score for +))012F"
    print("qual_score function is working")

    assert validate_base_seq("ATTCGATAG") == True, "Incorrect validation of DNA string ATTCGATAG"
    assert validate_base_seq("AUGACGU") == False, "RNA flag not passed in, AUGACGU did not pass DNA validation"
    assert validate_base_seq("AUGACGU", RNAflag=True) == True, "Incorrect validation of RNA string"
    assert validate_base_seq("ATDTTGCS") == False, "Non-DNA bases not detected"
    assert validate_base_seq("atgtgacta") == True, "Function is not case insensitive"
    print("validate_base_seq function is working")

    assert gc_content("GCGCGGGC") == 1.0, "GC content not calculated correctly for GCGCGGGC"
    assert gc_content("GCTAGTCA") == 0.5, "GC content not calculated correctly for GCTAGTCA"
    assert gc_content("gctagtca") == 0.5, "Function is not case insensitive"
    print("gc_content function is working")

    assert calc_median([10,16,25,45,95]) == 25.0, "Incorrect median for [25,45,10,16,95]"
    assert calc_median([1,2,3,4]) == 2.5, "Incorrect median for [1,2,3,4]"
    assert calc_median([1,1,1,1,1,1,1,2,4,6]) == 1, "Incorrect median for [1,1,1,2]"
    assert calc_median([4.6,7.9,8.5,56.2,67.0]) == 8.5, "Incorrect median for [56.2,4.6,8.5,7.9,67.0]"
    print("calc_median function is working")