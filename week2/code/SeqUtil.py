import csv
from collections import Counter
from heapq import nlargest
import operator

def dna2rna(dna):
    Dna2RnaMap = {'A':'A', 'T':'U', 'G':'G','C':'C'}
    rna = "".join([Dna2RnaMap[nucleotide] for nucleotide in dna])
    return rna
	
def reverseComplement(dna):
    complementMap = {'A':'T', 'T':'A', 'G':'C','C':'G'}
    complement = "".join([complementMap[nucleotide] for nucleotide in dna])
    reverseComplement = complement[::-1]
    return reverseComplement

def rna2peptide(rna, codonMap=None, mapFile=None):
    if codonMap==None and mapFile==None:
        raise ValueError('Either codonMap or mapFile must be provided')

    if codonMap==None:
        with open(mapFile) as f:
            freader = csv.reader(f, delimiter=' ')
            codonMap = {row[0]:row[1] for row in freader}
    
    codons = [rna[start:(start+3)] for start in range(0, len(rna)-2, 3)]
    peptide = "".join(["-" if codonMap[c]=="" else codonMap[c] for c in codons])
    peptides = peptide.split(sep="-")
    return peptides

def loadAminoAcidMass(file):
    with open(file) as f:
        freader = csv.reader(f, delimiter=' ')
        massMap = {row[0]:int(row[1]) for row in freader}
    return(massMap)

def massOfPeptide(peptide, AAMassMap):
    mass = [AAMassMap[p] for p in peptide]
    totalMass = sum(mass)
    return totalMass

def cyclospectrum(peptide, AAMassMap):
    lenPeptide = len(peptide)
    peptide2 = peptide * 2   # 2 copies of peptide to handle circular subpeptide

    spectrum = [0]          # spectrum storage
    for l in range(1,lenPeptide):     # loop for subpeptide length
        subpeptides = [peptide2[pos:pos+l] for pos in range(0, lenPeptide)]
        masses = [massOfPeptide(p, AAMassMap) for p in subpeptides]
        spectrum.extend(masses)

    spectrum.append(massOfPeptide(peptide, AAMassMap))
    spectrum = sorted(spectrum)
    return spectrum

def cyclospectrum2(peptideMassSeq):
    lenPeptide = len(peptideMassSeq)
    peptide2 = peptideMassSeq * 2   # 2 copies of peptideMassSeq to handle circular subpeptide

    spectrum = [0]          # spectrum storage
    for l in range(1,lenPeptide):     # loop for subpeptide length
        subpeptides = [peptide2[pos:pos+l] for pos in range(0, lenPeptide)]
        masses = [sum(p) for p in subpeptides]
        spectrum.extend(masses)

    spectrum.append(sum(peptideMassSeq))
    spectrum = sorted(spectrum)
    return spectrum

def linearspectrum2(peptideMassSeq):
    lenPeptide = len(peptideMassSeq)
    peptide2 = peptideMassSeq

    spectrum = [0]          # spectrum storage
    for l in range(1,lenPeptide+1):     # loop for subpeptide length
        subpeptides = [peptide2[pos:pos+l] for pos in range(0, lenPeptide-l+1)]
        masses = [sum(p) for p in subpeptides]
        spectrum.extend(masses)

    spectrum = sorted(spectrum)
    return spectrum

    
def cyclopeptideSequencing(spectrum, AAMasses):
    peptides = []  # candidate peptides
    sequences = []  # sequence storage

    spectrum = sorted(spectrum) # ensure that spectrum is sorted
    AAMasses = set(AAMasses)    # ensure unique mass

    spectrumSet = set(spectrum) # will be required later to calculate consistency
    
    # expand function
    def expand(peptides, AAMasses):
        expandedPeptides = []
        if len(peptides)==0:
            expandedPeptides = [[aa] for aa in AAMasses]
        else:
            expandedPeptides = [x+[y] for x in peptides for y in AAMasses]
        return(expandedPeptides)

    def isConsistent(peptide, spectrumSet):
        #return sum(peptide) in spectrumSet
        
        lenPeptide = len(peptide)
        for l in range(1,lenPeptide+1):   # loop for subpeptide length
            subpeptides = [peptide[pos:pos+l] for pos in range(0, lenPeptide-l+1)]
            masses = [sum(p) for p in subpeptides]
            # inconsistent if all masses are not present 
            massSet = set(masses)
            if((massSet & spectrumSet)!= massSet):
                return False

        return True
        
            
    
    parentMass = max(spectrum)  # parent mass is the max mass in spectrum

    loopCount = 0
    while True:
        loopCount = loopCount + 1
        
        peptides = expand(peptides, AAMasses)
        #print(loopCount, ":", len(peptides))
            
        for p in peptides[:]:
            mass = sum(p)
            if(mass == parentMass):         # p's mass = total mass
                if(cyclospectrum2(p)==spectrum):
                    sequences.append(p)
                peptides.remove(p)
            elif isConsistent(p, spectrumSet)==False:    # p is not consistent with spectrum
                peptides.remove(p)

        if len(peptides)==0:
            break

    return sequences


def peptideSpectrumScore(peptide, spectrum, cyclic=False):
    if cyclic==True:
        peptideSpectrum = cyclospectrum2(peptide)
    else:
        peptideSpectrum = linearspectrum2(peptide)

    peptideCounter = Counter(peptideSpectrum)
    spectrumCounter = Counter(spectrum)
    matchCount = [min([peptideCounter[k], spectrumCounter[k]]) for k in peptideCounter.keys()]
    score = sum(matchCount)
    return score
    

def leaderboardCyclopeptideSequencing(spectrum, N, AAMasses):
    leaderboard = []    # candidate peptides
    leaderPeptide = []  # sequence storage
    leaderScore = 0     # score of leader peptide

    spectrum = sorted(spectrum) # ensure that spectrum is sorted
    AAMasses = set(AAMasses)    # ensure unique mass

    # expand function
    def expand(peptides, AAMasses):
        expandedPeptides = []
        if len(peptides)==0:
            expandedPeptides = [[aa] for aa in AAMasses]
        else:
            expandedPeptides = [x+[y] for x in peptides for y in AAMasses]
        return(expandedPeptides)

    
    parentMass = max(spectrum)  # parent mass is the max mass in spectrum

    while True:
        # branch
        leaderboard = expand(leaderboard, AAMasses)
        
        # bound
        for p in leaderboard[:]:
            mass = sum(p)
            if mass == parentMass:         # p's mass = total mass
                score = peptideSpectrumScore(p, spectrum, cyclic=True)
                if score > leaderScore:
                    leaderPeptide = p
                    leaderScore = score
            elif mass > parentMass:
                leaderboard.remove(p)

        # if leaderboard is empty, break the loop
        if len(leaderboard)==0:
            break

        # calculate each peptide score and cut
        scores = [peptideSpectrumScore(p, spectrum, cyclic=True) for p in leaderboard]
        nthScore = min(nlargest(N,scores))
        leaderboard = [item[0] for item in list(zip(leaderboard, scores)) if item[1]>=nthScore]

    return leaderPeptide


def spectralConvolution(spectrum):
    lenSpectrum = len(spectrum)
    conv = [abs(spectrum[i]-spectrum[j]) for i in range(0, lenSpectrum-1) for j in range(i+1, lenSpectrum) if spectrum[i]!=spectrum[j]]
    return conv


def convolutionCyclopeptideSequencing(spectrum, M, N):
    '''
    M is the number of most-frequent amino acids with weight in [57, 200]
    N is the number of top scored peptides in leaderboard
    '''

    conv = spectralConvolution(spectrum)
    conv = [c for c in conv if 57<=c<=200]
    convCounter = Counter(conv)
    nthMultiplicity = min(nlargest(M, list(convCounter.values())))
    AAMasses = [mass for mass in convCounter if convCounter[mass] >= nthMultiplicity]
    peptideSeq = leaderboardCyclopeptideSequencing(spectrum, N, AAMasses)
    return peptideSeq

    
