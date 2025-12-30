# HorribleGeneAnalyzer
A horrible gene analyzer, still working on it 

# Example Usage

```python
from HorribleGeneAnalyzer import AnalyzeDNA, AnalyzeORFs

'''
    Horrible Gene Analyzer has two objects, both take a fasta file (genetic file) as a parameter

    Each object can then parse the file and build bitarrays for bases and codons in the AnalyzeDNA object, and orfs in AnalyzeORFs

    The objects can then return different outputs, as listed at the bottom

'''

# Load a file with python 
with open('file.fasta', 'f') as f:

  # Passes a fasta file as parameter
  bases = AnalyzeDNA(f)
  orfs = AnalyzeORFs(f)

  # Builds bitarray in memory
  bases.basebuild() # Can also do bases.reversebuild() for a reverse complement strand
  bases.aminosbuild()
  orfs.orfsbuild()

  pass

### Properties

# Accesses bases, codons and amino acids based of an index, like 10th and 11th, or with a slice 10:100 for 10 to 100 
bases.bases[index]
bases.codons[index]
bases.aminos[index]

# Accesses ORFs (potential genes) with the same indexing
orfs.orfs[index]

# Methods

# Returns ratio with the inputted fraction
bases.baseratio(form='at/gc')

# Returns the codon frequencies as a pandas series
bases.codonfreq()

# Returns the frequency of amino acids as a pandas series
bases.aminofreq()

```


