# ALIGNer: Advanced Large-scale Identification of Genomic ENtries

The goal of ALIGNer is to develop an efficient tool for translating DNA sequencing reads into proteins and aligning them across all six possible reading frames. Using dynamic programming and a custom alignment algorithm, the tool will compare DNA sequences against known protein databases. Designed with scalability in mind, GENE-AIM leverages parallel processing and distributed computing to handle large metagenomic datasets, providing a high-confidence filtering system for rapid identification of relevant protein matches. This tool aims to enhance functional annotations and comparative genomics in complex microbial communities, making it a valuable resource for large-scale metagenomic analysis.

## Step 1: Preprocessing
The initial input for the function is the input sequence that should be compared against multiple proteins and all the proteins that the user wants to compare the sequence against. Since the input sequence is given in a nucleotide sequence but the proteins are in amino acid sequences, the input sequence is converted to an amino acid sequence.

## Step X: Dynamic Programming
The amino acid sequences are used to match the input sequence to the proteins. The BLOSUM62 was used as the values in the dyanmic programming approach. 