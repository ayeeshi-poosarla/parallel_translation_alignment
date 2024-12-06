# ALIGNer: Advanced Large-scale Identification of Genomic ENtries

The goal of ALIGNer is to develop an efficient tool for translating DNA sequencing reads into proteins and aligning them across all six possible reading frames. Using dynamic programming and a custom alignment algorithm, the tool will compare DNA sequences against known protein databases. Designed with scalability in mind, GENE-AIM leverages parallel processing and distributed computing to handle large metagenomic datasets, providing a high-confidence filtering system for rapid identification of relevant protein matches. This tool aims to enhance functional annotations and comparative genomics in complex microbial communities, making it a valuable resource for large-scale metagenomic analysis.

## Step 1: Preprocessing
The initial input for the function is the input sequence that should be compared against multiple proteins and all the proteins that the user wants to compare the sequence against. Since the input sequence is given in a nucleotide sequence but the proteins are in amino acid sequences, the input sequence is converted to an amino acid sequence.

## Step X: Dynamic Programming
The amino acid sequences are used to match the input sequence to the proteins. The BLOSUM62 was used as the values in the dyanmic programming approach. 

## Step X: Parallel Processing 
This step combines the previous code and sends the proteins who have the 3 highest matched through kmer indexing to the dynamic programming algorithm, which then retrieves the protein that is the greatest match for that specific gene. These files are located in official_program. 

official_program contains tests called "official_tests_keep" that were used to test the functinality of the method. The final working method is located in final_program_manual.py. Here, the code is set-up, with one of our test cases (multiple_match) and retrieves the best protein for each gene in the terminal. It also calculates the parallel processing, sequential processing, and f1 score. 

To run this code, simply run final_program_manual.py. To test other test cases, final_program_automatic should be used. Here the location of the reference (fastq file) and proteins (fasta file) should be submitted in the locations reference_genomes_file =  and protein_file =. Some example fasta and fastq files to run are located officail_tests_keep. 

# Test Case Generation

Test cases are generated using `generate_test_cases.py`. The script creates synthetic DNA and protein sequences with controlled similarity levels:

```
test_cases/
├── exact_match/        # 100% similarity
├── high_similarity/    # 90% similarity  
├── moderate_similarity/# 70% similarity
└── low_similarity/     # 40% similarity
```

Each directory contains sequence files (`sequence.fasta`, `protein.fasta`) and `metadata.json` with expected scores. BioPython's `Seq` class generates valid sequences and reading frames.