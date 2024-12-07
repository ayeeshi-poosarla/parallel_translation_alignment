# ALIGNer: Advanced Large-scale Identification of Genomic Entries

The goal of ALIGNer is to develop an efficient tool for translating DNA sequencing reads into proteins and aligning them across all six possible reading frames. Using dynamic programming and a custom alignment algorithm, the tool will compare DNA sequences against known protein databases. Designed with scalability in mind, GENE-AIM leverages parallel processing and distributed computing to handle large metagenomic datasets, providing a high-confidence filtering system for rapid identification of relevant protein matches. This tool aims to enhance functional annotations and comparative genomics in complex microbial communities, making it a valuable resource for large-scale metagenomic analysis.

## Step 1: Preprocessing
The initial input for the function is the input sequence that should be compared against multiple proteins and all the proteins that the user wants to compare the sequence against. Since the input sequence is given in a nucleotide sequence but the proteins are in amino acid sequences, the input sequence is converted to an amino acid sequence.

## Step 2: Filtering  
In the filtering step, the input reference sequence is compared to a set of protein sequences using k-mer indexing. A hashmap is built where each k-mer from the reference sequence is mapped to its matching positions in each protein sequence. The proteins are then ranked based on the number of k-mer matches with the reference, and only the top-ranked proteins, those with the highest number of matches, are passed forward for dynamic programming alignment.

## Step 4: Parallel Processing 
This step combines the previous code and sends the proteins that have the 3 highest matched through kmer indexing to the dynamic programming algorithm, which then retrieves the protein that is the greatest match for that specific gene. These files are located in official_program. 

official_program contains tests called "official_tests_keep" that were used to test the functinality of the method. The final working method is located in final_program_manual.py. Here, the code is set-up, with one of our test cases (multiple_match) and retrieves the best protein for each gene in the terminal. It also calculates the parallel processing, sequential processing, and f1 score. 

## Steps to Reproduce the Results

### 1. Clone the GitHub Repository
1. Open a **Command Prompt** or **Terminal** in VS Code and run the following command:
2. Navigate to the directory where you want to store the project (for example, cd C:\Users\reach\OneDrive\Documents\2024-25\FALL\Computational Genomics)
3. Clone the repository using the git clone command: git clone https://github.com/ayeeshi-poosarla/parallel_translation_alignment.git
4. After cloning, navigate to the newly created project directory: cd <parallel_translation_alignment>

### 2. Set Up Python and Required Dependencies
1. Install Python (if not installed), download and install Python from python.org. During installation, check the option to Add Python to PATH.
2. Set Up a Virtual Environment (Optional but Recommended). Install virtualenv if you don't have it: pip install virtualenv
3. Create a virtual environment in your project directory: python -m venv venv
4. Activate the virtual environment: .\venv\Scripts\activate

### 3. Set Up Apache Spark
1. Download and Set Up Apache Spark from the official website https://spark.apache.org/downloads.html
2. Extract the Spark folder to the directory (for example, C:\Users\reach\OneDrive\Documents\2024-25\FALL\Computational Genomics\parallel_translation_alignment\spark)
3. Set Environment Variables
    a. Set the SPARK_HOME and JAVA_HOME environment variables by going to ... 
    b. Add C:\Users\reach\OneDrive\Documents\2024-25\FALL\Computational Genomics\parallel_translation_alignment\spark\bin to the PATH.
4. Test Spark Installation: Test that Spark is correctly set up by running: spark-shell
5. Open VS Code and open your project directory: code .
6. Install VS Code Extensions such as Python if you don't already have it
### 4. Run the Code
In VS Code’s terminal, run the Python script:
python <final_program_manual.py>
1. Running Spark Jobs
To run Spark with Python, set the PYSPARK_PYTHON environment variable:
set PYSPARK_PYTHON=python
2. Run the Spark job using the spark-submit command:
spark-submit <final_program_manual.py>
3. View results
4. To test other test cases, final_program_automatic.py should be used. Here the location of the reference (fastq file) and proteins (fasta file) should be submitted in the locations reference_genomes_file =  and protein_file =. Some example fasta and fastq files to run are located officail_tests_keep.py. 

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
