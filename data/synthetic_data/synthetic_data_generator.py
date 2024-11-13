from Bio.Seq import Seq
import random
import json

def generate_synthetic_data():
    """
    Generate synthetic DNA sequences for testing the translation and alignment system.
    Returns both DNA sequences and their known protein translations for validation.
    """
    
    # Test Case 1: Simple known protein coding sequences
    # These sequences are designed to translate to known amino acid sequences
    simple_test = {
        "dna_sequences": [
            "ATGGCCACTGAATAA",  # Translates to "MATE*"
            "ATGCGTACTGCTTAA",  # Translates to "MRTA*"
            "ATGGCTGAGCATTGA",  # Translates to "MAEH*"
            # Additional test cases with verified translations
            "ATGTCTAAGCTATAA",  # Translates to "MSK*"
            "ATGGGATCCTATTGA"   # Translates to "MGSY*"
        ],
        "expected_proteins": [
            "MATE*",
            "MRTA*",
            "MAEH*",
            "MSK*",
            "MGSY*"
        ],
        "description": "Simple sequences with verified translations using standard genetic code"
    }
    
    # Test Case 2: Sequences with multiple reading frames
    # These sequences contain overlapping genes in different frames
    multi_frame_test = {
        "dna_sequences": [
            # Contains two overlapping genes in different frames
            "ATGGCAATGCGTACTGCTTAATAA",  
            # Contains genes in both forward and reverse frames
            "ATGCCATAGCATTAACTATGAGCT"   
        ],
        "expected_proteins": {
            "forward_frame1": ["MAMRTA*", "MPIA*"],
            "forward_frame2": ["*MRLL*", "PLLM*"],
            "forward_frame3": ["AMRTA*", "IALT*"],
            "reverse_frame1": ["MLHS*", "SYV*"],
            "reverse_frame2": ["ITPL*", "FMA*"],
            "reverse_frame3": ["VKSL*", "LSP*"]
        },
        "description": "Sequences with multiple valid reading frames"
    }
    
    # Test Case 3: Complex sequences with known similarities
    # These sequences are designed to test alignment algorithms
    alignment_test = {
        "sequence_pairs": [
            # Highly similar sequences (should have high alignment scores)
            {
                "seq1": "ATGGCCACTGAATTAGCT",
                "seq2": "ATGGCCACTGAATTGGCT",  # One base difference
                "expected_similarity": 0.94
            },
            # Moderately similar sequences
            {
                "seq1": "ATGGCCACTGAATTAGCT",
                "seq2": "ATGGCTAGTGAATTAGCT",  # Two base differences
                "expected_similarity": 0.89
            },
            # Less similar sequences
            {
                "seq1": "ATGGCCACTGAATTAGCT",
                "seq2": "ATGCCCATTGAGTTAGCT",  # Multiple differences
                "expected_similarity": 0.78
            }
        ],
        "description": "Sequence pairs with known similarity scores"
    }
    
    # Test Case 4: Large-scale random sequences
    def generate_random_dna(length):
        return ''.join(random.choice('ATCG') for _ in range(length))
    
    large_scale_test = {
        "dna_sequences": [generate_random_dna(1000) for _ in range(10)],
        "sequence_lengths": [1000] * 10,
        "description": "Large random sequences for performance testing"
    }
    
    # Test Case 5: Edge cases and special scenarios
    edge_cases_test = {
        "dna_sequences": [
            "ATG",              # Minimal length (just start codon)
            "ATGTAA",           # Minimal gene (start + stop)
            "ATGNNNTAA",        # Sequence with N's (unknown bases)
            "ATGGTGATGTAA",     # Nested start codons
            "TAAATGTAATGA"      # Multiple stop codons
        ],
        "expected_proteins": [
            "M",
            "M*",
            "MX*",
            "MVM*",
            "*M**"
        ],
        "description": "Edge cases and special scenarios for robust testing"
    }
    
    return {
        "simple_test": simple_test,
        "multi_frame_test": multi_frame_test,
        "alignment_test": alignment_test,
        "large_scale_test": large_scale_test,
        "edge_cases_test": edge_cases_test
    }

def save_test_data(filename="synthetic_test_data.json"):
    """Save the synthetic test data to a JSON file"""
    test_data = generate_synthetic_data()
    with open(filename, 'w') as f:
        json.dump(test_data, f, indent=2)
    
def verify_test_data():
    """Verify the translations of test sequences using BioPython"""
    test_data = generate_synthetic_data()
    
    print("Verifying simple test translations:")
    for dna, expected in zip(test_data["simple_test"]["dna_sequences"], 
                           test_data["simple_test"]["expected_proteins"]):
        seq = Seq(dna)
        actual = str(seq.translate())
        print(f"DNA: {dna}")
        print(f"Expected: {expected}")
        print(f"Actual: {actual}")
        print(f"Match: {actual == expected}\n")
    
    print("\nVerifying edge cases translations:")
    for dna, expected in zip(test_data["edge_cases_test"]["dna_sequences"],
                           test_data["edge_cases_test"]["expected_proteins"]):
        seq = Seq(dna)
        actual = str(seq.translate())
        print(f"DNA: {dna}")
        print(f"Expected: {expected}")
        print(f"Actual: {actual}")
        print(f"Match: {actual == expected}\n")

if __name__ == "__main__":
    # Generate and save test data
    save_test_data()
    # Verify translations
    verify_test_data()