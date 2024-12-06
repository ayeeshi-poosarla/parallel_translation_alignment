from dataclasses import dataclass
from typing import List, Dict, Tuple
from Bio.Seq import Seq
import json
import random
import os

@dataclass
class TestCase:
    name: str
    description: str
    dna_sequence: str
    protein_match: str
    protein_sequence: str
    expected_translations: List[str]
    edge_case_type: str = ""
    reading_frames: Dict[str, str] = None
    expected_alignment_score: float = 0.0

def get_reading_frames(sequence: str) -> Dict[str, str]:
    """Get all six possible reading frames for a sequence."""
    seq = Seq(sequence)
    frames = {}
    
    # Forward frames
    for i in range(3):
        frame_seq = seq[i:len(seq) - ((len(seq)-i) % 3)]
        if len(frame_seq) >= 3:
            frames[f"frame{i+1}"] = str(frame_seq.translate())
    
    # Reverse complement frames
    rev_seq = seq.reverse_complement()
    for i in range(3):
        frame_seq = rev_seq[i:len(rev_seq) - ((len(rev_seq)-i) % 3)]
        if len(frame_seq) >= 3:
            frames[f"reverse_frame{i+1}"] = str(frame_seq.translate())
    
    return frames

def generate_protein_sequence(length: int) -> str:
    """Generate a protein sequence of specified length."""
    amino_acids = 'ARNDCEQGHILKMFPSTWYV'
    return ''.join(random.choice(amino_acids) for _ in range(length))

def create_matching_sequences(protein_length: int, similarity: float) -> Tuple[str, str]:
    """Create a protein sequence and a similar variant."""
    original = generate_protein_sequence(protein_length)
    variant = list(original)
    
    # Number of positions to mutate
    mutations = int(protein_length * (1 - similarity))
    positions = random.sample(range(protein_length), mutations)
    
    for pos in positions:
        amino_acids = 'ARNDCEQGHILKMFPSTWYV'.replace(original[pos], '')
        variant[pos] = random.choice(amino_acids)
    
    return original, ''.join(variant)

def create_test_case(case_type: str) -> TestCase:
    """Create a test case based on specified type."""
    
    if case_type == "exact_match":
        # Create perfectly matching sequences
        protein_seq = generate_protein_sequence(20)
        dna_seq = "ATG" + ''.join(random.choice(['GCT', 'GCC']) for _ in range(19)) + "TAA"
        alignment_score = 1.0
        
    elif case_type == "high_similarity":
        # Create highly similar sequences (90% similarity)
        protein_seq, variant = create_matching_sequences(20, 0.9)
        dna_seq = "ATG" + ''.join(random.choice(['GCT', 'GCC']) for _ in range(19)) + "TAA"
        alignment_score = 0.9
        
    elif case_type == "moderate_similarity":
        # Create moderately similar sequences (70% similarity)
        protein_seq, variant = create_matching_sequences(20, 0.7)
        dna_seq = "ATG" + ''.join(random.choice(['GCT', 'GCC']) for _ in range(19)) + "TAA"
        alignment_score = 0.7
        
    elif case_type == "low_similarity":
        # Create sequences with low similarity (40% similarity)
        protein_seq, variant = create_matching_sequences(20, 0.4)
        dna_seq = "ATG" + ''.join(random.choice(['GCT', 'GCC']) for _ in range(19)) + "TAA"
        alignment_score = 0.4
        
    else:
        raise ValueError(f"Unknown case type: {case_type}")

    # Get all reading frames
    frames = get_reading_frames(dna_seq)
    
    return TestCase(
        name=case_type,
        description=f"Test case for {case_type}",
        dna_sequence=dna_seq,
        protein_match=f"test_protein_{case_type}",
        protein_sequence=protein_seq,
        expected_translations=list(frames.values()),
        reading_frames=frames,
        expected_alignment_score=alignment_score
    )

def create_directory_structure():
    """Create directory structure for test cases."""
    base_dir = "test_cases"
    os.makedirs(base_dir, exist_ok=True)
    return base_dir

def save_test_case(test_case: TestCase, base_dir: str):
    """Save a single test case to its own directory."""
    # Create directory for this test case
    test_dir = os.path.join(base_dir, test_case.name)
    os.makedirs(test_dir, exist_ok=True)
    
    # Save DNA sequence
    with open(os.path.join(test_dir, "sequence.fasta"), 'w') as f:
        f.write(f">{test_case.name} {test_case.description}\n")
        f.write(f"{test_case.dna_sequence}\n")
    
    # Save protein sequence
    with open(os.path.join(test_dir, "protein.fasta"), 'w') as f:
        f.write(f">{test_case.protein_match}\n")
        f.write(f"{test_case.protein_sequence}\n")
    
    # Save metadata
    metadata = {
        "description": test_case.description,
        "protein_match": test_case.protein_match,
        "expected_alignment_score": test_case.expected_alignment_score,
        "reading_frames": test_case.reading_frames
    }
    
    with open(os.path.join(test_dir, "metadata.json"), 'w') as f:
        json.dump(metadata, f, indent=2)

def main():
    # Create directory structure
    base_dir = create_directory_structure()
    
    # Generate and save test cases
    case_types = ["exact_match", "high_similarity", "moderate_similarity", "low_similarity"]
    
    print("Generating test cases...")
    for case_type in case_types:
        test_case = create_test_case(case_type)
        save_test_case(test_case, base_dir)
        print(f"Created test case: {case_type}")
        
    print("\nTest case structure:")
    print("test_cases/")
    for case_type in case_types:
        print(f"├── {case_type}/")
        print("│   ├── sequence.fasta")
        print("│   ├── protein.fasta")
        print("│   └── metadata.json")

if __name__ == "__main__":
    main()