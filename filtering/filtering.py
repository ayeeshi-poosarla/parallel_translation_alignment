# filtering

def build_kmer_index(reference, protein, k):
    """
    Build a hashmap that maps each k-mer in the reference sequence
    to the positions where it matches in the protein sequence.
    """
    kmer_map = {}
    for i in range(len(reference) - k + 1):
        kmer = reference[i:i + k]
        match_positions = [
            j for j in range(len(protein) - k + 1) if protein[j:j + k] == kmer
        ]
        if match_positions:
            kmer_map[kmer] = match_positions
    return kmer_map


def get_top_kmer_match(kmer_map):
    """
    Find the k-mer with the highest number of matches in the k-mer map.
    """
    if not kmer_map:
        return None, 0
    top_kmer = max(kmer_map, key=lambda k: len(kmer_map[k]))
    return top_kmer, len(kmer_map[top_kmer])


def rank_proteins_by_matches(reference, proteins, k):
    """
    Rank proteins by the maximum number of k-mer matches they have with the reference.
    """
    protein_match_data = []

    for protein in proteins:
        kmer_map = build_kmer_index(reference, protein, k)
        _, max_matches = get_top_kmer_match(kmer_map)
        protein_match_data.append((protein, max_matches))

    # Sort proteins by the number of matches in descending order
    protein_match_data.sort(key=lambda x: x[1], reverse=True)

    # Return only the protein names in ranked order
    return [protein for protein, _ in protein_match_data]


def main():
    # Example input data
    sequences = [
        "ACGTGCTAGCTA",  # Reference sequence
        "TACGTGCTAGGCTAAGCTA",  # Protein 1
        "GTGCTAGCTACGTAGCTA",   
        "ACGTGCTAACGTAGCTAA"    
    ]

    # First sequence is the reference
    reference = sequences[0]
    proteins = sequences[1:]
    k = 3  # Length of k-mers to consider

    # Get the ranked protein names
    ranked_proteins = rank_proteins_by_matches(reference, proteins, k)

    # Print the ranked proteins
    print("Ranked Proteins:")
    print(ranked_proteins)


if __name__ == "__main__":
    main()
