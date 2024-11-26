import sys

def build_pbwt(sequences, N):
    """
    Constructs the positional prefix arrays (a_k) and divergence arrays (d_k) for PBWT.

    Parameters:
    - sequences: List of haplotype sequences (strings of '0's and '1's).
    - N: Number of variable sites in each sequence.

    Returns:
    - a: List of positional prefix arrays.
    - d: List of divergence arrays.
    """
    M = len(sequences)
    a = [[] for _ in range(N+1)]
    d = [[] for _ in range(N+1)]

    # Initialize a0 as the list of all sequence indices
    a[0] = list(range(M))
    d[0] = [0] * M

    for k in range(N):
        a_k = a[k]
        d_k = d[k]
        a_k_plus1 = []
        d_k_plus1 = []

        # Split sequences based on current bit at position k
        a_group = []
        d_group = []
        b_group = []
        e_group = []

        for idx in a_k:
            if sequences[idx][k] == '0':
                a_group.append(idx)
                d_group.append(d_k[idx])
            else:
                b_group.append(idx)
                e_group.append(d_k[idx])

        # Concatenate a and b groups to form a_{k+1}
        a_k_plus1 = a_group + b_group
        d_k_plus1 = d_group + e_group

        a[k+1] = a_k_plus1
        d[k+1] = d_k_plus1

    return a, d

def report_long_matches(a, d, sequences, N, L):
    """
    Identifies all maximal haplotype matches longer than a specified length L.

    Parameters:
    - a: List of positional prefix arrays.
    - d: List of divergence arrays.
    - sequences: List of haplotype sequences.
    - N: Number of variable sites.
    - L: Minimum match length.

    Returns:
    - matches: List of dictionaries containing match information.
    """
    M = len(sequences)
    matches = []

    for k in range(L, N+1):
        a_k = a[k]
        d_k = d[k]
        for i in range(1, M):
            # Check if the current and previous sequences diverged at or before k - L
            if d_k[i] <= k - L:
                seq1 = a_k[i-1]
                seq2 = a_k[i]
                start = d_k[i]
                end = k
                matches.append({
                    'seq1': seq1,
                    'seq2': seq2,
                    'start': start,
                    'end': end
                })

    return matches

def display_sample_matches(matches, sample_size=10):
    """
    Displays a sample of the identified matches.

    Parameters:
    - matches: List of dictionaries containing match information.
    - sample_size: Number of sample matches to display.
    """
    print(f"Total Maximal Haplotype Matches longer than 4: {len(matches)}\n")
    print("Sample Matches:")
    for idx in range(min(sample_size, len(matches))):
        match = matches[idx]
        seq1 = match['seq1']
        seq2 = match['seq2']
        start = match['start']
        end = match['end']
        print(f"Sequence {seq1} and Sequence {seq2} match from position {start} to {end}")

def main():
    """
    Main function to execute the PBWT-based haplotype matching.
    """
    # Provided dataset: list of 100 binary strings, each with 30 bits
    dataset_site_centric = [
        "100000100010110011110001001100",
        "000101001010111101100101101001",
        "111100001000010001110011011000",
        "111011110111010100001011010101",
        "011010001101011001000010111010",
        "001010011000101110011001101111",
        "101010000100000010111001101101",
        "100110110101000101111110100111",
        "011101000011100001110100011011",
        "100111011111001000010111011100",
        "001010001011100110100110111101",
        "110111011101111011110101011110",
        "011100110101001100000000101110",
        "000010101011001010010000111111",
        "000101001010111101100101101001",
        "001001010110001001100001100011",
        "001000011101010000110011000110",
        "100000101100011110011101001110",
        "010110010101001100010000001111",
        "010010000011101001110000000001",
        "110100101010011010100001011101",
        "100011101001111101100100011011",
        "110011110000111111001100001010",
        "010100011011001100011001110100",
        "111111100001101111000000110101",
        "001101110010000000100111111101",
        "001010011011010011010010101011",
        "011010001101011001000010111010",
        "101101100001101010010111100110",
        "010101111101111010110111010111",
        "101010010010011101111100011100",
        "010000010000000110110110000100",
        "011100000010111111011101100011",
        "100100010011101111101111100001",
        "001001101010100100110000010111",
        "001011011110100001010100010000",
        "110110110110111100100001110101",
        "001110110110100111010111110001",
        "011010001010100000110011011110",
        "100011101001111101100100011011",
        "011100111101110110110100111101",
        "111110010100001101010101010000",
        "100110010011011111011010010010",
        "110000011000101101001011000110",
        "101100001111011111101011011010",
        "111111000111111111110100101101",
        "100100001001000011010101000010",
        "100001010101111000010000011000",
        "101010011010011000101001101110",
        "111001111100000011111011110010",
        "000110000100011001101010001000",
        "010011001111101110101011111100",
        "111110110111110101110110001011",
        "011001011110100010100011001101",
        "001000001100000111001000111000",
        "001001010110001001100001100011",
        "101001000111001010100001111101",
        "101000101001101110000011110101",
        "101000111100010111000111110111",
        "011010101111110110011110000010",
        "001011011111110000111001001101",
        "110110001011100101110100000011",
        "011100101001111111110011011101",
        "100011001111101100101001100100",
        "000100101010111000011010001001",
        "100010101000010100100000100100",
        "110010101001101011110011001011",
        "010111100110111110000101011101",
        "000100101000111011000000111001",
        "001011011111110000111001001101",
        "001101100111010010100110111111",
        "110010011110001101110010000001",
        "000011111000011110011010000100",
        "100101111111101000110100011101",
        "011011000010010100100110110011",
        "100000110101111111100101010011",
        "111110110111110101110110001011",
        "000100011101101111101110010010",
        "111110111011110100011001101100",
        "110111100110001101010100100110",
        "011011010111010101111010010101",
        "100100010010000001010010011001",
        "000010111011100110100101001100",
        "010100010000010100110000010010",
        "001111000001010110000001100011",
        "110100100101101000011010100100",
        "000100011101000101010001101011",
        "100000011011100101001101100110",
        "111001010001100100011010100111",
        "011100100100010000110001111011",
        "011000100110100001101001101011",
        "110010000001011001011001111010",
        "001110001100000000101111010111",
        "000001111011101000001101001101",
        "111000010001110110011001010011",
        "111101101011010111100110111010",
        "010011100011010010000111010111",
        "101000111100010111000111110111",
        "011010101111110110011110000010",
        "010110111000010110111111010100"
    ]

    

    # Parameters
    L = 4  # Minimum match length

    # Step 1: Transform site-centric data to haplotype-centric data
    def transform_site_to_haplotype(site_centric_data, M=30, N=100):
        """
        Transforms site-centric data (list of site strings) to haplotype-centric data.

        Parameters:
        - site_centric_data: List of strings, each representing a site across M haplotypes.
        - M: Number of haplotypes (default 30).
        - N: Number of sites (default 100).

        Returns:
        - haplotypes: List of M haplotype strings, each of length N.
        """
        haplotypes = ['' for _ in range(M)]
        for site in site_centric_data:
            for hap_idx in range(M):
                haplotypes[hap_idx] += site[hap_idx]
        return haplotypes

    haplotypes = transform_site_to_haplotype(dataset_site_centric, M=30, N=100)

    # Validate dataset
    for idx, seq in enumerate(haplotypes):
        if len(seq) != 100:
            print(f"Error: Haplotype {idx} length is {len(seq)}, expected 100.")
            sys.exit(1)
        if any(c not in '01' for c in seq):
            print(f"Error: Haplotype {idx} contains invalid characters.")
            sys.exit(1)

    # Step 2: Build PBWT arrays
    a, d = build_pbwt(haplotypes, N=100)

    # Step 3: Identify maximal matches longer than L=4
    matches = report_long_matches(a, d, haplotypes, N=100, L=4)

    # Step 4: Display results
    display_sample_matches(matches, sample_size=10)

    # Optional: Display additional information
    # For example, visualize match distribution
    """
    import matplotlib.pyplot as plt

    def count_matches_per_position(matches, N):
        count = [0] * (N+1)
        for match in matches:
            end = match['end']
            count[end] +=1
        return count

    def plot_match_distribution(match_counts, L, N):
        positions = list(range(L, N+1))
        counts = [match_counts[k] for k in positions]

        plt.figure(figsize=(12,6))
        plt.plot(positions, counts, marker='o')
        plt.title('Distribution of Maximal Haplotype Matches by Position')
        plt.xlabel('Position')
        plt.ylabel('Number of Matches')
        plt.grid(True)
        plt.show()

    # Count matches per position
    match_counts = count_matches_per_position(matches, 100)

    # Plot the distribution
    plot_match_distribution(match_counts, L, 100)
    """

    # Optional: Save all matches to a file
    """
    with open("all_matches.txt", "w") as f:
        for match in matches:
            seq1 = match['seq1']
            seq2 = match['seq2']
            start = match['start']
            end = match['end']
            hap1 = haplotypes[seq1][start:end]
            hap2 = haplotypes[seq2][start:end]
            f.write(f"Sequence {seq1} and Sequence {seq2} match from position {start} to {end}:\n")
            f.write(f"  Haplotype {seq1}: {hap1}\n")
            f.write(f"  Haplotype {seq2}: {hap2}\n\n")
    print("\nAll matches have been saved to 'all_matches.txt'.")
    """

if __name__ == "__main__":
    main()
