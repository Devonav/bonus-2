README
PBWT-Based Haplotype Matching
This project implements a Positional Burrows-Wheeler Transform (PBWT) to identify and analyze maximal haplotype matches within a dataset of binary sequences. The PBWT efficiently computes prefix and divergence arrays, which are then used to identify matches of specified lengths across haplotype sequences.

Features
PBWT Construction: Efficiently constructs positional prefix arrays (a_k) and divergence arrays (d_k).
Haplotype Transformation: Converts site-centric data to haplotype-centric format.
Maximal Match Detection: Identifies all maximal matches longer than a specified length (L).
Sample Display: Displays a sample of the identified matches.
Validation: Ensures the dataset adheres to expected format and constraints.
Optional Features:
Visualize match distribution by position.
Save all identified matches to a file for further analysis.
Dataset Format
Input: A list of binary strings representing site-centric haplotype data.
Each string corresponds to one site across all haplotypes.
Example: ["100", "110", "011"] for 3 sites and 3 haplotypes.
How to Use
1. Dependencies
Ensure Python 3.x is installed. Optionally, matplotlib can be used for match distribution visualization.

Install additional libraries if needed:

bash
Copy code
pip install matplotlib
2. Execute the Script
Run the script with the provided dataset_site_centric or modify it with your data:

bash
Copy code
python pbwt_matching.py
3. Parameters
Modify the parameters directly in the script:

L: Minimum match length (default: 4).
M: Number of haplotypes (default: 30).
N: Number of sites (default: 100).
Output
Console Output:

Total matches identified.
A sample of maximal matches with sequence indices and positions.
Optional File Output:

Save all matches to all_matches.txt.
Optional Visualization:

Plot distribution of matches by position.
Code Overview
Functions
build_pbwt(sequences, N): Constructs positional prefix arrays and divergence arrays for the dataset.

report_long_matches(a, d, sequences, N, L): Identifies all maximal matches longer than L.

display_sample_matches(matches, sample_size): Displays a sample of matches in the console.

transform_site_to_haplotype(site_centric_data, M, N): Converts site-centric data into haplotype-centric format.

Main Script
Transforms site-centric data to haplotype-centric format.
Builds PBWT arrays.
Identifies and displays maximal matches longer than the specified length.
Example Output
plaintext
Copy code
Total Maximal Haplotype Matches longer than 4: 50

Sample Matches:
Sequence 3 and Sequence 7 match from position 5 to 12
Sequence 10 and Sequence 18 match from position 2 to 8
...
Optional: Match Visualization
To visualize the match distribution:

Uncomment the visualization section in the script.
Run the script to generate a plot.
