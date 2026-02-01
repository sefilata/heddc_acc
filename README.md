## heddc_acc
`heddc_acc` is an accelerated version of hEDDC (Kawahara & Morishita, 2025) which calculates the edit distance with duplications and contractions (EDDC) between tandem repeats (TRs), including complex TRs composed of multiple units.  
In this program, precomputations of distances between single units and unit sequences are omitted for long unit sequences, thus speeding up the entire program. Currently, the default parameters are 1.0 for mutations and indels (per base), and 0.5 * unit length for duplications/contractions.

## Usage
### Requirements and installation
After installing GCC with C++20 support, run the following commands in the directory.
```bash
git clone https://github.com/sefilata/heddc_acc.git
cd heddc_acc
make
```

### Running heddc_acc 
Run heddc_acc with the following command.
```bash
./heddc_acc -f <read_fasta> -u <unit_fasta> [-s <score_tsv>] [-v <variant_tsv>] [-t <time_txt>] [-e <encodings_txt>] [-p1 <parallel num (f-scores)>] [-p2 <parallel num (main-dp)>] [-p <parameter_pattern>] [-a <acc_pattern>]
```
An example of using heddc_acc to reconstruct large phylogenetic trees is shown in https://github.com/sefilata/analysis2_phylogenetic_trees

### Command line options and example files
- `-f` : Input FASTA file of TRs (Required)  
Each sequence must be in a single line.
```reads.fasta
> seq1
ACCACCACCACCACCACCAGC
> seq2
ACCACCACCACCACCAGC
> seq3
ACCACCACCACCACCAGCAGC
```
- `-u` : Input FASTA file of units (Required)  
Each sequence must be in a single line. You can obtain units from TRs using other software such as `uTR` (Masutani et.al., 2023).
```units.fasta
> unit1
ACC
> unit2
AGC
```
- `-s` : Output file of edit distance (Optional)  
Outputs a TSV matrix of edit distance with duplications and contractions. Rows and columns follow the order of the input TR FASTA file. You can change the decimal precision by changing `SCORE_PRECISION` variable in `hEDDC_cpp_count/main.cpp`.
If this option is omitted, the result is printed to stdout.
```scores.tsv
0	0.077151675	0.047619048
0.077151675	0	0.077151675
0.047619048	0.077151675	0
```
- `-v` : Output variations file (Optional)  
Outputs the edit operations (mutations, indels, duplications/contractions) which minimize the cost for each pair of input TRs. Duplications/contractions are shown in (unit sequence, number of dup/cons) format.
```variations.tsv
{mut:0, indel:0, dup:}	{mut:0, indel:0, dup:(ACC, 1)}	{mut:1, indel:0, dup:}
{mut:0, indel:0, dup:(ACC, 1)}	{mut:0, indel:0, dup:}	{mut:0, indel:0, dup:(AGC, 1)}
{mut:1, indel:0, dup:}	{mut:0, indel:0, dup:(AGC, 1)}	{mut:0, indel:0, dup:}
```
- `-t` : Output execution time file (Optional)  
Line 7 and 8 shows the execution time of the entire program (excluding string decomposer).
```time.txt
3 reads before compression
3 reads after compression
1 parallel (other than main dp)
1 parallel (main dp)
0 ACC_PATTERN (0:hEDDC, 1:heddc_acc)
0 PARAMETER_PATTERN (0:0.5*|unit|, 1:all 0.6)
5 msec
5085 usec
3 msec (c1,c2,valid_rules)
1 msec (f_scores)
0 msec (main dp)
3752 kB MAX (after c1,c2, valid_rules)
4072 kB MAX (after f_scores)
4072 kB MAX (after main dp)
```
- `-e` : Output encodings file (Optional)  
Outputs the result of string decomposer in FASTA-like format. The correspondence of unit IDs and their sequences is shown in the first line, in (unit ID, unit sequence) format.
```encodings.txt
# units: (0, ACC) (1, AGC) 
> seq1
0 0 0 0 0 0 1 
> seq2
0 0 0 0 0 1 
> seq3
0 0 0 0 0 1 1 
```
- `-p1` : Number of parallel threads for the first stage (Optional)  
Number of threads used to compute c1, c2, valid_rules, and f_scores.
The default value is 1.
- `-p2` : Number of parallel threads for the second stage (Optional)  
Number of threads used to compute the main DP.
This value should be smaller than p1, since the main DP requires more memory.
The default value is 1.
- `-p`  : Algorithm selection (Optional)  
Set to 0 to use the original hEDDC, and to 1 to use heddc_acc.
The default value is 1 (heddc_acc).
- `-a`  : Cost parameter selection (Optional)  
Set to 0 to use a duplication/contraction cost of 0.5 × unit length.
Set to 1 to use a fixed duplication/contraction cost of 0.6.
The default value is 1 (all costs = 0.6).
These values (e.g., 0.5 and 0.6) can be modified via the
PARAM_0_VALUE and PARAM_1_VALUE variables in heddc_acc_parallel.h.

## References
- Tamar Pinhas, Shay Zakov, Dekel Tsur, Michal Ziv-Ukelson, Efficient edit distance with duplications and contractions, *Algorithms Mol Biol* 8, 27, October 2013, [https://doi.org/10.1186/1748-7188-8-27](https://doi.org/10.1186/1748-7188-8-27)
- Tatiana Dvorkina, Andrey V Bzikadze, Pavel A Pevzner, The string decomposition problem and its applications to centromere analysis and assembly, *Bioinformatics* 36, Supplement_1, July 2020, [https://doi.org/10.1093/bioinformatics/btaa454](https://doi.org/10.1093/bioinformatics/btaa454)
- Riki Kawahara, Shinichi Morishita, Approximating edit distances between complex tandem repeats efficiently, *Bioinformatics* 41, 4, April 2025, [https://doi.org/10.1093/bioinformatics/btaf155](https://doi.org/10.1093/bioinformatics/btaf155)
- Bansho Masutani, Riki Kawahara, Shinichi Morishita, Decomposing mosaic tandem repeats accurately from long reads, *Bioinformatics* 39, 4 April 2023, [https://doi.org/10.1093/bioinformatics/btad185](https://doi.org/10.1093/bioinformatics/btad185)
