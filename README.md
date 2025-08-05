## heddc_acc
`heddc_acc` is an accelerated version of hEDDC (Kawahara & Morishita, 2025) which calculates the edit distance with duplications and contractions (EDDC) between tandem repeats (TRs), including complex TRs composed of multiple units.
In this program, precomputations of distances between single units and unit sequences are omitted for long unit sequences, thus speeding up the entire program. Currently, the default parameters are 1.0 for mutations and indels (per base), and 0.5 * unit length for duplications/contractions.

## Usage
### Requirements and installation
After installing GCC with C++20 support, run the following commands in the directory.
```bash
git clone https://github.com/sefilata/heddc_acc
cd heddc_acc
make
```

### Running heddc_acc 
Run heddc_acc with the following command.
```bash
./heddc_acc -f reads.fasta -u units.fasta -s scores.tsv -v variations.tsv -t time.txt -e encodings.txt
```

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
Outputs the edit operations (mutations, indels, duplications/contractions) which minimize the cost for each pair of input TRs. Duplications/contractions are shown in (unit length, number of operations) format.
```variations.tsv
{mut:0, indel:0, dup:}	{mut:0, indel:0, dup:(3, 1)}	{mut:1, indel:0, dup:}
{mut:0, indel:0, dup:(3, 1)}	{mut:0, indel:0, dup:}	{mut:0, indel:0, dup:(3, 1)}
{mut:1, indel:0, dup:}	{mut:0, indel:0, dup:(3, 1)}	{mut:0, indel:0, dup:}
```
- `-t` : Output execution time file (Optional)  
The first line shows the total execution time of the entire program (excluding string decomposer).
```time.txt
70 msec
62 msec (c1,c2,valid_rules)
5 msec (f_scores)
2 msec (main dp)
```
- `-e` : output encodings file (Optional)  
Outputs the result of string decomposer in FASTA-like format. The correspondence of unit IDs and the their sequences is shown in the first line, in (number, unit sequence) format.
```encodings.txt
# units: (0, ACC) (1, AGC) 
> seq1
0 0 0 0 0 0 1 
> seq2
0 0 0 0 0 1 
> seq3
0 0 0 0 0 1 1 
```

## References
- Tamar Pinhas, Shay Zakov, Dekel Tsur, Michal Ziv-Ukelson, Efficient edit distance with duplications and contractions, *Algorithms Mol Biol* 8, 27, October 2013, [https://doi.org/10.1186/1748-7188-8-27](https://doi.org/10.1186/1748-7188-8-27)
- Tatiana Dvorkina, Andrey V Bzikadze, Pavel A Pevzner, The string decomposition problem and its applications to centromere analysis and assembly, *Bioinformatics* 36, Supplement_1, July 2020, [https://doi.org/10.1093/bioinformatics/btaa454](https://doi.org/10.1093/bioinformatics/btaa454)
- Riki Kawahara, Shinichi Morishita, Approximating edit distances between complex tandem repeats efficiently, *Bioinformatics* 41, 4, April 2025, [https://doi.org/10.1093/bioinformatics/btaf155](https://doi.org/10.1093/bioinformatics/btaf155)
- Bansho Masutani, Riki Kawahara, Shinichi Morishita, Decomposing mosaic tandem repeats accurately from long reads, *Bioinformatics* 39, 4 April 2023, [https://doi.org/10.1093/bioinformatics/btad185](https://doi.org/10.1093/bioinformatics/btad185)
