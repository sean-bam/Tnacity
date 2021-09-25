This repository hosts the code to identify boundaries of Tn7-like tranpsosons, as described in [Benler et al. 2021](https://doi.org/10.1101/2021.08.23.457393). The source data can be obtained on the [NCBI FTP](https://ftp.ncbi.nih.gov/pub/yutinn/benler_2021/Tn7/source_data/).

# Installation

1. (Optional) It's always a good idea to work within virtual environment so there aren't any dependency clashes:
```
python3 -m venv environments/testenv
source activate environments/testenv/bin/activate
```

2. Clone this directory and get an executable version of the code, with python dependencies
```
git clone https://github.com/sean-bam/Tnacity 
cd Tnacity
chmod +x bin/tnacity.py
pip install pandas biopython
```

The code was tested on pandas `pandas == 1.2.3` and `biopython == 1.78`

3. Install the external dependencies
 - BLAST (tested on v2.12.0)
 - [MEME](https://meme-suite.org/meme/doc/install.html?man_type=web) (tested  on v5.3.0)

Make sure that both `meme/bin/` and `/meme/libexec/meme-5.3.0` are on path. Specifically, you should be able to call the programs `meme` and `fasta-get-markov` from the command line. If not, see the installation instructions of MEME, which tell you to add both of these folders to your path: 
>export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.3.0:$PATH

4. Test the installation from the Tnacity main directory;
```
cd Tnacity/
python tests/UnitTest.py 
```

If all went well, there should be a folder named `tests/Tn7/` and no error messages.

# Usage

Tnacity requires the following parameters to be set
1. `-i` : a fasta file or genbank file. 
2. `-f` : if fasta file is provided, Tnacity also requires a five-column, comma-separated table defining the open reading frames. The columns are `accession,protein_id,start,stop,feature`. An example can be found in the folder `/tests/41359_0_94.ft`.
3. `-c` : A nucleotide coordinate. This will serve as the initial starting point "anchor" to begin the search.
4. `-a` : An accession that matches the header of the fasta file or genbank file exactly. 
5. `-o` : A directory to store the output, which will be created if not present already. Each input should have its own output directory.

Tnacity can also be run with the following optional parameters:
 - `-k` : Keep all intermediate files produced during the search.
 - `-t` : Identify motifs 'de-novo' (default), or 'search' the genome using known motifs 
 - `-m`: If searching a genome using input motifs, specify the /path/to/motifs file here. Currently only accepting motifs in MEME XML format (all the motifs in the file must have a unique name)
 - `-dust` : Value for DUST low complexity filtering of input sequences (30)
 - `-bf_minlength` : Using BLASTn to finding inverted repeats, this is the minimum length of the alignment (15 bp).
 - `-bf_maxlength` : Using BLASTn to finding inverted repeats, this is the maximum length of the alignment (30 bp).
  - `-bf_maxgaps` : Using BLASTn to finding inverted repeats, this is the maximum number of gaps allowed (2).
  - `-bf_maxmm` : Using BLASTn to finding inverted repeats, this is the maximum number of mismatches allowed (6).
  - `-bf_leash` : Using BLASTn to finding inverted repeats, this is the maximum distance one end can be from the provided coordinate
  - `-m_leash` : Same as `bf_leash`, but for motifs detected by MEME/Fimo, i.e., the maximum distance one end detected can be from the provided coordinate
  - `-m_FDR` : Limits hits to those with the false-discovery rates below X, as reported by Fimo. See the Fimo [docs](https://meme-suite.org/meme/doc/fimo.html?man_type=web) for more about the the -q_value. (0.1)
  - `-m_cFDR` : The product of all motif instances q-values per end sequence must be lower than X (0.01)
  - `-m_space` : Motifs that are greater than X bp from the next closest motif are discarded (75)
  - `-m_foot` : The total footprint of all motifs on a single end sequence must be less than X (120)

## Examples 

Run Tnacity using all the default parameters. This will first run Blast to detect inverted repeats (IRs) in intergenic sequences around `103000` (the start coordinate of the Tn7 TnsB transposase), then identify motifs with MEME on all pairs of sequences with IRs. 
>bin/tnacity.py -i tests/NZ_CP067307.1.gb -infmt genbank -o tmp -c 5042301 -a NZ_CP067307.1

The output directory contains two files: 
1. `predicted_ends.csv` : a table of the predicted coordinates of the transposase binding sites, and 
2. `tns_motifs.xml` : a [meme-formatted](meme-suite.org/meme/doc/meme-format.html?man_type=web) XML file of the detected motifs. 

The motif XML file can be used to search other genomes for the presence of this motif. For example, the ends of Tn7 are not detected by Tnacity when it is located on the plasmid. This is likely due to the short length of the plasmid compared to a whole genome, which changes the markov model of the input DNA sequence. But, by providing the previously produced XML file via `-m`, one can correctly identify the boundaries using the `-search` option:
>bin/tnacity.py -i tests/NC_002525.gb -infmt genbank -o tmp2/ -c 18254 -a NC_002525.1 -t search -m tmp/tns_motifs.xml 

The [NCBI FTP](https://ftp.ncbi.nih.gov/pub/yutinn/benler_2021/Tn7/source_data/) has all the motifs reported in the paper (n ~600). You can use all of these to search a genome of interest (but it may take awhile to process all of the motifs).

# Citation
Please consider citing the paper if you find this code useful:
 - [Benler et al. 2021](https://doi.org/10.1101/2021.08.23.457393)
