##### CoMPoSE custom probe design

##### Purpose
Workflow for custom probe design for transgenes to be spiked into single cell probe-based detection with 10x flex. 
Requires input FASTA of target sequences. These are used to nominate 50 bp probe binding sites meeting sequence composition constraints. 
Scressn the resulting 25 bp lfet and right had arms agains the human (or mouse if modified) transcriptome for off-target hybridization. 
Selects non-overlapping sets to construct an IDT pool. Includes optional functions to restrict sequences or plot qc figures. 

#### Author
Written by Jason Swinderman
Reviewed by Aidan Winters
Revised and error-tested with Claude Opus 5

#### Date
Written November, 8th 2023\
Last revised: October, 6th 2026

#### Contents

```
probe_design.py              nomination, off-target handling, selection, oPool construction
probe_design_utils.py        input helpers 
probe_design_qc.py           QC plots
probe_design.ipynb           worked examples for domain and gapfill ORF probe sets in the paper
scripts/
  make_hs_refseq_blastdb.sh  build the human RefSeq RNA BLAST database
  probe_blast.sh             BLAST probe arms and filter by alignment length
references/                  input tables and third-party reference files
probe_nomination/            pipeline outputs
opools/                      IDT order sheets
```

#### Requirements

```bash
conda env create -f environment.yml
conda activate compose-probe-design
```

#### Reference files you need to download

1. **Chromium Human Transcriptome Probe Set v1.0.1 (GRCh38-2020-A)**, from the
   10x Genomics support site. Place it at
   `references/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv`
   — this is the default path in `gene_to_ensembl_wta_filter()` and
   `construct_probe_reference()`.
   https://www.10xgenomics.com/support/software/cell-ranger/downloads

2. **Human RefSeq RNA BLAST database**:

   ```bash
   ./scripts/make_hs_refseq_blastdb.sh -o references
   ```

### NCBI BLAST of RHS and LHS probe candidates against reference transcriptome
#### lhs_rhs_probe_binding_site_split() 
    • Constructs FASTA of the nominated probe-binding sites for the lhs and rhs sequences. 
    • This is the input for the below shell script(s). It's easiest to run this in commandline. 
#### make_hs_refseq_blastdb.sh 
    • Constructs a human refseq blast database 
#### domain_probe_blast.sh 
    • BLASTs nominated LHS and RHS probe binding sites and records those with more than 20 bp recognition of sequences in the reference transcriptome

## Filter off-target binding probes from BLAST results
### read_blast_results() 
    • imports and subsets blast off-target hits annotating lhs or rhs probes with >20 offtarget priming 
    • annotates each probe with the maximum off-target gene 
### process_lhs_rhs() 
    • concatenates rhs and lhs off-target binding infromation from blast with the probe annotation table from nominate_probes()
    
## Nominating top sets of probes for each target. 
### nominate_top_probe_set() calls:
    • calculate_priority_score() to calculate priority score for each probe pair 
    • nominate_non_overlapping_sets() to generate all unique sets of nonverlapping probe pairs
    • select_highest_priority_sets() to assign the set of probe-pairs with the greatest total priority score. 

#### Citation
Scalable probe-based single-cell transcriptional profiling for virtual cell perturbation mapping and synthetic biology phenotyping
Jason T. Swinderman, Po-Yuan Tung, Aidan Winters, Laine Goudy, Caroline M. Wilson, Lexi R. Bounds, Noam Teyssier, Ayush Agrawal, Alex Dobin, Tony Hua, Hani Goodarzi, Felix Y. Feng, Alex Marson, Dave P. Burke, Patrick D. Hsu, Yusuf H. Roohani, Silvana Konermann, Michael Kosicki, Nianzhen Li, Luke A. Gilbert
bioRxiv 2026.02.04.703058; doi: https://doi.org/10.64898/2026.02.04.703058

#### License

MIT License

Copyright (c) 2026 Jason Swinderman

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


