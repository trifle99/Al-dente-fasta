```
·················································································
:     _     _         ____                _          _____            _         :
:    / \   | |       |  _ \   ___  _ __  | |_  ___  |  ___|__ _  ___ | |_  __ _ :
:   / _ \  | | _____ | | | | / _ \| '_ \ | __|/ _ \ | |_  / _` |/ __|| __|/ _` |:
:  / ___ \ | ||_____|| |_| ||  __/| | | || |_|  __/ |  _|| (_| |\__ \| |_| (_| |:
: /_/   \_\|_|       |____/  \___||_| |_| \__|\___| |_|   \__,_||___/ \__|\__,_|:
:                                                                               :
·················································································
```

<!-- About the project -->
## About the project

This is a python project I started to help practice my python coding and improve my understanding of sequencing data.
I created a python version of FastQC where this package 'Al-Dente Fasta' parses fastq files and allows you to perform quality control tests on sequencing datasets.

Currently, this tool only works on fastq files where all reads are the same length.
The analysis it provides are:
* Basic summary (total base counts, total read counts, total gc%)
* Per base position quality score across each read box plot
* Per base position average quality scores across each read plot
* Per base read content plot 
* Per read GC% content plot
* Per tile read avg quality scores plot (for illumina sequences only as they contain extra information on their sequence id's)

<!-- Getting started -->
## Getting started

To use this package, you can first install it by using: 

```pip install al-dente-fasta==0.1.0```

### Running quality control analysis:
1. [Per base position quality score across each read](/images/Figure_1.png)
   
   (Box plot analysis > shows the distribution of quality scores of each base position across each read in a box plot)
   ```
   import al-dente_fasta.preprocess as afp
   afp.qual_out('path/to/file_in/fastq.gz', 'path/to/file_out/qual_out')
   afp.qual_convert('path/to/file_out/qual_out', 'path/to/file_out/qual_convert')
   afp.quals_boxplot('path/to/file_out/qual_convert')
   ```

2. [Per read average quality score count](/images/Figure_2.png)
   
   (Takes the average quality score per read and shows a plot to visualize the counts distribution)
   ```
   import al-dente_fasta.counts_graphs as afcg
   afcg.avg_quals('path/to/file_out/qual_convert', 'path/to/file_out/avg_quals'
   #kde_or_histo takes in values of 0/1/2, enter 0 to view a histogram plot, 1 for a kernal density estimate plot and 2 for a combination of both.
   afcg.counts_graph('path/to/file_out/avg_quals', histo_or_kde, 'qual')
   
   ```

3. [Per read GC% content count](/images/Figure_3.png)
   
   (Takes the GC% content per read and shows a plot to visualize the counts distribution)
   ```
   import al-dente_fasta.counts_graphs as afcg
   afcg.gc_per_read(path/to/file_in/fastq.gz', 'path/to/file_out/gc_per_read'
   #kde_or_histo takes in values of 0/1/2, enter 0 to view a histogram plot, 1 for a kernal density estimate plot and 2 for a combination of both.
   afcg.counts_graph('path/to/file_out/gc_per_read', histo_or_kde, 'gc')
   ```

4. [Per base position DNA Base% content across all reads](/images/Figure_4.png)
   
   (Plot showing the distribution of each nucleotide base across base positions for each read)
   ```
   import al-dente_fasta.ntd_bases as afnb
   afnb.base_graphs('path/to/file_in/fastq.gz'
   

   ```
   
5. [Per base position N base% content across all reads](/images/Figure_5.png)
    
   (Plot showing the distribution of N (unknown base where sequencer cannot determine base or encounters an error) across base positions for each read)
   ```
   import al-dente_fasta.ntd_bases as afnb
   afnb.n_graph('path/to/file_in/fastq.gz'
   ```
   
6. [Per base position average quality score across each tile heatmap](/images/Figure_6.png)
    
   (Illumina sequencers only: shows heatmap of average quality scores of base positions for each read according to their tile number.)
   ```
   import al-dente_fasta.tile_seqs as afts
   afts.tile_num('path/to/file_in/fastq.gz', 'path/to/file_out/tile_num')
   #enc is the value used to encode quality scores. 
   #most sequencers use a 33 encoding code, as a guide for Sanger and Illumina 1.8+ sequencers > use 33, for Illumina 1.0-1.8 sequencers > use 64.
   result=afts.tile_link('path/to/file_out/tile_num', 'path/to/file_out/qual_out', enc')
   afts.heat_map(result)
   ```

### Benchmarks

Fastq files are very large files even when compressed as they can contain information about millions of reads, therefore processing them can take some time.
I have added a [benchmark](benchmarks.txt) file which gives guidelines how long each module/function can take for a [fastq](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001182785.1/)([SRR2093871_1.fastq.gz](https://www.ebi.ac.uk/ena/browser/view/PRJNA288953)) file containing ~7million reads.


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

[FastQC Documenation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)

<!-- PyPi -->
## PyPi
[PyPi](https://pypi.org/project/al-dente-fasta/0.1.0/#description)
