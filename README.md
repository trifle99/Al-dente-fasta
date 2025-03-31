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
I created a python version of FastQC where this package 'Al-Dente Fasta' contains modules and functions for quality control checking on fastq sequencing datasets.

Currently, this tool only works on fastq files where all reads are the same length.
The analysis it provides are:
* Basic summary (total base counts, total read counts, total gc%)
* Per base read quality scores box plot
* Per read average quality scores plot
* Per base read content plot 
* Per read GC% content plot
* Per tile read avg quality scores plot (for illumina sequences only as they contain extra information on their sequence id's)

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* [FastQC Documenation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)
