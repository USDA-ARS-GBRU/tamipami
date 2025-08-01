# Tamipami: A web application and command line interface for finding the TAM and PAM sites of novel endonucleases

__Carlos Orosco<sup>1</sup>, Pyush Jain<sup>1</sup>, Adam R. Rivers<sup>2</sup>__

1. University of Florida, Department of Chemical Engineering, Gainesville, FL USA
2. United States Department of Agriculture, Agricultural Research Service, Genomics and Bioinformatics Research Unit, Gainesville, FL USA


<br>

When a new Cas or TnpB guided endonuclease is discovered or engineered, one of the first tasks is identifying its PAM or TAM recognition site. This can be done by treating a pool of plasmid DNA containing a target site adjacent to a random region. The random regions containing  the PAM/TAM site are recognized and cut in the presence of the endonuclease and a guide RNA. These become depleted in the sequencing library. By comparing an uncut control library  to a cut  experimental library it is possible to identify the PAM/TAM site.  The method was introduced by [Walton et al. 2021]( https://doi.org/10.1038/s41596-020-00465-2).
That work deposited the plasmid pools with [Addgene]( https://www.addgene.org/pooled-library/kleinstiver-ht-pamda/), making the lab protocol accessible.   

This web application and command line application builds on the work by creating software that simplifies the analysis of the sequencing data and adds rich interactive visualizations for selecting the PAM/TAM site.

## Installation

The application can be installed into a conda environemnt using Pip and conda.

```{bash}
conda create -n tamipamienv -c conda-forge python=3.13
conda activate tamipamienv
conda install -c bioconda bbmap
git clone git@github.com:USDA-ARS-GBRU/tamipami.git
cd tamipami
pip install .

```

## Web Application 

The web application can be found at [https://tamipami.che.ufl.edu](https://tamipami.che.ufl.edu)

To launch a the streamlit web application locally at http://localhost:8501 run this commmand:

```{bash}
streamlit run tamipami/app.py 
```


Tamipami Inputs | Tamipami results
----------------|------------------
![Tamipami App](/tamipami/assets/app_screenshot1.png) | ![Tamipami App output](/tamipami/assets/app_screenshot2.png)

### Web application use

1. Load the forward and reverse FASTQ files for a single control and experimental library in the four input boxes on the sidebar.
2. Select the Addgene library used for the experiment. Or...
3. If no library is selected, enter the target sequence and the orientation (5prime is PAM-Target, 3prime is Target-PAM (like spCas9)).
4. Select the maximum length to analyze. You can compare all smaller lengths once you have analyzed the data
5. Hit 'Submit'. This will process your files.
6. After you have hit 'Submit' you can explore your data interactively. The key interface is the z-score slider bar. Moving that will set the 
cutoff value to separate kmers that cut from those that did not. This will update the histogram of sequence z-scores, the table of reads the degenerate sequences created and the sequence motif.
7. Be sure to explore each length tab, to select the PAM/TAM site best supported by your data.
8. Export the raw run data.

## Command line program

To facilitate automated analysis we have also released a command line version of TamiPami.  The CLI has two subcommands: `process` and `predict`. The first step is to process the raw FASTQ data into counts of kmers. This data can be sent ot a JSON file or STDOUT.  The second step is to predict the PAM Sequences given a specific cutoff value. The commands can be piped together.

```
TamiPami: a CLI application to parse High throughput PAM/TAM site sequencing data

positional arguments:
  {process,predict}
    process          Sub-command "process": used to process FASTQ data into a summarized json output
    predict          Subcommand "predict" use to predict PAMs/TAMs and summary data for a selected length and cutoff value

options:
  -h, --help         show this help message and exit
  --log LOG          Log file
  -V, --version      show program's version number and exit
```

### Process

```
usage: tamipami process [-h] --cont1 CONT1 --cont2 CONT2 --exp1 EXP1 --exp2 EXP2 [--outfile OUTFILE] [--library {RTW554,RTW555,RTW572,RTW574}] [--spacer SPACER] [--orientation {3prime,5prime}]
                        [--length [3-8]]

options:
  -h, --help            show this help message and exit
  --cont1, -c CONT1     A forward .fastq, .fq, .fastq.gz or .fq.gz file. .
  --cont2, -c2 CONT2    A reverse .fastq, .fq, .fastq.gz or .fq.gz file.
  --exp1, -e EXP1       A forward .fastq, .fq, .fastq.gz or .fq.gz file.
  --exp2, -e2 EXP2      A reverse .fastq, .fq, .fastq.gz or .fq.gz file.
  --outfile OUTFILE     A path to a hdf5 file of the results, if missing binary data will be sent to STDOUT.
  --library {RTW554,RTW555,RTW572,RTW574}
                        The Addgene library pool. For custom pools use the --spacer and --orientation flags
  --spacer SPACER       The spacer sequence for the guide RNA. Not needed if ---library is used
  --orientation {3prime,5prime}
                        the side of the spacer the PAM/TAM is on
  --length [3-8]        The length of the PAM or TAM sequences
  ```

  ### Predict

```
usage: tamipami predict [-h] [--input INPUT] [--cutoff CUTOFF] --predict_out PREDICT_OUT

options:
  -h, --help            show this help message and exit
  --input INPUT         An file containing the data from a TamiPami process run or downloaded data from the web ap, if not input is provided STDIN will
                        be assumed
  --cutoff CUTOFF       A json string containing the kmer lengths and the Zscore cutoff values above which kmers are considered part of the PAM/TAM.
                        Single and double quotes are required. If no cutoff is provided it will be automatically calculated using univariate k means
                        clustering. Example input: --cutoff '{"3": 2, "4": 2, "5": 2, "6": 2}'
  --predict_out PREDICT_OUT
                        A file directory containing the PAM/TAM and degenerate sequences identified
```

## Example data

Example datasets are available under NCBI Bioproject [PRJNA1298332 : Tamipami: Software for determining the PAM and TAM sites of new CRISPR/Cas and TnpB nucleases ](https://www.ncbi.nlm.nih.gov/bioproject/1298332). For more information see [`data/readme.md`](`data/readme.md).


## Performance

Detailed performance metrics for 7 datasets are found at [`benchmark/readme.md`](benchmark/readme.md). TL;DR most datasets take 15-30 seconds to process.

## Need more help?

If you need help or encounter errors please request support on the [Tamipami Github issues page](https://github.com/USDA-ARS-GBRU/tamipami/issues)   
    

## Server notes

Docker deployment command:

```
docker run -d \
  -p 8501:8501 \
  --name tamipami-app \
  streamlit-app:latest \
  streamlit run /app/tamipami/app.py --server.port=8501 --server.address=0.0.0.0
```

## License information

This software is a work of the United States Department of Agriculture,
Agricultural Research Service and is not copyrightable under U.S. Code Title 17, Section 105. It is released under a Creative Commons CC0
public domain attribution.