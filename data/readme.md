# Data for Tamipami paper

Data for this paper are available under NCBI Bioproject [PRJNA1298332 : Tamipami: Software for determining the PAM and TAM sites of new CRISPR/Cas and TnpB nucleases ](https://www.ncbi.nlm.nih.gov/bioproject/1298332)

## Download from NCBI

To download using [SRA toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit):

```{bash}
mkdir -p fastq

for i in $(seq 34761000 34761011); do
  SRR="SRR$i"
  fasterq-dump $SRR --split-files --outdir fastq
  gzip fastq/${SRR}_1.fastq
  gzip fastq/${SRR}_2.fastq
done
```

To use wget on, for example, SRR34761000:
```
wget -O SRR34761000.fastq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq\?acc\=SRR34761000
```


## Metadata files

- [`Model.organism.animal.1.0.xlsx`](metadata/Model.organism.animal.1.0.xlsx) the NCBI Biosample metadata file
- [`SRA-metadata-15476171-processed-ok.xlsx`](metadata/SRA-metadata-15476171-processed-ok.xlsx) The SRA metadata file with an the NCBI identifiers
- [`SRA_metadata.xlsx`](metadata/SRA_metadata.xlsx) SRA metadatafiel uploaded to NCBI for SRA deposition
- [`experiments.xlsx`](metadata/experiments.xlsx) A file describing which Sample files ere paired as paired and experimental input files for each experiment. This lists both the original sample names and the names of sequence files after re-downloading the data from SRA

## Simulation

- `sim_data.py` - a script to generate synthetic count data conforming to a negative binomial distribution. Used for initial development.