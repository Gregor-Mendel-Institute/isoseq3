# isoseq3
This is a [nextflow](https://github.com/nextflow-io/nextflow) implementation of [Pacific Biosciences IsoSeq3.1 pipeline](https://github.com/PacificBiosciences/IsoSeq3/blob/master/README_v3.1.md).

## Setup
### Requirements
This pipeline requires anaconda[https://anaconda.org/]. The required dependencies will then be installed by nextflow into a conda virtual environment

### Configuration
The current pipeline is designed to run on the Mendel cluster of the [GMI Vienna](https://www.gmi.oeaw.ac.at/). To make it run for your group edit the `projectName` & `fasta` parameters in the the mendel.config file to fit to your group project and needs. To make it run on another infrastructure simply add a new nextflow config file in the conf folder and source via the nextflow.config file. See [here](https://www.nextflow.io/docs/latest/config.html) for more information. 

## Workflow of the pipeline
1. circular consensus calling
2. primer removal
3. refine reads
4. merge samples (optional)
3. cluster reads 
4. polish reads
5. align transcripts (optional)


## Run the pipeline
To run the pipeline simply run e.g:

```bash
nextflow run_isoseq3.nf --input "/lustre/scratch/users/falko.hofmann/isoseq/samples/*/" --output "/lustre/scratch/users/falko.hofmann/isoseq/results/*/ --primer_type default
