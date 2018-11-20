# isoseq3
This is a [nextflow](https://github.com/nextflow-io/nextflow) implementation of [Pacific Biosciences IsoSeq3 pipeline](https://github.com/PacificBiosciences/IsoSeq3).

## Setup
### Requirements
This pipeline requires:
1. [IsoSeq3](https://github.com/PacificBiosciences/IsoSeq3)
2. [STAR](https://github.com/alexdobin/STAR/)
3. [samtools](https://github.com/samtools/samtools) 
4. [bedtools2](https://github.com/arq5x/bedtools2)

### Configuration
The current pipeline is designed to run on the Mendel cluster of the [GMI Vienna](https://www.gmi.oeaw.ac.at/). To make it run for your group edit the `projectName`, `annotation` & `fasta` parameters in the the mendel.config file to fit to your group project and needs. To make it run on another infrastructure simply add a new nextflow config file in the conf folder and source via the nextflow.config file. See [here](https://www.nextflow.io/docs/latest/config.html) for more information. 

## Workflow of the pipeline
1. run ccs
2. run lima
3. cluster reads
4. polish reads
5. align with STARlong
6. convert to bed12

## Run the pipeline
To run the pipeline simply run e.g:

```bash
nextflow run_isoseq3.nf --input "/lustre/scratch/users/falko.hofmann/isoseq/*/"
