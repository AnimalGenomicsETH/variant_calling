# Variant calling

Some different workflows related to variant calling from existing genome alignments.

## Testing

Running with snakemake and polars

Example of snakemake_env.yaml

```
name: snake8
channels:
  - conda-forge
  - bioconda
  - nodefaults
dependencies:
  - snakemake
  - polars
```

We can then create and activate this conda environment, and test the overall workflow should run with

```
conda env create -f snakemake_env.yaml
conda activate snake8
cd .test
snakemake -s ../Snakefile --configfile config.yaml -n
```

This should print out that <TBD> jobs would be run, covering most pathways used in variant calling.
The snakemake rules were designed to run using our [snakemake profiles for Euler](https://github.com/AnimalGenomicsETH/euler_profiles).

## Usage

Many of the individual tools are assumed to be installed on $PATH, and so will almost certainly fail out of the box.
This is a work in progress to be more reproducible.

#### Detailed pipeline steps

![rulegraph](rulegraph.svg)

## Config

### Example snakemake command

Assuming you have the Euler profiles installed and are using snakemake v8+ with the slurm executor plugin, we can run

```
snakemake -s <path>/variant_calling/Snakefile --configfile config.yaml --profile "slurm/v8" --executor slurm --dry-run
```

## Roadmap

Several ideas for improvement.
Feel free to make pull requests with other solutions.

 - [ ] containerize almost every rule
 - [ ] add more general approach for imputation/filtering
 - [ ] re-integrate mendelian analyses
 
