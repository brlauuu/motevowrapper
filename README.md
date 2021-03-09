![Tests](https://github.com/brlauuu/motevowrapper/workflows/Tests/badge.svg) 
[![Upload Python Package](https://github.com/brlauuu/motevowrapper/actions/workflows/python-publish-to-pypi.yml/badge.svg?event=release)](https://github.com/brlauuu/motevowrapper/actions/workflows/python-publish-to-pypi.yml)

# motevowrapper

Simple Python parser for MotEvo files.

To install, run:

```bash
pip install motevowrapper
```

## MotEvo

[MotEvo](https://pubmed.ncbi.nlm.nih.gov/22334039/) (Arnold et al. 2012) is a Bayesian probabilistic model for prediction of transcription factor binding sites (TFBSs) for a given set of position weight matrices (PWMs) and DNA sequences. It was developed by van Nimwegen lab at the Biozentrum (University of Basel, Switzerland) and it can be acquired [here](https://swissregulon.unibas.ch/sr/software).

This repository contains the source code for a simple Python package that allows you to:

1. Run MotEvo with given parameters
2. Parse MotEvo output files
3. Visualize visualize site density per motif

## Installing MotEvo

MotEvo source code can be downloaded from the [SwissRegulon website](https://swissregulon.unibas.ch/sr/software). You can either download the source and compile it, or download binaries for MacOS or Linux. Don't forget to add path to executables to your `.bashrc` or `.bash_profile`. You can do it by simply running

```bash
export PATH=$PATH:/path/to/motevo/bin
```

## Running MotEvo from MotevoWrapper

Method for running MotEvo is `run_motevo(...)`. Method description is the following:

```python
sites_file, priors_file = mw.run_mote(
    sequences_file=None,        # Or alignments file
    wm_path=None,               # Path to PWM file
    working_directory="./",     # (optional) Working directory where MotEvo will be ran
    mode="TFBS",                # Mode
    tree=None,                  # (optional) Tree if available, otherwise no other species are assumed
    ref_species=None,           # Reference species
    em_prior=None,              # Expectation Maximization prior
    ufe_wm_prior=None,          # (optional) UFE weight matrix prior. Essentially number of other possible PWMs. Tied to other 'ufe' parameters
    ufe_wm_file=None,           # (optional) UFE file. Generated by running "runUFE" on a given tree and nucleotide likelihood. Tied to other 'ufe' parameters
    ufe_wm_len=None,            # (optional) UFE weight matrix length. Length of other 'competing' motifs on given sequences. If set to 'auto', given motif length is used.
    background_prior=None,      # Background prior. 1-background_prior == probability that a given motif has a site on a given sequence
    bgA=0.25,                   # Probability of randomly seeing nucleotide 'A' in a given sequence
    bgT=0.25,                   # Probability of randomly seeing nucleotide 'T' in a given sequence
    bgG=0.25,                   # Probability of randomly seeing nucleotide 'G' in a given sequence
    bgC=0.25,                   # Probability of randomly seeing nucleotide 'C' in a given sequence
    sites_file=None,            # (optional) Path to storing 'sites' file. Default path is {working_directory}/sites_{motif}.wm
    priors_file=None,           # (optional) Path to storing 'priors' file. Default path is {working_directory}/priors_{motif}.wm
    print_site_als=1,           # (optional) Outputting alignments in the sites file.
    minposterior=0.1,           # Minimal posterior allowed. Minimal probability that a motif will bind on a a given sequence. Only 1 site!
    try_until_succeeding=False, # Option that allows user to run MotEvo until it works. Sometimes MotEvo breaks due to memory allocation which is where this option might be useful.
    verbose=False,              # Verbose
)
```

For example, in order to use it you can use the following example:

```python
import motevowrapper.motevowrapper as mw

sites_file, priors_file = run_motevo(
    sequences_file="zebrafish_promoters.fa",
    working_directory="./",
    wm_path="REST.wm",
    tree="(danRer11: 1.0);",
    ref_species="danRer11",
    em_prior=0,
    background_prior=0.8,
)

```

## Parsing MotEvo files from MotevoWrapper

MotEvo produces 2 files: `sites` and `priors` file. Usage of the package is simple. For a given MotEvo sites file stored at `/path/to/sites_MOTIF.wm` by calling:

```python
import motevowrapper.motevowrapper as mw

df_sites = mw.parse_sites('/path/to/sites_file') # Motif binding sites
df_priors = mw.parse_sites('/path/to/priors_file') # Final file with priors

```

you get a Pandas data frame containing parsed data from the MotEvo run. Further manipulation with the dataframe allows getting motif binding density on all sequences, number of binding sites, number of different species from alignment used, etc.

## Visualizing site density per motif using MotevoWrapper

```python
df = mw.parse_sites("sites_REST.wm")
mw.plot_site_distribution("REST", df)
```

## References

1. Arnold, Phil, et al. "MotEvo: integrated Bayesian probabilistic methods for inferring regulatory sites and motifs on multiple alignments of DNA sequences." Bioinformatics 28.4 (2012): 487-494.
