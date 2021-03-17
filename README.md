# motevowrapper

Simple Python parser for MotEvo files.

To install, run:

```bash
pip install motevowrapper
```

[![Tests Ubuntu](https://github.com/brlauuu/motevowrapper/workflows/Tests/badge.svg)](https://github.com/brlauuu/motevowrapper/actions/workflows/python-tests.yml)
[![Tests MacOS](https://github.com/brlauuu/motevowrapper/actions/workflows/python-tests-mac.yml/badge.svg)](https://github.com/brlauuu/motevowrapper/actions/workflows/python-tests-mac.yml)
[![Upload Python Package](https://github.com/brlauuu/motevowrapper/actions/workflows/python-publish-to-pypi.yml/badge.svg)](https://github.com/brlauuu/motevowrapper/actions/workflows/python-publish-to-pypi.yml)

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

## Running MotEvo from `motevowrapper`

Method for running MotEvo is `run_motevo(...)`. Parameter description can be found in the [MotEvo source code](https://swissregulon.unibas.ch/sr/software). I copied it here for better visibility.

```python
sites_file, priors_file = mw.run_motevo(
    # Command line parameters
    sequences_file,             # Sequences or alignments file
    wm_path,                    # Path to the position weight matrix (PWM) of a given motif
    working_directory="./",     # Working directory

    # General
    Mode="TFBS",                # (word) Mode of running. This can be TFBS (TFBS predictions; default), ENH (enhancer finding), or WMREF (weight matrix refinement)
    refspecies,                 # (word) The identifier of the reference species (as found in the sequence identifier and in the phylogenetic tree).
    TREE,                       # (tree string) Phylogenetic tree in Newick format.
    restrictparses,             # (binary) When 1 only use sites that have a reference weight matrix score bigger than 0. Default: 0. Only used for testing.
    singlestrand,               # (binary) When 1 predict sites only on the positive strand.

    # Priors
    bgprior,                    # (real number) Prior probability for putting down a background at each position.
    EMprior=0,                  # (binary) Use the expectation maximization algorithm to find the priors that maximize the probability of the observed alignment.
    priordiff,                  # (real value) Convergence criterion for prior estimation, e.g. at 0.01 iteration stops when priors change by less than 1%.
    UFEwmprior,                 # (real number) The prior weight of the UFE model relative to the other weight matrices.

    # Background model
    markovorderBG,              # (integer) Markov order of the background model.
    bgA=0.25,                   # (real number) background probability for A (for the zeroth order model)
    bgC=0.25,                   # (real number) background probability for C (for the zeroth order model)
    bgG=0.25,                   # (real number) background probability for G (for the zeroth order model)
    bgT=0.25,                   # (real number) background probability for T (for the zeroth order model)
    mybgfile,                   # (file name) Input file containing a higher order background model.

    # UFE
    UFEwmfile,                  # (file name) Input file containing the UFE model (run 'runUFE' to create it for a given tree and background model.)
    UFEwmlen,                   # (integer) The length of UFE model sites.
    UFEprint,                   # (binary) When set to zero UFE sites are not reported in the site file.
    UFEwmproffile,              # (file name) Output file containing UFE model probabilities at each position.

    # TFBS output
    sitefile,                   # (file name) Output file name of the file containing the predicted sites.
    priorfile,                  # (file name) Output file containing information like site density, final priors, and the total number of sites for each WM.
    loglikfile,                 # (file name) Output file containing log-likelihood of each sequence (or alignment) in the input data.
    minposterior=0.1,           # (real number) When printing sites, only print sites with a posterior bigger than this cut-off.
    printsiteals=1,             # (binary) When set to zero sequence alignments are not printed in the output file.

    # WM refinement
    minposteriorWM,             # (real number) When doing weight matrix refinement, only include sites to refine that have a minimal posterior bigger than this cut-off.
    wmdiff,                     # (real number) Convergence criterium for WM refinement, e.g. at 0.01 iteration stops when WM entries change by less than 1%

    # Enhancer prediction
    CRMfile,                    # (file name) Output file containing the results when running MotEvo in the enhancer prediction mode.
    winlen,                     # (integer) Length of the enhancer window used in enhancer prediction mode.
    steplen,                    # (integer) Number of positions by which the window is moved at each step during enhancer prediction.

    # Additional parameters
    try_until_succeeding=False, # Run MotEvo until there `sites` and `priors` files are created
    verbose=False,              # Print more details during MotEvo run
)
```

You can note two parameters were added, `try_until_succeeding` and `verbose`. These were added for the needs of this Python wrapper.

Parameters that have default value set, will be used for sure, including:

* `TREE` which is set to species tree in case phylogenetic tree is not provided.
* `UFEwmlen` which is set the length of PWM in use in case `"auto"` is passed to this parameter.

For example, in order to use it you can use the following example:

```python
import motevowrapper.motevowrapper as mw

sites_file, priors_file = run_motevo(
    sequences_file="zebrafish_promoters.fa",
    working_directory="./",
    wm_path="REST.wm",
    tree="(danRer11: 1.0);",
    ref_species="danRer11",
    background_prior=0.8,
)

```

For more information on how to use all of MotEvo's options, please check out [MotEvo source code](https://swissregulon.unibas.ch/sr/software) and [MotEvo paper](https://pubmed.ncbi.nlm.nih.gov/22334039/).

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
