![Tests](https://github.com/brlauuu/motevowrapper/workflows/Tests/badge.svg) ![Upload to PyPi](https://github.com/brlauuu/motevowrapper/workflows/Upload%20to%20PyPi/badge.svg)

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

## Running MotEvo from MotevoWrapper

TODO

## Parsing MotEvo files from MotevoWrapper

MotEvo produces 2 files: `sites` and `priors` file. Usage of the package is simple. For a given MotEvo sites file stored at `/path/to/sites_MOTIF.wm` by calling:

```python
import motevowrapper as mw

df_sites = mw.parse_sites('/path/to/sites_MOTIF.wm') # Motif binding sites
df_priors = mw.parse_sites('/path/to/priors_MOTIF.wm') # Final file with priors

```

you get a Pandas data frame containing parsed data from the MotEvo run. Further manipulation with the dataframe allows getting motif binding density on all sequences, number of binding sites, number of different species from alignment used, etc.

## Visualizing site density per motif using MotevoWrapper

TODO

## References

1. Arnold, Phil, et al. "MotEvo: integrated Bayesian probabilistic methods for inferring regulatory sites and motifs on multiple alignments of DNA sequences." Bioinformatics 28.4 (2012): 487-494.
