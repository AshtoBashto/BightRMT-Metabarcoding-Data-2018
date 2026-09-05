# Bight '18 RMP zooplankton DNA metabarcoding

Data and analysis code for Bandy et al., "Zooplankton DNA metabarcoding supports
regional ocean acidification monitoring and identifies novel candidate
bioindicators in the Southern California Bight" (PLOS ONE, PONE-D-26-23798).

CO1 (Leray primer) metabarcoding of mesozooplankton (>200 µm bongo net tows)
across 20 Southern California Bight sites and four seasons (Spring, Summer, Fall
2019; Winter 2020), collected under the Southern California Bight 2018 Regional
Monitoring Program.

## Contents

- **`code/`** – R analysis scripts that take the REVAMP ASV table through to
  the manuscript figures, tables, and statistics. See
  [`code/README.md`](code/README.md) for the pipeline, how to run it, and a
  description of every data file.
- **`data/`** – the processed data the scripts read and write.
- **`figures/`** – output location for the scripts (starts empty).
- **`KronaPlot.html`** – interactive Krona plot of the mesozooplankton taxonomy,
  viewable at
  <https://ashtobashto.github.io/bandy-etal-2026-zooplankton-metabarcoding-oa/KronaPlot.html>
  (made with <https://github.com/marbl/Krona>).

## Carbonate chemistry

The aragonite saturation state values used throughout the analysis
(`avgArag` and related columns in `data/Metadata_aligned.csv`) were estimated by
co-author C.A. Frieder from CTD temperature and oxygen using a regression trained
on discrete bottle samples, following Alin et al. (2012), with additional
errors-in-variables uncertainty propagation by N. Lombardo. The estimation code
and the underlying CTD and bottle-chemistry data (Southern California Bight 2018
RMP monitoring data) are archived separately; contact the corresponding author or
SCCWRP for access.

## Reference libraries and raw sequences

Raw sequence reads and FAIRe-formatted sample and molecular metadata are being
deposited separately (NCBI SRA and/or OBIS); accession information will be
added here once that submission is complete.

## License

Code and data in this repository are released under the MIT License (see
[`LICENSE`](LICENSE)).

## Citation

If accepted, please cite the paper. This repository is archived on Zenodo:
[10.5281/zenodo.22314208](https://doi.org/10.5281/zenodo.22314208) (resolves to
the latest version; `v1.0.0` is
[10.5281/zenodo.22314209](https://doi.org/10.5281/zenodo.22314209)).
