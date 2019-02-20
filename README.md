# AnnotateSpliceogenicity

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Description

This script uses the [ENIGMA][1] thresholds to add spliceogenicity predictions to your [VEP][2] output.

For this to work correctly, the input must be in [Tab-delimited format][3] and have been annotated with the `MaxEntScan` and `HGVSIntronOffset` [VEP plugins][4].

## Usage

```
usage: annotate_spliceogenicity.py [-h] [<in.tsv>] [<out.tsv>]

optional arguments:
  -h, --help  show this help message and exit

required arguments:
  <in.tsv>    Read tab-delimited values from FILE [stdin]
  <out.tsv>   Write tab-delimited values to FILE [stdout]
```

## Citation

Jannah Shamsani, Stephen H Kazakoff, Irina M Armean, Will McLaren, Michael T Parsons, Bryony A Thompson, Tracy A O’Mara, Sarah E Hunt, Nicola Waddell, Amanda B Spurdle; A plugin for the Ensembl Variant Effect Predictor that uses MaxEntScan to predict variant spliceogenicity, _Bioinformatics_, , bty960, <https://doi.org/10.1093/bioinformatics/bty960>

## Contact

- [Jannah Shamsani](mailto:Jan.Shamsani@qimrberghofer.edu.au)
- [Stephen Kazakoff](mailto:Stephen.Kazakoff@qimrberghofer.edu.au)

[1]: https://enigmaconsortium.org
[2]: https://www.ensembl.org/vep
[3]: https://ensembl.org/info/docs/tools/vep/vep_formats.html#tab
[4]: https://github.com/Ensembl/VEP_plugins
