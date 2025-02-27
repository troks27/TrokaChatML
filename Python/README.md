# TrokaChatML

## Analysis Pipeline for Single-Cell RNAseq to Drug Target Discovery

TrokaChatML is an integrated computational framework designed for analyzing single-cell RNA sequencing (scRNAseq) data to uncover perturbed cell--cell communications and identify potential therapeutic targets. The pipeline combines differential expression analysis, machine learning-based noise filtration and perturbed communication prediction, tensor decomposition, and knowledge graph--driven drug target discovery, providing a robust and scalable approach from raw data processing to actionable biological insights.

## Installation

```bash
# Install from GitHub
pip install git+https://github.com/miketroka/TrokaChatML.git#subdirectory=Python
```

## Modules

TrokaChatML consists of the following modules:

- **trokachat_pathway_drug_discovery**: Tools for pathway analysis and drug discovery
- **trokachat_tensor_decomp**: Tensor decomposition methods for multi-dimensional data analysis
- **trokachat_perturb_model**: Perturbation modeling and analysis tools
- **trokachat_noise_model**: Noise filtering and modeling tools
- **trokachat_computation**: Core computational utilities

## Usage

```python
# Import specific modules
from trokachatml import trokachat_pathway_drug_discovery
from trokachatml import trokachat_tensor_decomp

# Use specific functions
from trokachatml.trokachat_pathway_drug_discovery import analyze_pathway
```

Refer to the documentation in the `docs/` directory for detailed usage examples and vignettes.

## License

This software is provided under the PENN ACADEMIC SOFTWARE LICENSE AGREEMENT. See the LICENSE file for details.

## Authors

- Michael Troka (troks27@upenn.edu)
- Dana Graves (University of Pennsylvania)
- Michael Gonzalez (University of Pennsylvania)

## Citation

If you use TrokaChatML in your research, please cite:

```
Troka M, Graves D, Gonzalez M. TrokaChatML: An integrated computational framework for single-cell RNA sequencing analysis and drug target discovery. [Journal/Conference details]
```

## Contact

For questions, feedback, or issues, please contact:
Michael Troka (troks27@upenn.edu)