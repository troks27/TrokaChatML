"""
TrokaChatML: Analysis Pipeline for Single-Cell RNAseq to Drug Target Discovery

TrokaChatML is an integrated computational framework designed for analyzing 
single-cell RNA sequencing (scRNAseq) data to uncover perturbed cell--cell
communications and identify potential therapeutic targets. The pipeline combines
differential expression analysis, machine learning-based noise filtration and
perturbed communication prediction, tensor decomposition, and knowledge
graph--driven drug target discovery, providing a robust and scalable approach
from raw data processing to actionable biological insights.
"""

__version__ = "0.0.0.9000"
__author__ = "Michael Troka, Dana Graves, Michael Gonzalez"

# Import submodules to make them accessible
try:
    from . import trokachat_pathway_drug_discovery
except ImportError:
    pass

try:
    from . import trokachat_tensor_decomp
except ImportError:
    pass

try:
    from . import trokachat_perturb_model
except ImportError:
    pass

try:
    from . import trokachat_noise_model
except ImportError:
    pass

try:
    from . import trokachat_computation
except ImportError:
    pass
