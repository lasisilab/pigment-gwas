import os
import sys

# make the scripts/ modules importable as top-level (gwas_catalog, harmonize, vizhelpers)
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "scripts"))
