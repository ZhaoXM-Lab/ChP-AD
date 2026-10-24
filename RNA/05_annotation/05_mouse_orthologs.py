#!/usr/bin/env python3
"""Cache BioMart homologues for the human CPEC marker panel."""

import os
from pathlib import Path
import pandas as pd
import json
from gseapy import Biomart

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
out = proj_dir / 'data/RNA/output/annotation/human_mouse_orthologs.tsv'
definition = json.loads((proj_dir / 'RNA/00_reference/cell_types.json').read_text())
genes = sorted({g for panel in definition['markers'].values() for g in panel})
result = Biomart().query(dataset='hsapiens_gene_ensembl', attributes=['external_gene_name','mmusculus_homolog_associated_gene_name'],
                         filters={'external_gene_name':genes})
result.dropna().drop_duplicates().to_csv(out, sep='\t', index=False)
