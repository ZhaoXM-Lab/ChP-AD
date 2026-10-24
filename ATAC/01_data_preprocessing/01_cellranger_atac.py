#!/usr/bin/env python3
"""Quantify retained ATAC libraries with Cell Ranger ATAC 2.2.0."""

import os
from pathlib import Path
import subprocess
import pandas as pd

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
samples = pd.read_csv(proj_dir / 'ATAC/00_reference/samples.tsv', sep='\t', dtype=str)
binary = os.environ.get('CELLRANGER_ATAC', 'cellranger-atac')
work = proj_dir / 'data/ATAC/cellranger'
work.mkdir(parents=True, exist_ok=True)
if samples.empty:
    raise ValueError('ATAC sample manifest is empty.')
for row in samples.loc[samples.include == '1'].itertuples():
    subprocess.run([binary,'count','--id='+row.sample,'--sample='+row.fastq_sample,
                    '--reference='+str(proj_dir / row.reference),'--fastqs='+str(proj_dir / row.fastq_dir),
                    '--localcores='+os.environ.get('THREADS','8'),'--localmem='+os.environ.get('MEM_GB','64')], cwd=work, check=True)
