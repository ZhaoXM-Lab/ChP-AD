#!/usr/bin/env python3
"""Quantify retained human and mouse libraries with Cell Ranger 9.0.0."""

import os
from pathlib import Path
import subprocess
import pandas as pd

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
samples = pd.read_csv(proj_dir / 'RNA/00_reference/samples.tsv', sep='\t', dtype=str).fillna('')
binary = os.environ.get('CELLRANGER', 'cellranger')
cores = os.environ.get('THREADS', '16')
memory = os.environ.get('MEM_GB', '128')
if samples.empty:
    raise ValueError('RNA sample manifest is empty.')
for row in samples.loc[samples.include == '1'].itertuples():
    work = proj_dir / 'data/RNA/cellranger' / row.species
    work.mkdir(parents=True, exist_ok=True)
    subprocess.run([binary, 'count', '--id=' + row.sample,
                    '--sample=' + row.fastq_sample,
                    '--transcriptome=' + str(proj_dir / row.transcriptome),
                    '--fastqs=' + str(proj_dir / row.fastq_dir),
                    '--create-bam=true', '--include-introns=true',
                    '--localcores=' + cores, '--localmem=' + memory], cwd=work, check=True)
