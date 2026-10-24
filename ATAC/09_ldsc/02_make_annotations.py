#!/usr/bin/env python3
"""Create binary SNP annotations in exactly the BIM row order."""

import os
from pathlib import Path
import pandas as pd
from pybedtools import BedTool

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
reference = proj_dir / 'data/reference/ldsc/1000G_EUR_Phase3_plink'
base = proj_dir / 'data/ATAC/output/ldsc'
out = base / 'annotations'
out.mkdir(parents=True, exist_ok=True)
beds = sorted((base / 'beds').glob('*.bed'))
if not beds:
    raise ValueError('Run 01_export_annotations.R first.')
for chrom in range(1,23):
    bim = pd.read_csv(reference / f'1000G.EUR.QC.{chrom}.bim', sep=r'\s+', header=None,
                      names=['CHR','SNP','CM','BP','A1','A2'])
    intervals = BedTool.from_dataframe(pd.DataFrame({'chr':'chr'+bim.CHR.astype(str),'start':bim.BP-1,'end':bim.BP,'row':bim.index}))
    for bed in beds:
        rows = [int(record[3]) for record in intervals.intersect(BedTool(str(bed)), u=True)]
        annot = bim[['CHR','BP','SNP','CM']].copy()
        annot['ANNOT'] = 0
        annot.loc[rows,'ANNOT'] = 1
        annot.to_csv(out / f'{bed.stem}.{chrom}.annot.gz', sep='\t', index=False, compression='gzip')
