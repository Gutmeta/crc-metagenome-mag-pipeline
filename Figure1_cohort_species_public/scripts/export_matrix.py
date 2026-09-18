#!/usr/bin/env python3
"""Convert a frozen long-format cohort profile to a species-by-sample matrix."""
from pathlib import Path
import argparse
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cohort',required=True,help='cohort_id from COHORT_TABLES.tsv')
    parser.add_argument('--output',required=True,type=Path,help='New TSV or TSV.GZ output path')
    args=parser.parse_args()
    registry=pd.read_csv(ROOT/'COHORT_TABLES.tsv',sep='\t')
    found=registry[registry.cohort_id.eq(args.cohort)]
    if len(found)!=1:
        parser.error('Unknown or ambiguous cohort_id')
    if args.output.exists():
        parser.error('Output exists; choose a new path to preserve frozen inputs')
    frame=pd.read_csv(ROOT/found.iloc[0].profile_file,sep='\t')
    # This is the upstream loading convention: absent sample-species pairs are zero,
    # and repeated sample-species rows are summed. No renormalization or filtering.
    matrix=frame.pivot_table(index='clade_name',columns='sample_alias',values='rel_abund',aggfunc='sum',fill_value=0.)
    matrix.index.name='species'
    args.output.parent.mkdir(parents=True,exist_ok=True)
    matrix.to_csv(args.output,sep='\t',compression='infer')
    print(f'Wrote {matrix.shape[0]} species × {matrix.shape[1]} samples.')


if __name__=='__main__':
    main()
