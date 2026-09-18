#!/usr/bin/env python3
"""Validate frozen cohort species profiles and minimal analysis-group labels."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def main():
    manifest=pd.read_csv(ROOT/'MANIFEST.tsv',sep='\t')
    for row in manifest.itertuples():
        path=ROOT/row.file
        if not path.is_file() or len(path.read_bytes())!=row.bytes or hashlib.sha256(path.read_bytes()).hexdigest()!=row.sha256:
            raise ValueError(f'Checksum mismatch: {row.file}')
    registry=pd.read_csv(ROOT/'COHORT_TABLES.tsv',sep='\t',keep_default_na=False)
    arms=pd.read_csv(ROOT/'ANALYSIS_GROUP_COUNTS.tsv',sep='\t',keep_default_na=False)
    if len(registry)!=20 or registry.cohort_id.duplicated().any():
        raise ValueError('Expected 20 independent source cohorts')
    total_samples=0; total_rows=0
    for row in registry.itertuples():
        p=ROOT/row.profile_file
        if hashlib.sha256(p.read_bytes()).hexdigest()!=row.released_profile_sha256:
            raise ValueError(f'Profile hash mismatch: {row.cohort_id}')
        profile=pd.read_csv(p,sep='\t',dtype={'sample_alias':str,'clade_name':str})
        if profile.columns.tolist()!=['sample_alias','clade_name','rel_abund']:
            raise ValueError('Unexpected profile fields')
        if not profile.sample_alias.str.fullmatch('SMP_[0-9a-f]{20}').all():
            raise ValueError('Non-release sample identifier found')
        if not profile.clade_name.str.startswith('s__').all():
            raise ValueError('Non-species rank found')
        if not np.isfinite(profile.rel_abund).all() or profile.rel_abund.lt(0).any():
            raise ValueError('Invalid abundance value')
        groups=pd.read_csv(ROOT/row.sample_groups_file,sep='\t')
        allowed={'sample_id','CRA','CRC','CRC_CRA','IBD','IBS'}
        if not set(groups.columns)<=allowed or groups.sample_id.duplicated().any():
            raise ValueError('Unexpected grouping fields or duplicate identifiers')
        if set(profile.sample_alias)!=set(groups.sample_id):
            raise ValueError('Profile/group sample mismatch')
        for arm in arms[arms.cohort_id.eq(row.cohort_id)].itertuples():
            counts=groups[arm.disease].value_counts()
            if counts.get('case',0)!=arm.case_samples or counts.get('control',0)!=arm.control_samples:
                raise ValueError('Analysis group counts mismatch')
            if not set(groups[arm.disease])<={'case','control','not_in_comparison'}:
                raise ValueError('Unexpected sample group')
        if len(profile)!=row.released_profile_rows or profile.sample_alias.nunique()!=row.released_profile_samples or profile.clade_name.nunique()!=row.released_profile_species:
            raise ValueError('Profile dimensions mismatch')
        total_samples+=profile.sample_alias.nunique(); total_rows+=len(profile)
    expected=json.loads((ROOT/'VALIDATION_SUMMARY.json').read_text())
    if total_samples!=expected['released_profile_samples'] or total_rows!=expected['released_profile_rows']:
        raise ValueError('Total dimensions mismatch')
    print(f'Validated 20 cohorts; {total_samples} profiled samples; {total_rows} abundance rows.')


if __name__=='__main__':
    main()
