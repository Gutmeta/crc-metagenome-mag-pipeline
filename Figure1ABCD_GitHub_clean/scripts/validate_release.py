#!/usr/bin/env python3
"""Validate the release, its scientific summaries and an optional reproduction."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET
from zipfile import ZipFile

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
from PIL import Image, ImageChops
from scipy.stats import false_discovery_control

ROOT = Path(__file__).resolve().parents[1]
SOURCE_TABLES = [f"Source_Data_Figure1{panel}.tsv" for panel in "ABCD"] + ["Source_Data_Figure1A_heterogeneity.tsv"]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_manifest():
    manifest = pd.read_csv(ROOT / "MANIFEST.tsv", sep="\t")
    if manifest.file.duplicated().any():
        raise ValueError("Duplicate manifest entries")
    for row in manifest.itertuples():
        path = ROOT / row.file
        if not path.is_file() or path.stat().st_size != row.bytes or sha256(path) != row.sha256:
            raise ValueError(f"Checksum/size mismatch: {row.file}")
    covered = set(manifest.file)
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
              if p.is_file() and p.name != 'MANIFEST.tsv'
              and not any(part in ('reproduced', '__pycache__', '.venv', '.git') for part in p.relative_to(ROOT).parts)}
    if covered != actual:
        raise ValueError("Unmanifested or missing release files")


def check_portability():
    # Detect concrete personal absolute paths and common credential formats;
    # pattern definitions themselves contain no actual personal identifiers.
    patterns = [r'/(?:home|mnt|Users)/[A-Za-z0-9_.-]+/',
                r'[A-Za-z]:\\(?:Users|Documents)\\',
                r'gh[pousr]_[A-Za-z0-9]{30,}', r'AKIA[A-Z0-9]{16}',
                r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----']
    paths = [p for p in ROOT.rglob('*') if p.is_file() and not any(
        part in ('reproduced', '__pycache__', '.venv', '.git') for part in p.relative_to(ROOT).parts)]
    for path in paths:
        payloads = []
        if path.suffix in ('.xlsx', '.docx'):
            with ZipFile(path) as archive:
                payloads = [(name, archive.read(name).decode('utf-8', errors='ignore'))
                            for name in archive.namelist() if name.endswith(('.xml', '.rels'))]
        elif path.suffix == '.zip':
            with ZipFile(path) as archive:
                for name in archive.namelist():
                    if name.endswith(('.xlsx', '.docx')):
                        with ZipFile(io.BytesIO(archive.read(name))) as office:
                            payloads.extend((f'{name}:{member}', office.read(member).decode('utf-8', errors='ignore'))
                                            for member in office.namelist() if member.endswith(('.xml', '.rels')))
                    elif name.endswith(('.md', '.tsv', '.py', '.json', '.txt')):
                        payloads.append((name, archive.read(name).decode('utf-8', errors='ignore')))
                    elif name.endswith('.tsv.gz'):
                        payloads.append((name, gzip.decompress(archive.read(name)).decode('utf-8', errors='ignore')))
        elif path.suffix in ('.py', '.md', '.tsv', '.txt', '.json', '.svg', '.pdf') or path.name == '.gitignore':
            payloads = [(path.name, path.read_bytes().decode('utf-8', errors='ignore'))]
        for name, text in payloads:
            if any(re.search(pattern, text) for pattern in patterns):
                raise ValueError(f"Portability/privacy review needed: {path.relative_to(ROOT)}:{name}")


def check_figure(folder):
    svg_path = folder / 'figures/Figure1_ABCD.svg'
    svg = svg_path.read_text()
    match = re.search(r'width="([0-9.]+)pt" height="([0-9.]+)pt"', svg)
    if not match:
        raise ValueError("SVG dimensions missing")
    np.testing.assert_allclose([float(v) for v in match.groups()], [170/25.4*72,55/25.4*72], atol=1e-5)
    pdf = (folder / 'figures/Figure1_ABCD.pdf').read_bytes().decode('latin1')
    media = re.search(r'/MediaBox\s*\[([^\]]+)\]', pdf)
    if not media:
        raise ValueError("PDF page size missing")
    np.testing.assert_allclose([float(v) for v in media.group(1).split()], [0,0,170/25.4*72,55/25.4*72],atol=1e-5)
    for extension, size, dpi in [('png',(2007,649),300),('tiff',(4015,1299),600)]:
        if extension == 'tiff' and not (folder / f'figures/Figure1_ABCD.{extension}').exists():
            continue
        with Image.open(folder / f'figures/Figure1_ABCD.{extension}') as image:
            if image.size != size:
                raise ValueError(f"Incorrect {extension} size")
            np.testing.assert_allclose(np.asarray(image.info['dpi'],dtype=float),[dpi,dpi],atol=.01)
    def inspect(element, inherited):
        style = dict(inherited)
        for key in ('fill','font-family','font-weight','font-style'):
            if key in element.attrib:
                style[key] = element.attrib[key]
        for item in element.get('style','').split(';'):
            if ':' in item:
                key,value = item.split(':',1)
                style[key.strip()] = value.strip()
        if element.tag.endswith('}text'):
            if style.get('font-family') != "'Arial'" or style.get('fill','black') not in ('black','#000000'):
                raise ValueError("Text must be black Arial")
            if style.get('font-weight','normal') not in ('normal','400') or style.get('font-style','normal') != 'normal':
                raise ValueError("Text must be Regular")
        for child in element:
            inspect(child,style)
    inspect(ET.fromstring(svg),{})
    if any(x in svg.lower() for x in ('#f3f6f8','#f4f7f9','labels: observed / expected')):
        raise ValueError("Figure background or annotation does not match the specified layout")
    for text in ['I² = 4.3%', 'I² = 79.4%', 'I² = 83.5%', '1.45×', '1.35×', '1.42×', '2.75×']:
        if text not in svg:
            raise ValueError(f"Missing figure annotation: {text}")


def check_ibs_selection(selected, backgrounds):
    evidence = pd.read_csv(ROOT/'input/ibs_species_selection.tsv', sep='\t')
    required = {
        'study_code', 'species', 'case_samples', 'control_samples', 'detected_samples',
        'detection_fraction', 'mean_relative_abundance', 'case_mean_relative_abundance',
        'control_mean_relative_abundance', 'log2_fold_change', 'p_value', 'bh_q_value',
        'rf_top50_selection_frequency', 'direction_consistency', 'tested_comparisons',
        'eligible_for_selection', 'in_figure1_set',
    }
    if set(evidence.columns) != required or len(evidence) != 465 or evidence.species.duplicated().any():
        raise ValueError('Incorrect IBS selection-evidence dimensions or fields')
    if not evidence.study_code.eq('Mars_2020_IBS').all() or not evidence.species.str.startswith('s__').all():
        raise ValueError('Unexpected IBS study or taxon identifiers')
    if set(evidence.species) != backgrounds['IBS_Mars_2020']:
        raise ValueError('IBS evidence and detection background do not match')
    if not evidence.case_samples.eq(323).all() or not evidence.control_samples.eq(151).all():
        raise ValueError('Incorrect IBS comparison sample counts')
    samples = evidence.case_samples + evidence.control_samples
    np.testing.assert_allclose(evidence.detection_fraction, evidence.detected_samples/samples, atol=1e-12)
    if not evidence.detected_samples.between(0, samples).all():
        raise ValueError('Invalid detection counts')
    passes_filter = ((evidence.detection_fraction >= .05) | (evidence.detected_samples >= 10)) & (
        evidence.mean_relative_abundance >= 1e-5)
    if not passes_filter.all():
        raise ValueError('IBS background includes species that fail the abundance filter')
    for field in ('p_value', 'bh_q_value', 'rf_top50_selection_frequency', 'direction_consistency'):
        if not evidence[field].between(0, 1).all():
            raise ValueError(f'Invalid IBS evidence values: {field}')
    np.testing.assert_allclose(false_discovery_control(evidence.p_value), evidence.bh_q_value, rtol=1e-12, atol=1e-12)
    effect = np.log2((evidence.case_mean_relative_abundance+1e-6)/(evidence.control_mean_relative_abundance+1e-6))
    np.testing.assert_allclose(effect, evidence.log2_fold_change, rtol=1e-12, atol=1e-12)
    frequency = evidence.rf_top50_selection_frequency
    np.testing.assert_allclose(frequency*5, np.round(frequency*5), atol=1e-12)
    eligible = (evidence.direction_consistency >= .60) & (evidence.tested_comparisons >= 1) & (frequency > 0)
    for field in ('eligible_for_selection', 'in_figure1_set'):
        if evidence[field].dtype != bool:
            raise ValueError(f'IBS selection flag must be boolean: {field}')
    np.testing.assert_array_equal(eligible, evidence.eligible_for_selection)
    if int(eligible.sum()) != 81 or set(evidence.loc[eligible, 'species']) != selected['IBS_Mars_2020']:
        raise ValueError('IBS eligibility and figure membership do not match')
    np.testing.assert_array_equal(evidence.in_figure1_set, eligible)


def check_statistics(folder):
    source = folder / 'source_data'
    workbook = source / 'Source_Data_Figure1.xlsx'
    for name in SOURCE_TABLES:
        table = pd.read_csv(source/name,sep='\t')
        panel = name.removeprefix('Source_Data_Figure1').removesuffix('.tsv')
        sheet = 'Figure 1A heterogeneity' if panel == 'A_heterogeneity' else f'Figure 1{panel}'
        pd.testing.assert_frame_equal(table, pd.read_excel(workbook,sheet_name=sheet), check_dtype=False)
    auc = pd.read_csv(source/'Source_Data_Figure1A.tsv',sep='\t')
    original_input = pd.read_csv(ROOT/'input/panel_a_lodo_auc.tsv',sep='\t')
    columns = ['disease','project_accession']
    pd.testing.assert_frame_equal(auc.sort_values(columns).reset_index(drop=True),
                                  original_input.sort_values(columns).reset_index(drop=True))
    if len(auc) != 24 or auc.project_accession.nunique() != 19:
        raise ValueError('Unexpected AUC dimensions')
    memberships = pd.read_csv(ROOT/'input/species_sets.tsv',sep='\t')
    backgrounds = pd.read_csv(ROOT/'input/species_detection_backgrounds.tsv',sep='\t')
    if memberships.duplicated().any() or backgrounds.duplicated().any():
        raise ValueError('Duplicate species memberships')
    selected = {n:set(g.species) for n,g in memberships.groupby('analysis_set')}
    bg = {n:set(g.species) for n,g in backgrounds.groupby('analysis_set')}
    if {k:len(v) for k,v in selected.items()} != {'CRC_CRA':300,'IBD':300,'IBS_Mars_2020':81}:
        raise ValueError('Incorrect species-set sizes')
    if {k:len(v) for k,v in bg.items()} != {'CRC_CRA':1813,'IBD':1357,'IBS_Mars_2020':465}:
        raise ValueError('Incorrect backgrounds')
    for name in selected:
        if not selected[name] <= bg[name]:
            raise ValueError('Selected species outside their background')
    check_ibs_selection(selected, bg)
    overlap = pd.read_csv(source/'Source_Data_Figure1C.tsv',sep='\t')
    for row in overlap.itertuples():
        names = row.comparison_key.split('|')
        common = len(set.intersection(*(bg[n] for n in names)))
        expected = common*np.prod([len(selected[n])/len(bg[n]) for n in names])
        np.testing.assert_allclose(row.theoretical_expected,expected)
        if row.observed_shared_species != len(set.intersection(*(selected[n] for n in names))):
            raise ValueError('Observed overlap mismatch')
        if row.own_background_sizes != '|'.join(str(len(bg[n])) for n in names):
            raise ValueError('Background sizes mismatch')
    if not overlap.permutation_iterations.eq(100000).all() or not overlap.random_seed.eq(20260903).all():
        raise ValueError('Unexpected simulation settings')
    with np.load(folder/'analysis/overlap/Overlap_null_counts.npz') as archive:
        counts = archive['counts']
        if counts.shape != (100000,4) or list(archive['comparisons']) != overlap.comparison.tolist():
            raise ValueError('Null-count dimensions/order mismatch')
    np.testing.assert_allclose(counts.mean(axis=0),overlap.permutation_expected)
    np.testing.assert_allclose(np.quantile(counts,.95,axis=0),overlap.permutation_95th_percentile)
    exceed = (counts >= overlap.observed_shared_species.to_numpy()).sum(axis=0)
    np.testing.assert_array_equal(exceed,overlap.exceedance_count)
    np.testing.assert_allclose((1+exceed)/100001,overlap.permutation_p)
    np.testing.assert_allclose(false_discovery_control(overlap.permutation_p),overlap.permutation_bh_q)
    np.testing.assert_array_equal((overlap.permutation_bh_q<.05)&(overlap.observed_shared_species>overlap.permutation_95th_percentile),overlap.significant_enrichment)
    summaries = pd.read_csv(source/'Source_Data_Figure1B.tsv',sep='\t')
    for row in summaries.itertuples():
        values = auc.loc[auc.disease.eq(row.disease),'auc'].to_numpy()
        np.testing.assert_allclose(row.standard_deviation_auc,values.std(ddof=1))
        np.testing.assert_allclose(row.interquartile_range,np.subtract(*np.percentile(values,[75,25])))
    annotated = pd.read_csv(source/'Source_Data_Figure1A_heterogeneity.tsv',sep='\t')
    full = pd.read_csv(folder/'analysis/auc_heterogeneity/AUC_exploratory_heterogeneity.tsv',sep='\t')
    pd.testing.assert_frame_equal(annotated,full[full.effect_scale.eq('logit_auc')].reset_index(drop=True))
    from reproduce import exploratory_statistics, pairwise_matrix
    recomputed,_ = exploratory_statistics(auc)
    pd.testing.assert_frame_equal(full,recomputed,check_exact=False,rtol=1e-12,atol=1e-12)
    pd.testing.assert_frame_equal(pd.read_csv(source/'Source_Data_Figure1D.tsv',sep='\t'),pairwise_matrix(selected))
    from docx import Document
    legend = (folder/'docs/Figure1_legend.md').read_text()
    documents = [legend]
    word_path = folder/'docs/Figure1_methods_and_legend.docx'
    if word_path.exists():
        documents.append('\n'.join(p.text for p in Document(word_path).paragraphs))
    for text in documents:
        if 'exploratory' not in text or '100,000' not in text:
            raise ValueError('Method description missing')
        for row in overlap.itertuples():
            if f'{row.permutation_expected:.2f}' not in text or f'{row.permutation_bh_q:.6f}' not in text:
                raise ValueError('Document statistics mismatch')


def compare_reproduction(folder):
    for name in ['figures/Figure1_ABCD.tiff', 'docs/Figure1_methods_and_legend.docx']:
        if not (folder/name).is_file():
            raise ValueError(f'Missing generated artifact: {name}')
    check_figure(folder)
    check_statistics(folder)
    for name in SOURCE_TABLES:
        pd.testing.assert_frame_equal(pd.read_csv(ROOT/'source_data'/name,sep='\t'),
                                      pd.read_csv(folder/'source_data'/name,sep='\t'),
                                      check_exact=False,rtol=1e-12,atol=1e-12)
    with np.load(ROOT/'analysis/overlap/Overlap_null_counts.npz') as a, np.load(folder/'analysis/overlap/Overlap_null_counts.npz') as b:
        np.testing.assert_array_equal(a['counts'],b['counts'])
    for ext in ['png','tiff']:
        if not (ROOT/f'figures/Figure1_ABCD.{ext}').exists():
            continue
        with Image.open(ROOT/f'figures/Figure1_ABCD.{ext}') as a, Image.open(folder/f'figures/Figure1_ABCD.{ext}') as b:
            if ImageChops.difference(a.convert('RGB'),b.convert('RGB')).getbbox() is not None:
                raise ValueError(f'{ext} pixels differ; check dependency and Arial versions')
    print('Reproduction matches: source tables, 100,000 null draws and reference-image pixels.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reproduced',type=Path,help='Optional regenerated output directory.')
    args=parser.parse_args()
    check_manifest()
    check_portability()
    check_figure(ROOT)
    check_statistics(ROOT)
    if args.reproduced:
        compare_reproduction(args.reproduced.resolve())
    print('Release checks passed: checksums, portability, statistics, typography, dimensions and documents.')


if __name__ == '__main__':
    main()
