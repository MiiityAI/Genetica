#!/usr/bin/env python3
"""
Universal SNP Extractor
=======================
Combines ClinVar + GWAS Catalog data with smart merge by rsID.

Features:
- Merges data from both sources (no duplicates)
- MAF filter (default ≥5%)
- ClinVar significance filters (pathogenic, risk_factor, drug_response)
- Full annotation: OR, p-value, SIFT, PolyPhen, frequencies, positions

Usage:
    python snp_extractor_v2.py --genes TYR,MC1R -o results.csv
    python snp_extractor_v2.py --kegg pathway.xml --min-maf 0.01 -o results.xlsx
    python snp_extractor_v2.py --genes APOE --pathogenic --risk-factor -o apoe.csv
"""

import argparse
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Dict, Optional, Set

import pandas as pd
import requests
from tqdm import tqdm

try:
    import mygene
    HAS_MYGENE = True
except ImportError:
    HAS_MYGENE = False


class UniversalSNPExtractor:
    """Combined ClinVar + GWAS extractor with smart merge"""
    
    GWAS_API = "https://www.ebi.ac.uk/gwas/rest/api"
    MYVARIANT_API = "https://myvariant.info/v1"
    
    def __init__(self, 
                 min_maf: float = 0.05,
                 p_threshold: float = 5e-8,
                 include_pathogenic: bool = False,
                 include_risk_factor: bool = False,
                 include_drug_response: bool = False,
                 verbose: bool = True):
        self.min_maf = min_maf
        self.p_threshold = p_threshold
        self.include_pathogenic = include_pathogenic
        self.include_risk_factor = include_risk_factor
        self.include_drug_response = include_drug_response
        self.verbose = verbose
        self.session = requests.Session()
        self.session.headers.update({"Accept": "application/json"})
    
    def log(self, msg: str):
        if self.verbose:
            print(msg)

    # =========================================
    # INPUT PARSING
    # =========================================
    
    def parse_kegg_xml(self, xml_path: str) -> List[str]:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        entrez_ids = []
        for entry in root.findall(".//entry[@type='gene']"):
            name = entry.get('name', '')
            ids = re.findall(r'hsa:(\d+)', name)
            entrez_ids.extend(ids)
        entrez_ids = list(set(entrez_ids))
        self.log(f"Found {len(entrez_ids)} Entrez IDs in {xml_path}")
        return self.entrez_to_symbols(entrez_ids)
    
    def parse_kegg_folder(self, folder_path: str) -> List[str]:
        all_genes = []
        for xml_file in Path(folder_path).glob("*.xml"):
            genes = self.parse_kegg_xml(str(xml_file))
            all_genes.extend(genes)
        return list(set(all_genes))
    
    def entrez_to_symbols(self, entrez_ids: list) -> List[str]:
        if not HAS_MYGENE:
            raise ImportError("mygene required")
        mg = mygene.MyGeneInfo()
        result = mg.querymany(entrez_ids, scopes='entrezgene',
                              fields='symbol', species='human',
                              as_dataframe=True, silent=True)
        if 'symbol' in result.columns:
            return result['symbol'].dropna().unique().tolist()
        return []
    
    def parse_gene_list(self, genes_input: str) -> List[str]:
        if os.path.isfile(genes_input):
            with open(genes_input) as f:
                genes = [line.strip() for line in f if line.strip()]
        else:
            genes = [g.strip() for g in genes_input.split(',')]
        return [g.upper() for g in genes if g]

    # =========================================
    # CLINVAR FETCHING
    # =========================================
    
    def fetch_clinvar(self, genes: List[str]) -> Dict[str, Dict]:
        """Fetch ClinVar data, returns dict keyed by rsID"""
        results = {}
        
        for gene in tqdm(genes, desc="Fetching ClinVar", disable=not self.verbose):
            variants = self._fetch_clinvar_gene(gene)
            for v in variants:
                rsid = v.get('rsid')
                if rsid and rsid.startswith('rs'):
                    # If already exists, keep the one with more info
                    if rsid not in results or self._has_more_info(v, results[rsid]):
                        results[rsid] = v
            time.sleep(0.1)
        
        self.log(f"ClinVar: {len(results)} unique rsIDs")
        return results
    
    def _fetch_clinvar_gene(self, gene: str) -> List[Dict]:
        url = f"{self.MYVARIANT_API}/query"
        params = {
            'q': f'clinvar.gene.symbol:{gene}',
            'fields': ','.join([
                'clinvar.rcv.clinical_significance',
                'clinvar.rcv.conditions',
                'dbsnp.rsid',
                'dbsnp.ref',
                'dbsnp.alt',
                'dbsnp.chrom',
                'dbsnp.hg19',
                'dbsnp.hg38',
                'cadd.phred',
                'cadd.1000g',
                'gnomad_genome.af',
                'gnomad_exome.af',
                'dbnsfp.sift.pred',
                'dbnsfp.sift.score',
                'dbnsfp.polyphen2.hdiv.pred',
                'dbnsfp.polyphen2.hdiv.score',
            ]),
            'size': 1000
        }
        
        try:
            resp = self.session.get(url, params=params, timeout=30)
            resp.raise_for_status()
            hits = resp.json().get('hits', [])
            
            results = []
            for h in hits:
                parsed = self._parse_clinvar_hit(h, gene)
                if parsed:
                    results.append(parsed)
            return results
        except Exception as e:
            self.log(f"  ClinVar error for {gene}: {e}")
            return []
    
    def _parse_clinvar_hit(self, h: dict, gene: str) -> Optional[Dict]:
        # Get rsID
        rsid = h.get('dbsnp', {}).get('rsid', '')
        if not rsid:
            _id = h.get('_id', '')
            if _id.startswith('rs'):
                rsid = _id
        
        if not rsid:
            return None
        
        # Clinical significance
        clin_sig = ''
        conditions = []
        rcv = h.get('clinvar', {}).get('rcv', [])
        if isinstance(rcv, dict):
            rcv = [rcv]
        if isinstance(rcv, list):
            for r in rcv:
                if isinstance(r, dict):
                    sig = r.get('clinical_significance', '')
                    if sig and not clin_sig:
                        clin_sig = sig
                    cond = r.get('conditions', {})
                    if isinstance(cond, dict) and cond.get('name'):
                        conditions.append(cond['name'])
                    elif isinstance(cond, list):
                        for c in cond:
                            if isinstance(c, dict) and c.get('name'):
                                conditions.append(c['name'])
        
        # dbSNP info
        dbsnp = h.get('dbsnp', {})
        chrom = dbsnp.get('chrom', '')
        ref = dbsnp.get('ref', '')
        alt = dbsnp.get('alt', '')
        if isinstance(alt, list):
            alt = ','.join(str(a) for a in alt)
        
        # Positions
        pos_hg19 = self._extract_pos(dbsnp.get('hg19'))
        pos_hg38 = self._extract_pos(dbsnp.get('hg38'))
        
        # MAF from various sources
        maf = self._get_maf(h)
        
        # SIFT/PolyPhen
        dbnsfp = h.get('dbnsfp', {})
        sift_pred = self._first_val(dbnsfp.get('sift', {}).get('pred'))
        sift_score = self._first_val(dbnsfp.get('sift', {}).get('score'))
        polyphen_pred = self._first_val(dbnsfp.get('polyphen2', {}).get('hdiv', {}).get('pred'))
        polyphen_score = self._first_val(dbnsfp.get('polyphen2', {}).get('hdiv', {}).get('score'))
        
        # CADD
        cadd = self._first_val(h.get('cadd', {}).get('phred'))
        
        # Population frequencies from 1000G
        cadd_1000g = h.get('cadd', {}).get('1000g', {})
        
        return {
            'rsid': rsid,
            'gene': gene,
            'chr': chrom,
            'position_hg19': pos_hg19,
            'position_hg38': pos_hg38,
            'ref_allele': ref,
            'alt_allele': alt,
            'clinical_significance': clin_sig,
            'conditions': '; '.join(conditions[:3]),
            'maf_global': maf,
            'maf_eur': self._first_val(cadd_1000g.get('eur_af')),
            'maf_afr': self._first_val(cadd_1000g.get('afr_af')),
            'maf_eas': self._first_val(cadd_1000g.get('eas_af')),
            'maf_amr': self._first_val(cadd_1000g.get('amr_af')),
            'maf_sas': self._first_val(cadd_1000g.get('sas_af')),
            'sift_pred': self._map_sift(sift_pred),
            'sift_score': sift_score,
            'polyphen_pred': self._map_polyphen(polyphen_pred),
            'polyphen_score': polyphen_score,
            'cadd_phred': cadd,
            'source': 'clinvar',
        }
    
    def _get_maf(self, h: dict) -> Optional[float]:
        """Extract MAF from multiple sources"""
        # Try gnomAD genome
        gnomad = h.get('gnomad_genome', {}).get('af', {})
        if isinstance(gnomad, dict):
            af = gnomad.get('af')
            if af is not None:
                return self._first_val(af)
        
        # Try gnomAD exome
        gnomad = h.get('gnomad_exome', {}).get('af', {})
        if isinstance(gnomad, dict):
            af = gnomad.get('af')
            if af is not None:
                return self._first_val(af)
        
        # Try 1000G from CADD
        cadd_1000g = h.get('cadd', {}).get('1000g', {})
        if isinstance(cadd_1000g, dict):
            af = cadd_1000g.get('af')
            if af is not None:
                return self._first_val(af)
        
        return None

    # =========================================
    # GWAS FETCHING
    # =========================================
    
    def fetch_gwas(self, genes: List[str]) -> Dict[str, Dict]:
        """Fetch GWAS data, returns dict keyed by rsID"""
        results = {}
        
        for gene in tqdm(genes, desc="Fetching GWAS", disable=not self.verbose):
            associations = self._fetch_gwas_gene(gene)
            for a in associations:
                rsid = a.get('rsid')
                if rsid and rsid.startswith('rs'):
                    # Keep best p-value
                    if rsid not in results:
                        results[rsid] = a
                    else:
                        if a.get('p_value', 1) < results[rsid].get('p_value', 1):
                            results[rsid] = a
            time.sleep(0.2)
        
        self.log(f"GWAS: {len(results)} unique rsIDs")
        return results
    
    def _fetch_gwas_gene(self, gene: str) -> List[Dict]:
        url = f"{self.GWAS_API}/singleNucleotidePolymorphisms/search/findByGene"
        
        try:
            resp = self.session.get(url, params={"geneName": gene}, timeout=30)
            if resp.status_code == 404:
                return []
            resp.raise_for_status()
            
            snps = resp.json().get("_embedded", {}).get("singleNucleotidePolymorphisms", [])
            
            results = []
            for snp in snps:
                rsid = snp.get("rsId", "")
                if not rsid:
                    continue
                
                # Get associations for this SNP
                assocs = self._fetch_snp_associations(rsid, gene)
                results.extend(assocs)
            
            return results
        except Exception as e:
            return []
    
    def _fetch_snp_associations(self, rsid: str, gene: str) -> List[Dict]:
        url = f"{self.GWAS_API}/singleNucleotidePolymorphisms/{rsid}/associations"
        
        try:
            resp = self.session.get(url, timeout=30)
            if resp.status_code == 404:
                return []
            resp.raise_for_status()
            
            assocs = resp.json().get("_embedded", {}).get("associations", [])
            
            results = []
            for a in assocs:
                parsed = self._parse_gwas_association(a, rsid, gene)
                if parsed:
                    results.append(parsed)
            
            return results
        except:
            return []
    
    def _parse_gwas_association(self, a: dict, rsid: str, gene: str) -> Optional[Dict]:
        # P-value
        p_value = None
        for pval in a.get("pvalues", []):
            try:
                mantissa = float(pval.get("mantissa", 1))
                exponent = int(pval.get("exponent", 0))
                p_value = mantissa * (10 ** exponent)
                break
            except:
                pass
        
        if p_value is None:
            try:
                p_value = float(a.get("pvalue", 1))
            except:
                p_value = 1.0
        
        # Filter by p-value
        if p_value > self.p_threshold:
            return None
        
        # OR/beta
        or_val = a.get("orPerCopyNum")
        beta = a.get("betaNum")
        ci_text = a.get("range", "")
        ci_lower, ci_upper = self._parse_ci(ci_text)
        
        # Risk allele
        risk_alleles = a.get("riskAlleles", [])
        risk_allele = ""
        risk_allele_freq = None
        if risk_alleles:
            ra = risk_alleles[0]
            risk_allele = ra.get("riskAlleleName", "").split("-")[-1]
            try:
                risk_allele_freq = float(ra.get("riskFrequency", ""))
            except:
                pass
        
        # Trait
        traits = []
        for t in a.get("efoTraits", []):
            trait_name = t.get("trait", "")
            if trait_name:
                traits.append(trait_name)
        
        # Study info
        study = a.get("study", {})
        pubmed_id = study.get("pubmedId", "")
        pub_date = study.get("publicationDate", "")
        study_year = pub_date[:4] if pub_date else ""
        sample_size = self._parse_sample_size(study.get("initialSampleSize", ""))
        
        return {
            'rsid': rsid,
            'gene': gene,
            'or': or_val,
            'beta': beta,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'p_value': p_value,
            'risk_allele': risk_allele,
            'risk_allele_freq': risk_allele_freq,
            'trait': '; '.join(traits[:3]),
            'pubmed_id': pubmed_id,
            'study_year': study_year,
            'sample_size': sample_size,
            'source': 'gwas',
        }
    
    def _parse_ci(self, ci_text: str) -> tuple:
        if not ci_text:
            return None, None
        try:
            match = re.search(r'\[?([\d.]+)\s*[-–]\s*([\d.]+)\]?', ci_text)
            if match:
                return float(match.group(1)), float(match.group(2))
        except:
            pass
        return None, None
    
    def _parse_sample_size(self, text: str) -> Optional[int]:
        if not text:
            return None
        try:
            numbers = re.findall(r'([\d,]+)\s*(?:individuals|cases|samples)', text, re.I)
            if numbers:
                return int(numbers[0].replace(',', ''))
            numbers = re.findall(r'[\d,]+', text)
            if numbers:
                return max(int(n.replace(',', '')) for n in numbers)
        except:
            pass
        return None

    # =========================================
    # ENRICHMENT (dbSNP for GWAS entries)
    # =========================================
    
    def enrich_with_dbsnp(self, gwas_data: Dict[str, Dict]) -> Dict[str, Dict]:
        """Add dbSNP/MyVariant info to GWAS entries"""
        rsids = list(gwas_data.keys())
        if not rsids:
            return gwas_data
        
        batch_size = 100
        for i in tqdm(range(0, len(rsids), batch_size), 
                      desc="Enriching GWAS with dbSNP", disable=not self.verbose):
            batch = rsids[i:i+batch_size]
            enriched = self._fetch_dbsnp_batch(batch)
            
            for rsid, info in enriched.items():
                if rsid in gwas_data:
                    gwas_data[rsid].update(info)
            
            time.sleep(0.1)
        
        return gwas_data
    
    def _fetch_dbsnp_batch(self, rsids: List[str]) -> Dict[str, Dict]:
        url = f"{self.MYVARIANT_API}/variant"
        fields = [
            'dbsnp.rsid', 'dbsnp.ref', 'dbsnp.alt', 'dbsnp.chrom',
            'dbsnp.hg19', 'dbsnp.hg38',
            'cadd.phred', 'cadd.1000g',
            'gnomad_genome.af', 'gnomad_exome.af',
            'dbnsfp.sift.pred', 'dbnsfp.sift.score',
            'dbnsfp.polyphen2.hdiv.pred', 'dbnsfp.polyphen2.hdiv.score',
        ]
        
        try:
            resp = self.session.post(url, json={"ids": rsids},
                                     params={"fields": ",".join(fields)},
                                     timeout=30)
            resp.raise_for_status()
            data = resp.json()
            
            results = {}
            if isinstance(data, list):
                for item in data:
                    if not isinstance(item, dict):
                        continue
                    rsid = item.get('dbsnp', {}).get('rsid', item.get('_id', ''))
                    if rsid:
                        results[rsid] = self._parse_dbsnp_info(item)
            
            return results
        except Exception as e:
            return {}
    
    def _parse_dbsnp_info(self, item: dict) -> Dict:
        dbsnp = item.get('dbsnp', {})
        dbnsfp = item.get('dbnsfp', {})
        
        ref = dbsnp.get('ref', '')
        alt = dbsnp.get('alt', '')
        if isinstance(alt, list):
            alt = ','.join(str(a) for a in alt)
        
        # MAF
        maf = self._get_maf(item)
        
        # 1000G frequencies
        cadd_1000g = item.get('cadd', {}).get('1000g', {})
        
        return {
            'chr': dbsnp.get('chrom', ''),
            'position_hg19': self._extract_pos(dbsnp.get('hg19')),
            'position_hg38': self._extract_pos(dbsnp.get('hg38')),
            'ref_allele': ref,
            'alt_allele': alt,
            'maf_global': maf,
            'maf_eur': self._first_val(cadd_1000g.get('eur_af')),
            'maf_afr': self._first_val(cadd_1000g.get('afr_af')),
            'maf_eas': self._first_val(cadd_1000g.get('eas_af')),
            'maf_amr': self._first_val(cadd_1000g.get('amr_af')),
            'maf_sas': self._first_val(cadd_1000g.get('sas_af')),
            'sift_pred': self._map_sift(self._first_val(dbnsfp.get('sift', {}).get('pred'))),
            'sift_score': self._first_val(dbnsfp.get('sift', {}).get('score')),
            'polyphen_pred': self._map_polyphen(self._first_val(dbnsfp.get('polyphen2', {}).get('hdiv', {}).get('pred'))),
            'polyphen_score': self._first_val(dbnsfp.get('polyphen2', {}).get('hdiv', {}).get('score')),
            'cadd_phred': self._first_val(item.get('cadd', {}).get('phred')),
        }

    # =========================================
    # MERGE & FILTER
    # =========================================
    
    def merge_data(self, clinvar: Dict[str, Dict], gwas: Dict[str, Dict]) -> List[Dict]:
        """Merge ClinVar and GWAS by rsID"""
        all_rsids = set(clinvar.keys()) | set(gwas.keys())
        self.log(f"Merging: {len(clinvar)} ClinVar + {len(gwas)} GWAS = {len(all_rsids)} unique rsIDs")
        
        merged = []
        for rsid in all_rsids:
            cv = clinvar.get(rsid, {})
            gw = gwas.get(rsid, {})
            
            if cv and gw:
                # Both sources - merge
                record = self._merge_records(cv, gw)
                record['source'] = 'both'
            elif cv:
                record = cv.copy()
                record['source'] = 'clinvar'
            else:
                record = gw.copy()
                record['source'] = 'gwas'
            
            merged.append(record)
        
        return merged
    
    def _merge_records(self, cv: Dict, gw: Dict) -> Dict:
        """Merge ClinVar and GWAS records for same rsID"""
        # Start with ClinVar (has more variant info)
        record = cv.copy()
        
        # Add GWAS fields
        gwas_fields = ['or', 'beta', 'ci_lower', 'ci_upper', 'p_value',
                       'risk_allele', 'risk_allele_freq', 'trait',
                       'pubmed_id', 'study_year', 'sample_size']
        
        for field in gwas_fields:
            if field in gw and gw[field]:
                record[field] = gw[field]
        
        # Fill missing position/allele info from GWAS
        for field in ['chr', 'position_hg19', 'position_hg38', 'ref_allele', 'alt_allele']:
            if not record.get(field) and gw.get(field):
                record[field] = gw[field]
        
        # Fill missing MAF
        if not record.get('maf_global') and gw.get('maf_global'):
            record['maf_global'] = gw['maf_global']
        
        return record
    
    def filter_results(self, data: List[Dict]) -> List[Dict]:
        """Apply MAF and clinical significance filters"""
        filtered = []
        
        for record in data:
            # MAF filter
            maf = record.get('maf_global')
            if maf is None:
                # No MAF info - check if it passes clinical significance filters
                if not self._passes_clinical_filter(record):
                    continue
            else:
                try:
                    if float(maf) < self.min_maf:
                        # Below MAF threshold - check clinical significance
                        if not self._passes_clinical_filter(record):
                            continue
                except:
                    if not self._passes_clinical_filter(record):
                        continue
            
            filtered.append(record)
        
        self.log(f"After filtering: {len(filtered)} variants")
        return filtered
    
    def _passes_clinical_filter(self, record: Dict) -> bool:
        """Check if record passes clinical significance filters"""
        # If from GWAS (has OR/p-value), keep it
        if record.get('or') or record.get('p_value'):
            return True
        
        # Check ClinVar significance
        clin_sig = str(record.get('clinical_significance', '')).lower()
        
        if self.include_pathogenic and 'pathogenic' in clin_sig:
            return True
        if self.include_risk_factor and 'risk factor' in clin_sig:
            return True
        if self.include_drug_response and 'drug response' in clin_sig:
            return True
        
        # If no clinical filters enabled, reject
        if not any([self.include_pathogenic, self.include_risk_factor, self.include_drug_response]):
            return False
        
        return False

    # =========================================
    # HELPERS
    # =========================================
    
    def _extract_pos(self, pos_data) -> str:
        if isinstance(pos_data, dict):
            return str(pos_data.get('start', ''))
        elif isinstance(pos_data, list) and pos_data:
            if isinstance(pos_data[0], dict):
                return str(pos_data[0].get('start', ''))
        return ''
    
    def _first_val(self, val):
        if val is None:
            return None
        if hasattr(val, '__iter__') and not isinstance(val, (str, dict)):
            val = list(val)
            return val[0] if val else None
        return val
    
    def _map_sift(self, pred) -> str:
        if not pred:
            return ''
        mapping = {'D': 'Deleterious', 'T': 'Tolerated'}
        return mapping.get(str(pred), str(pred))
    
    def _map_polyphen(self, pred) -> str:
        if not pred:
            return ''
        mapping = {'D': 'Probably damaging', 'P': 'Possibly damaging', 'B': 'Benign'}
        return mapping.get(str(pred), str(pred))
    
    def _has_more_info(self, new: Dict, old: Dict) -> bool:
        """Check if new record has more info than old"""
        new_score = sum(1 for v in new.values() if v)
        old_score = sum(1 for v in old.values() if v)
        return new_score > old_score

    # =========================================
    # MAIN PIPELINE
    # =========================================
    
    def run(self, genes: List[str]) -> pd.DataFrame:
        self.log(f"Processing {len(genes)} genes...")
        self.log(f"Filters: MAF ≥ {self.min_maf}, p-value ≤ {self.p_threshold}")
        
        # Fetch ClinVar
        clinvar_data = self.fetch_clinvar(genes)
        
        # Fetch GWAS
        gwas_data = self.fetch_gwas(genes)
        
        # Enrich GWAS with dbSNP info
        gwas_data = self.enrich_with_dbsnp(gwas_data)
        
        # Merge
        merged = self.merge_data(clinvar_data, gwas_data)
        
        # Filter
        filtered = self.filter_results(merged)
        
        if not filtered:
            return pd.DataFrame()
        
        # Build DataFrame
        columns = [
            'gene', 'rsid', 'source', 'chr', 'position_hg19', 'position_hg38',
            'ref_allele', 'alt_allele', 'risk_allele',
            'maf_global', 'maf_eur', 'maf_afr', 'maf_eas', 'maf_amr', 'maf_sas',
            'or', 'beta', 'ci_lower', 'ci_upper', 'p_value',
            'clinical_significance', 'conditions', 'trait',
            'sift_pred', 'sift_score', 'polyphen_pred', 'polyphen_score', 'cadd_phred',
            'pubmed_id', 'study_year', 'sample_size',
        ]
        
        df = pd.DataFrame(filtered)
        
        # Reorder columns
        existing_cols = [c for c in columns if c in df.columns]
        other_cols = [c for c in df.columns if c not in columns]
        df = df[existing_cols + other_cols]
        
        # Format numeric columns
        for col in ['maf_global', 'maf_eur', 'maf_afr', 'maf_eas', 'maf_amr', 'maf_sas', 
                    'sift_score', 'polyphen_score']:
            if col in df.columns:
                df[col] = df[col].apply(lambda x: f"{float(x):.4f}" if pd.notna(x) and x != '' else '')
        
        for col in ['or', 'beta', 'ci_lower', 'ci_upper', 'cadd_phred']:
            if col in df.columns:
                df[col] = df[col].apply(lambda x: f"{float(x):.3f}" if pd.notna(x) and x != '' else '')
        
        if 'p_value' in df.columns:
            df['p_value'] = df['p_value'].apply(lambda x: f"{float(x):.2e}" if pd.notna(x) and x != '' else '')
        
        # Sort by source (both first), then p-value
        source_order = {'both': 0, 'gwas': 1, 'clinvar': 2}
        df['_sort'] = df['source'].map(source_order)
        df = df.sort_values(['_sort', 'gene', 'rsid']).drop(columns=['_sort'])
        
        return df
    
    def save(self, df: pd.DataFrame, path: str):
        if path.endswith('.xlsx'):
            df.to_excel(path, index=False, freeze_panes=(1, 0))
        elif path.endswith('.tsv'):
            df.to_csv(path, sep='\t', index=False)
        else:
            df.to_csv(path, index=False)
        self.log(f"Saved {len(df)} variants to {path}")


def main():
    parser = argparse.ArgumentParser(
        description='Universal SNP Extractor: ClinVar + GWAS combined',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --genes TYR,MC1R -o results.csv
  %(prog)s --kegg pathway.xml --min-maf 0.01 -o results.xlsx
  %(prog)s --genes APOE --pathogenic --risk-factor -o apoe.csv
  %(prog)s --genes BRCA1 --min-maf 0 --pathogenic -o rare_pathogenic.csv
        """
    )
    
    # Input
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--genes', '-g', help='Comma-separated gene symbols')
    input_group.add_argument('--gene-file', '-f', help='File with genes (one per line)')
    input_group.add_argument('--kegg', '-k', help='KEGG pathway XML')
    input_group.add_argument('--kegg-folder', '-K', help='Folder with KEGG XMLs')
    
    # Filters
    parser.add_argument('--min-maf', type=float, default=0.05,
                        help='Minimum MAF (default: 0.05)')
    parser.add_argument('--p-threshold', type=float, default=5e-8,
                        help='GWAS p-value threshold (default: 5e-8)')
    parser.add_argument('--pathogenic', action='store_true',
                        help='Include ClinVar pathogenic (even if MAF < threshold)')
    parser.add_argument('--risk-factor', action='store_true',
                        help='Include ClinVar risk factors')
    parser.add_argument('--drug-response', action='store_true',
                        help='Include ClinVar drug response')
    
    # Output
    parser.add_argument('--output', '-o', default='snp_results.csv',
                        help='Output file (.csv, .tsv, .xlsx)')
    parser.add_argument('--quiet', '-q', action='store_true')
    
    args = parser.parse_args()
    
    extractor = UniversalSNPExtractor(
        min_maf=args.min_maf,
        p_threshold=args.p_threshold,
        include_pathogenic=args.pathogenic,
        include_risk_factor=args.risk_factor,
        include_drug_response=args.drug_response,
        verbose=not args.quiet,
    )
    
    # Get genes
    if args.genes:
        genes = extractor.parse_gene_list(args.genes)
    elif args.gene_file:
        genes = extractor.parse_gene_list(args.gene_file)
    elif args.kegg:
        genes = extractor.parse_kegg_xml(args.kegg)
    elif args.kegg_folder:
        genes = extractor.parse_kegg_folder(args.kegg_folder)
    
    if not genes:
        print("No genes found!")
        sys.exit(1)
    
    # Run
    results = extractor.run(genes)
    
    if results.empty:
        print("No variants found!")
        sys.exit(0)
    
    # Save
    extractor.save(results, args.output)
    
    # Summary
    print(f"\n{'='*50}")
    print("SUMMARY")
    print(f"{'='*50}")
    print(f"Genes: {len(genes)}")
    print(f"Variants: {len(results)}")
    
    if 'source' in results.columns:
        print(f"\nBy source:")
        for src, cnt in results['source'].value_counts().items():
            print(f"  {src}: {cnt}")
    
    print(f"\nOutput: {args.output}")


if __name__ == '__main__':
    main()