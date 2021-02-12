from tqdm.auto import tqdm
import json
import pickle
from collections import defaultdict

import zipfile
import io
import gzip
import re

import chembl_structure_pipeline as csp  # https://github.com/chembl/ChEMBL_Structure_Pipeline.git
import rdkit

class GsrsIndex():
    """
    Data downloaded from https://gsrs.ncats.nih.gov/#/
    """
    
    def __init__(self, data_dir='.'):
        rdkit.RDLogger.DisableLog('rdApp.*')
        
        self.data_dir = data_dir
        
        try:
            self.load_indexes()
        except:
            self.gsrs_dict = None
            self.gsrs_inchis = None
        
    def fetch_data(self, fn, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
        
        gsrs_records = []
        with open(f'{data_dir}/{fn}', 'rb') as f:
            gsrs_zipfile = zipfile.ZipFile(f, 'r')
            for fp in gsrs_zipfile.filelist:
                if fp.filename.split('/')[-1] == 'smallSeedData.gsrs':
                    gsrs_datafile = gsrs_zipfile.read(fp)

            gsrs_data = gzip.open(io.BytesIO(gsrs_datafile)).read().decode('latin-1')

            for line in tqdm(re.split('\n', gsrs_data)):
                if line:
                    line = '{' + '{'.join(line.split('{')[1:])
                    gsrs_records.append(json.loads(line))
        
        self.gsrs_dict = {}
        for gr in gsrs_records:
            if gr['substanceClass'] in {'concept'}:
                continue
            if not gr['status'] in {'approved'}:
                continue

            unii = gr['approvalID']
            substance_class = gr['substanceClass']
            definition_level = gr['definitionLevel']
            names = [(n['name'], n['type']) for n in gr['names']]
            codes = defaultdict(set)
            for n in gr['codes']:
                try:
                    codes[n['codeSystem']].add((n['type'], n['code']))
                except KeyError as e:
                    pass
                    
            try: 
                structure = {
                    'molfile': gr['structure']['molfile'],
                    'atropisomerism': gr['structure']['atropisomerism'],
                    'stereoCenters': gr['structure']['stereoCenters'],
                    'definedStereo': gr['structure']['definedStereo'],
                    'ezCenters': gr['structure']['ezCenters'],
                    'charge': gr['structure']['charge'],
                    'stereochemistry': gr['structure']['stereochemistry']
                }
            except: 
                structure = {}

            relationships = []
            for r in gr['relationships']:
                if 'approvalID' in r['relatedSubstance']:
                    relationships.append((r['relatedSubstance']['approvalID'], r['type']))

            self.gsrs_dict[unii] = {
                'substance_class': substance_class, 
                'definition_level': definition_level, 
                'names': names, 
                'codes': {k:list(vs) for k,vs in codes.items()}, 
                'structure': structure, 
                'relationships': relationships
            }
    
    def gen_inchi_index(self):
        def molblock2inchi(molblock):
            try:
                mol = rdkit.Chem.MolFromMolBlock(molblock)
                return rdkit.Chem.MolToInchi(mol)
            except:
                return None
    
        gsrs_mol = {}
        no_struct = []
        for k,v in tqdm(self.gsrs_dict.items(), leave=True, position=0, desc='Standardising MolFiles'):
            issues = None
            if 'molfile' in v['structure']:
                molblock = v['structure']['molfile']
            else:
                no_struct.append(k)
                continue

            issues = csp.checker.check_molblock(molblock)

            standardised_molblock = csp.standardizer.standardize_molblock(molblock)  # standardise
            parent_molblock = csp.standardizer.get_parent_molblock(molblock)  # get parent

            gsrs_mol[k] = {
                'raw': molblock,
                'standardised': standardised_molblock,
                'parent': parent_molblock,
                'issues': issues
            }
        
        self.gsrs_inchis = {}
        for unii,d in tqdm(gsrs_mol.items(), leave=True, position=0, desc='Converting to InChI'):
            # filter the issues out
            if max([0]+[i[0] for i in d['issues']]) > 2: 
                continue

            # generate InChiKeys
            self.gsrs_inchis[unii] = {k:molblock2inchi(v) for k,v in d.items()}

        for unii,d in self.gsrs_inchis.items():
             self.gsrs_inchis[unii] = {k:v for k,v in d.items() if (k in {'raw', 'standardised', 'parent'}) and v}
    
    def save_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f"{data_dir}/gsrs_dict.json", 'wt') as f:
            json.dump(self.gsrs_dict, f)
        with open(f"{data_dir}/gsrs_inchis.json", 'wt') as f:
            json.dump(self.gsrs_inchis, f)
        
    def load_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f"{data_dir}/gsrs_dict.json", 'rt') as f:
            self.gsrs_dict = json.load(f)
        with open(f"{data_dir}/gsrs_inchis.json", 'rt') as f:
            self.gsrs_inchis = json.load(f)
            
    def query(self, unii):
        return self.gsrs_dict[unii]
    
    def get_inchi(self, unii):
        if unii in self.gsrs_inchis:
            return self.gsrs_inchis[unii]
        
