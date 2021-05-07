import pickle
import json
import sqlalchemy as sa
import cx_Oracle
from inchicompare import inchicompare as ic  # https://github.com/timrozday/inchicompare.git
from collections import defaultdict
import rdkit
from tqdm.auto import tqdm

import chembl_ident

class ChemblStructureIndex():
    def __init__(self, data_dir='.', chembl_index=None):
        self.data_dir = data_dir
        self.chembl_db = None
        if chembl_index is None:
            self.chembl_index = chembl_ident.ChemblIndexes(data_dir=self.data_dir)
        else:
            self.chembl_index = chembl_index
        
        rdkit.RDLogger.DisableLog('rdApp.*')
        
        try:
            self.load_indexes()
        except:
            pass
        
    def connect_to_chembl(self, username, password, url):
        self.chembl_db = cx_Oracle.connect(username, password, url)
#         print(chembl_db.dsn, "running version", chembl_db.version)
        
        return self.chembl_db

    def load_inactive_compounds(self, data_dir=None):
                
        if data_dir is None:
            data_dir = self.data_dir
            
        # load inactive data
        with open(f"{data_dir}/exclude_inchis.json", 'rt') as f:
            self.exclude_inchis = set(json.load(f))
        with open(f"{data_dir}/split_inactive_inchis.json", 'rt') as f:
            self.split_inactive_inchis = set(json.load(f))
            
        with open(f"{data_dir}/salts.smi", 'rt') as f:
            cols = ['name', 'smiles'] 
            salts = [{cols[i]:v for i,v in enumerate(l.split('\t'))} for l in f]

        with open(f"{data_dir}/solvents.smi", 'rt') as f:
            cols = ['name', 'smiles'] 
            solvents = [{cols[i]:v for i,v in enumerate(l.split('\t'))} for l in f]
        
        salts_inchi = set()
        for d in salts:
            try:
                mol = rdkit.Chem.MolFromSmiles(d['smiles'])
                inchi = rdkit.Chem.MolToInchi(mol)
                salts_inchi.add(inchi)
            except Exception as e:
                pass
        
        solvents_inchi = set()
        for d in solvents:
            try:
                mol = rdkit.Chem.MolFromSmiles(d['smiles'])
                inchi = rdkit.Chem.MolToInchi(mol)
                solvents_inchi.add(inchi)
            except Exception as e:
                pass
        
        self.split_inactive_inchis.update(salts_inchi)
        self.split_inactive_inchis.update(solvents_inchi)
        
        self.conn_split_inactive_inchis = {ic.inchi_conn_layer(i) for i in self.split_inactive_inchis}
        
        
    def fetch_data(self, data_dir=None):
        
        if data_dir is None:
            data_dir = self.data_dir
        
        self.load_inactive_compounds(data_dir=data_dir)
        
        chembl_cursor = self.chembl_db.cursor()
        
        # fetch structures
        
        self.compound_inchis = {}
        for molregno, inchi in tqdm(chembl_cursor.execute("select MOLREGNO, STANDARD_INCHI from CHEMBL.COMPOUND_STRUCTURES"), leave=True, position=0, desc='ChEMBL structures'):
            if inchi is None:
                continue
            
            ci = self.chembl_index.get_chembl_ident(molregno=molregno)
            self.compound_inchis[ci] = inchi
        
        structure_id2drugbase_id = {}
        for drugbase_id, mrn, structure_id in tqdm(chembl_cursor.execute('select ID, MOLREGNO, MOLECULE_STRUCTURE_ID from DRUGBASE.MOLECULE_DICTIONARY'), leave=True, position=0, desc='Drugbase structure IDs'):
            ci = self.chembl_index.get_chembl_ident(drugbase_id=drugbase_id)
            structure_id2drugbase_id[structure_id] = ci
        
        for structure_id, inchi in tqdm(chembl_cursor.execute("select MOLECULE_STRUCTURE_ID, INCHI from DRUGBASE.MOLECULE_STRUCTURE"), leave=True, position=0, desc='Drugbase structures'):
            if inchi is None:
                continue
            inchi = inchi.read()
            if inchi is None:
                continue
            if structure_id in structure_id2drugbase_id:
                ci = structure_id2drugbase_id[structure_id]
                self.compound_inchis[ci] = inchi
        
        self.inchi_index = defaultdict(set)
        for k,v in self.compound_inchis.items():
            self.inchi_index[v].add(k)
        
        # split InChIs describing multiple molecules
        
        self.split_inchi_index = defaultdict(set)
        for ci, inchi in tqdm(self.compound_inchis.items(), leave=True, position=0, desc='Split InChIs'):    
            mol = rdkit.Chem.MolFromInchi(inchi)
            if mol is None:
                continue
            
            try:
                for m in rdkit.Chem.rdmolops.GetMolFrags(mol, asMols=True):  # split
                    i = rdkit.Chem.MolToInchi(m)
                    c = ic.inchi_conn_layer(i)
                    if not c in self.conn_split_inactive_inchis:
                        self.split_inchi_index[i].add(ci)  # add to index
            except:
                pass
            
        self.split_inchi_index = dict(self.split_inchi_index)   
        
        # extract connectivity layer of InChI and build index
        
        self.inchi_connectivity_index = defaultdict(set)
        for ci, inchi in tqdm(self.compound_inchis.items(), leave=True, position=0, desc='Connectivity layer InChI index'):
            try:
                inchi = ic.strip_inchi(inchi, exclude_inchis=self.conn_split_inactive_inchis)  # filter inchi
            except:
                pass

            if inchi:
                c = ic.inchi_conn_layer(inchi)
                self.inchi_connectivity_index[c].add(ci)
        self.inchi_connectivity_index = dict(self.inchi_connectivity_index)
        
        self.inchi_split_connectivity_index = defaultdict(set)
        for inchi, mols in tqdm(self.split_inchi_index.items(), leave=True, position=0, desc='Split connectivity layer InChI index'):
            c = ic.inchi_conn_layer(inchi)
        #     if not c in conn_split_inactive_inchis:
            self.inchi_split_connectivity_index[c].update(mols)
        
        self.inchi_split_connectivity_index = {k:v for k,v in self.inchi_split_connectivity_index.items() if not k in self.conn_split_inactive_inchis}  # remove innactive InChIs
        
        chembl_cursor.close()
        
    def save_indexes(self, data_dir=None):
        
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f"{data_dir}/compound_inchis.pkl", 'wb') as f:
            pickle.dump({k.__tuple__():vs for k,vs in self.compound_inchis.items()}, f)
        with open(f"{data_dir}/inchi_index.pkl", 'wb') as f:
            pickle.dump({k:{v.__tuple__() for v in vs} for k,vs in self.inchi_index.items()}, f)
        with open(f"{data_dir}/split_inchi_index.pkl", 'wb') as f:
            pickle.dump({k:{v.__tuple__() for v in vs} for k,vs in self.split_inchi_index.items()}, f)
        with open(f"{data_dir}/inchi_connectivity_index.pkl", 'wb') as f:
            pickle.dump({k:{v.__tuple__() for v in vs} for k,vs in self.inchi_connectivity_index.items()}, f)
        with open(f"{data_dir}/inchi_split_connectivity_index.pkl", 'wb') as f:
            pickle.dump({k:{v.__tuple__() for v in vs} for k,vs in self.inchi_split_connectivity_index.items()}, f)
    
    def load_indexes(self, data_dir=None):
        
        if data_dir is None:
            data_dir = self.data_dir
        
        self.load_inactive_compounds(data_dir=data_dir)
        
        with open(f"{data_dir}/compound_inchis.pkl", 'rb') as f:
            self.compound_inchis = {chembl_ident.ChemblIdent(*k):vs for k,vs in pickle.load(f).items()}
        with open(f"{data_dir}/inchi_index.pkl", 'rb') as f:
            self.inchi_index = {k:{chembl_ident.ChemblIdent(*v) for v in vs} for k,vs in pickle.load(f).items()}
        with open(f"{data_dir}/split_inchi_index.pkl", 'rb') as f:
            self.split_inchi_index = {k:{chembl_ident.ChemblIdent(*v) for v in vs} for k,vs in pickle.load(f).items()}
        with open(f"{data_dir}/inchi_connectivity_index.pkl", 'rb') as f:
            self.inchi_connectivity_index = {k:{chembl_ident.ChemblIdent(*v) for v in vs} for k,vs in pickle.load(f).items()}
        with open(f"{data_dir}/inchi_split_connectivity_index.pkl", 'rb') as f:
            self.inchi_split_connectivity_index = {k:{chembl_ident.ChemblIdent(*v) for v in vs} for k,vs in pickle.load(f).items()}
    

    def get_structure(self, obj=None, drugbase_id=None, molregno=None, chembl_id=None):
        if obj is None:
            obj = self.chembl_index.get_chembl_ident(drugbase_id=drugbase_id, molregno=molregno, chembl_id=chembl_id)
        if obj:
            if obj in self.compound_inchis:
                return self.compound_inchis[obj]
    
    def query_inchi_conn(self, inchi, strip=True):
        if strip:
            try:
                inchi = ic.strip_inchi(inchi, exclude_inchis=self.conn_split_inactive_inchis)  # filter inchi
            except:
                pass
        
        inchi_conn = ic.inchi_conn_layer(inchi)
        
        r = set()
        if inchi_conn in self.inchi_connectivity_index:
            r.update(self.inchi_connectivity_index[inchi_conn])
        if inchi_conn in self.inchi_split_connectivity_index:
            r.update(self.inchi_split_connectivity_index[inchi_conn])
        return r
            
    def query_inchi(self, inchi, strip=True):
        if strip:
            try:
                inchi = ic.strip_inchi(inchi, exclude_inchis=self.conn_split_inactive_inchis)  # filter inchi
            except:
                pass
        
        r = set()
        if inchi in self.inchi_index:
            r.update(self.inchi_index[inchi])
        if inchi in self.split_inchi_index:
            r.update(self.split_inchi_index[inchi])
        return r
    
    def is_active(self, inchi):
        stripped_inchi = ic.strip_inchi(inchi, exclude_inchis=self.conn_split_inactive_inchis)  # filter inchi
        return bool(stripped_inchi)
    
    def query(self, inchi, connectivity=True, strip=True, consistency=True, split=True, inactive=False, filter_layers={'h','f','p','q','i','t','b','m','s'}):
        if strip:
            try:
                stripped_inchi = ic.strip_inchi(inchi, exclude_inchis=self.conn_split_inactive_inchis)  # filter inchi
                assert bool(stripped_inchi)
                inchi = stripped_inchi
            except:
                if inactive == False:
                    return None
                else:
                    pass
        
        if split:
            mol = rdkit.Chem.MolFromInchi(inchi)
            if not mol is None:
                try:
                    results = None
                    for m in rdkit.Chem.rdmolops.GetMolFrags(mol, asMols=True):  # split
                        i = rdkit.Chem.MolToInchi(m)
                        r = self.query(i, connectivity=connectivity, strip=False, consistency=consistency, split=False)
                        if results:
                            results.update(r)
                        else:
                            results = r
                    return results
                except:
                    pass
        
        if connectivity:
            r = self.query_inchi_conn(inchi, strip=False)
        else:
            r = self.query_inchi(inchi, strip=False)
            
            
        if consistency:
            consistencies = {}
            for i in r:
                r_inchi = self.get_structure(i)  # get inchi
#                 consistencies[i] = ic.compare_consistent(inchi, r_inchi, filter_layers=filter_layers)  # compare with query to get consistency
                consistencies[i] = self.inchi_overlap(inchi, r_inchi, strip=strip, consistency=consistency, filter_layers=filter_layers)
            return consistencies
        else:
            return r
    
    def inchi_overlap(self, inchi1, inchi2, strip=True, consistency=True, filter_layers={'h','f','p','q','i','t','b','m','s'}):
        def split_strip(inchi, strip=True):
            mol = rdkit.Chem.MolFromInchi(inchi)
            split_inchis = set()
            for m in rdkit.Chem.rdmolops.GetMolFrags(mol, asMols=True):  # split
                i = rdkit.Chem.MolToInchi(m)

                if strip:
                    c = ic.inchi_conn_layer(i)
                    if not c in self.conn_split_inactive_inchis:
                        split_inchis.add(i)
                else:
                    split_inchis.add(i)

            return split_inchis

        if isinstance(inchi1, list):
            inchi1 = set(inchi1)
        if isinstance(inchi2, list):
            inchi2 = set(inchi2)

        if isinstance(inchi1, set):
            split_inchis1 = set(inchi1)
        else:
            split_inchis1 = split_strip(inchi1, strip=strip)

        if isinstance(inchi2, set):
            split_inchis2 = set(inchi2)
        else:
            split_inchis2 = split_strip(inchi2, strip=strip)

        # get conn layer
        conn_split_inchis1 = defaultdict(set)
        for i in split_inchis1:
            conn_split_inchis1[ic.inchi_conn_layer(i)].add(i)
        conn_split_inchis2 = defaultdict(set)
        for i in split_inchis2:
            conn_split_inchis2[ic.inchi_conn_layer(i)].add(i)

        # compare conn layers to get candidate matches
        candidates = set()
        for k in set(conn_split_inchis1.keys()) & set(conn_split_inchis2.keys()):
            candidates.update(it.product(conn_split_inchis1[k], conn_split_inchis2[k]))

        # check these candidate matches by chacking consistency
        matches = set()
        for i1,i2 in candidates:
            if consistency:
                c,s = ic.compare_consistent(i1, i2, filter_layers=filter_layers)
                if c:
                    matches.add((i1,i2))
            else:
                if i1==i2:
                    matches.add((i1,i2))

        # return overlap between the split inchi fragments
        i1_matches = {i1 for (i1,i2) in matches}
        i2_matches = {i2 for (i1,i2) in matches}
        overlap = i1_matches | i2_matches

        return split_inchis1-overlap, matches, split_inchis2-overlap  # remainder from inchi1, paired matches, remainder from inchi2
