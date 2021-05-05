# import json
import pickle
import re
from collections import defaultdict
from tqdm.auto import tqdm
import cx_Oracle
import sqlalchemy as sa
import sqlite3

import chembl_ident
from .chembl_name_db import ChemblNameDB


class ChemblNameIndex():
    def __init__(self, data_dir='.', chembl_index=None):
        self.data_dir = data_dir
        
        self.chembl_db = None
        self.name_store = None
        self.name_session = None
        
        if chembl_index is None:
            self.chembl_index = chembl_ident.ChemblIndexes(data_dir=self.data_dir)
        else:
            self.chembl_index = chembl_index
        
        try:
            self.load_query_index()
        except:
            self.name2substances = None
            self.substance2names = None
            self.filtered_name2substances = None
    
    def connect_to_chembl(self, username, password, url):
        self.chembl_db = cx_Oracle.connect(username, password, url)
#         print(chembl_db.dsn, "running version", chembl_db.version)
        
        return self.chembl_db
    
    def connect_to_namestore(self, create=False):
        self.name_store = ChemblNameDB(f"sqlite:///{self.data_dir}/chembl_names.sqlite")
        if create:
             self.name_store.create()
        
        self.name_session = self.name_store.SessionMaker()
        
        return self.name_store

    def fetch_chembl_data(self):
        if self.chembl_db is None:
            self.connect_to_chembl()
            
        chembl_cursor = self.chembl_db.cursor()
        
        sql_query = (   
                        "select MD.ID, MD.MOLREGNO, MD.PREF_NAME "
                        "from DRUGBASE.MOLECULE_DICTIONARY MD "
                        "where MD.DELETED = 0 "
                    )
        chembl_cursor.execute(sql_query)
        
        self.pref_names = []
        for drugbase_id, molregno, pref_name in tqdm(chembl_cursor, desc='Drugbase preferred names'):
            if pref_name == "NA":
                continue
            self.pref_names.append((drugbase_id, molregno, pref_name))

        sql_query = (   
                        "select MD.ID, MD.MOLREGNO, MS.NAME, MST.NAME "
                        "from DRUGBASE.MOLECULE_SYNONYM MS "
                        "left join DRUGBASE.MOLECULE_SYNONYM_TYPE MST "
                        "on MS.MOLECULE_SYNONYM_TYPE_ID = MST.ID "
                        "left join DRUGBASE.MOLECULE_DICTIONARY MD "
                        "on MS.MOLECULE_DICTIONARY_ID = MD.ID "
                        "where MD.DELETED = 0 "
                    )

        chembl_cursor.execute(sql_query)
        
        self.mol_syns = []
        for drugbase_id, molregno, name, name_type in tqdm(chembl_cursor, desc='Drugbase molecule synonyms'):
            self.mol_syns.append((drugbase_id, molregno, name, name_type))

        sql_query = (   
                        "select MD.ID, MD.MOLREGNO, MCN.NAME "
                        "from DRUGBASE.MOLECULE_CHEMICAL_NAME MCN "
                        "left join DRUGBASE.MOLECULE_DICTIONARY MD "
                        "on MCN.MOLECULE_DICTIONARY_ID = MD.ID "
                        "where MD.DELETED = 0 "
                    )
        chembl_cursor.execute(sql_query)
        
        self.chem_names = []
        for drugbase_id, molregno, name in tqdm(chembl_cursor, desc='Drugbase molecule chemical names'):
            if name:
                name = name.read()
            self.chem_names.append((drugbase_id, molregno, name))

        sql_query = (   
                        "select DC.MOLECULE_DICTIONARY_ID, DC.MOLREGNO, DC.DAILYMED_INGREDIENT "
                        "from DRUGBASE.DAILYMED_COMPOUNDS DC "
                    )

        chembl_cursor.execute(sql_query)
        
        self.dm_names = []
        for drugbase_id,molregno,dm_ingr in tqdm(chembl_cursor, desc='Drugbase DailyMed ingredient names'):
            self.dm_names.append((drugbase_id,molregno,dm_ingr))

        sql_query = (   
                        "select MD.CHEMBL_ID, MD.MOLREGNO, MD.PREF_NAME "
                        "from CHEMBL.MOLECULE_DICTIONARY MD "
                    )
        chembl_cursor.execute(sql_query)
        
        self.chembl_pref_names = []
        for chembl_id,mrn,pref_name in tqdm(chembl_cursor, desc='ChEMBL preferred names'):
            if pref_name == "NA":
                continue
            self.chembl_pref_names.append((chembl_id,mrn,pref_name))

        sql_query = (   
                        "select MD.CHEMBL_ID, MS.MOLREGNO, MS.SYNONYMS "
                        "from CHEMBL.MOLECULE_SYNONYMS MS "
                        "left join CHEMBL.MOLECULE_DICTIONARY MD "
                        "on MS.MOLREGNO = MD.MOLREGNO"
                    )
        chembl_cursor.execute(sql_query)
        
        self.chembl_mol_syns = []
        for chembl_id,mrn,synonyms in tqdm(chembl_cursor, desc='ChEMBL molecule synonyms'):
            self.chembl_mol_syns.append((chembl_id,mrn,synonyms))

        sql_query = (   
                        "select MD.CHEMBL_ID, CR.MOLREGNO, CR.COMPOUND_NAME "
                        "from CHEMBL.COMPOUND_RECORDS CR "
                        "left join CHEMBL.MOLECULE_DICTIONARY MD "
                        "on CR.MOLREGNO = MD.MOLREGNO"
                    )
        chembl_cursor.execute(sql_query)
        
        self.chembl_compound_names = []
        for chembl_id,mrn,compound_name in tqdm(chembl_cursor, desc='ChEMBL compound names'):
            self.chembl_compound_names.append((chembl_id,mrn,compound_name))

        sql_query = (   
                        "select MD.CHEMBL_ID, FO.MOLREGNO, PR.TRADE_NAME "
                        "from CHEMBL.FORMULATIONS FO left join CHEMBL.PRODUCTS PR "
                        "on FO.PRODUCT_ID = PR.PRODUCT_ID "
                        "left join CHEMBL.MOLECULE_DICTIONARY MD "
                        "on FO.MOLREGNO = MD.MOLREGNO"
                    )
        chembl_cursor.execute(sql_query)
        
        self.chembl_trade_names = []
        for chembl_id,mrn,trade_name in tqdm(chembl_cursor, desc='ChEMBL trade names'):
            self.chembl_trade_names.append((chembl_id,mrn,trade_name))
            
        chembl_cursor.close()
            
    def save_to_db(self, batch_size=1000):
        def batch(iterable, n=1):
            l = len(iterable)
            for ndx in range(0, l, n):
                yield iterable[ndx:min(ndx + n, l)]
                
        if self.name_store is None:
            self.connect_to_namestore()
        
        name_session = self.name_store.SessionMaker()
        
        # import all substances
        
        def update_substance_index(substance_index, g):
            for db_id,mrn,chembl_id in g:
                ci = self.chembl_index.get_chembl_ident(drugbase_id=db_id, molregno=mrn)
                if not ci in substance_index:
                    substance_index[ci] = len(substance_index)
            return substance_index
        
        substance_index = {self.chembl_index.get_chembl_ident(drugbase_id=db_id, molregno=mrn, chembl_id=None) for db_id,mrn,name in self.pref_names}
        substance_index |= {self.chembl_index.get_chembl_ident(drugbase_id=db_id, molregno=mrn, chembl_id=None) for db_id,mrn,name,name_type in self.mol_syns}
        substance_index |= {self.chembl_index.get_chembl_ident(drugbase_id=db_id, molregno=mrn, chembl_id=None) for db_id,mrn,name in self.chem_names}
        substance_index |= {self.chembl_index.get_chembl_ident(drugbase_id=db_id, molregno=mrn, chembl_id=None) for db_id,mrn,name in self.dm_names}
        substance_index |= {self.chembl_index.get_chembl_ident(drugbase_id=None, molregno=mrn, chembl_id=chembl_id) for chembl_id,mrn,name in self.chembl_pref_names}
        substance_index |= {self.chembl_index.get_chembl_ident(drugbase_id=None, molregno=mrn, chembl_id=chembl_id) for chembl_id,mrn,name in self.chembl_mol_syns}
        substance_index |= {self.chembl_index.get_chembl_ident(drugbase_id=None, molregno=mrn, chembl_id=chembl_id) for chembl_id,mrn,name in self.chembl_compound_names}
        substance_index |= {self.chembl_index.get_chembl_ident(drugbase_id=None, molregno=mrn, chembl_id=chembl_id) for chembl_id,mrn,name in self.chembl_trade_names}
        substance_index -= {None}
        substance_index = {ci:i for i,ci in enumerate(substance_index)}
        
        import_data = [{'id': i, 'molregno': ci.molregno, 'drugbase_id': ci.drugbase_id, 'chembl_id': ci.chembl_id} for ci,i in substance_index.items()]
        for b in tqdm(batch(import_data, batch_size), total=int(len(import_data)/batch_size)+1, position=0, leave=True, desc='Substances'):
            name_session.execute(self.name_store.Substance.__table__.insert(), b)
            name_session.flush()
        
        
        name_index = {(name,"drugbase pref",None) for db_id,mrn,name in self.pref_names if name}
        name_index |= {(name,"drugbase syn",name_type) for db_id,mrn,name,name_type in self.mol_syns if name}
        name_index |= {(name,"drugbase chem",None) for db_id,mrn,name in self.chem_names if name}
        name_index |= {(name,"drugbase dm",None) for db_id,mrn,name in self.dm_names if name}
        name_index |= {(name,"chembl pref",None) for chembl_id,mrn,name in self.chembl_pref_names if name}
        name_index |= {(name,"chembl syn",None) for chembl_id,mrn,name in self.chembl_mol_syns if name}
        name_index |= {(name,"chembl compound",None) for chembl_id,mrn,name in self.chembl_compound_names if name}
        name_index |= {(name,"chembl trade",None) for chembl_id,mrn,name in self.chembl_trade_names if name}
        name_index = {n:i for i,n in enumerate(name_index)}
        
        import_data = [{'id': i, 'name': name, 'table': table, 'type': name_type, 'lower': name.lower() if name else None} for (name, table, name_type),i in name_index.items()]
        for b in tqdm(batch(import_data, batch_size), total=int(len(import_data)/batch_size)+1, position=0, leave=True, desc='Names'):
            name_session.execute(self.name_store.Name.__table__.insert(), b)
            name_session.flush()
        
        
        import_data = set()
        
        for db_id,mrn,name in self.pref_names:
            if not name:
                continue
            ci = self.chembl_index.get_chembl_ident(drugbase_id=db_id, molregno=mrn, chembl_id=None)
            if ci is None:
                continue
            substance_id = substance_index[ci]
            name_id = name_index[(name,"drugbase pref",None)]
            import_data.add((substance_id, name_id))
            
        allowed_syn_types = {'USAN', 'USAN_R', 'INN', 'INN_R', 'USP', 'BAN', 'ATC', 'FDA', 'NF', 'MI', 'JAN', 'DCF', 'WHO-DD', 'BN_USP', 'BNF', 'CTGOV', 'TN', 'BN'}
        for db_id,mrn,name,name_type in self.mol_syns:
            if not name:
                continue
            if not name_type in allowed_syn_types:
                continue
            ci = self.chembl_index.get_chembl_ident(drugbase_id=db_id, molregno=mrn, chembl_id=None)
            if ci is None:
                continue
            substance_id = substance_index[ci]
            name_id = name_index[(name,"drugbase syn",name_type)]
            import_data.add((substance_id, name_id))
            
        for db_id,mrn,name in self.chem_names:
            if not name:
                continue
            ci = self.chembl_index.get_chembl_ident(drugbase_id=db_id, molregno=mrn, chembl_id=None)
            if ci is None:
                continue
            substance_id = substance_index[ci]
            name_id = name_index[(name,"drugbase chem",None)]
            import_data.add((substance_id, name_id))
            
        for db_id,mrn,name in self.dm_names:
            if not name:
                continue
            ci = self.chembl_index.get_chembl_ident(drugbase_id=db_id, molregno=mrn, chembl_id=None)
            if ci is None:
                continue
            substance_id = substance_index[ci]
            name_id = name_index[(name,"drugbase dm",None)]
            import_data.add((substance_id, name_id))
            
        for chembl_id,mrn,name in self.chembl_pref_names:
            if not name:
                continue
            ci = self.chembl_index.get_chembl_ident(drugbase_id=None, molregno=mrn, chembl_id=chembl_id)
            if ci is None:
                continue
            substance_id = substance_index[ci]
            name_id = name_index[(name,"chembl pref",None)]
            import_data.add((substance_id, name_id))
            
        for chembl_id,mrn,name in self.chembl_mol_syns:
            if not name:
                continue
            ci = self.chembl_index.get_chembl_ident(drugbase_id=None, molregno=mrn, chembl_id=chembl_id)
            if ci is None:
                continue
            substance_id = substance_index[ci]
            name_id = name_index[(name,"chembl syn",None)]
            import_data.add((substance_id, name_id))
            
        for chembl_id,mrn,name in self.chembl_compound_names:
            if not name:
                continue
            ci = self.chembl_index.get_chembl_ident(drugbase_id=None, molregno=mrn, chembl_id=chembl_id)
            if ci is None:
                continue
            substance_id = substance_index[ci]
            name_id = name_index[(name,"chembl compound",None)]
            import_data.add((substance_id, name_id))
            
        for chembl_id,mrn,name in self.chembl_trade_names:
            if not name:
                continue
            ci = self.chembl_index.get_chembl_ident(drugbase_id=None, molregno=mrn, chembl_id=chembl_id)
            if ci is None:
                continue
            substance_id = substance_index[ci]
            name_id = name_index[(name,"chembl trade",None)]
            import_data.add((substance_id, name_id))
        
        import_data = [{'id': i, 'substance_id': substance_id, 'name_id': name_id} for i,(substance_id,name_id) in enumerate(import_data)]
        for b in tqdm(batch(import_data, batch_size), total=int(len(import_data)/batch_size)+1, position=0, leave=True, desc='Assocs'):
            name_session.execute(self.name_store.Substance2Name.__table__.insert(), b)
            name_session.flush()

        name_session.commit()
        name_session.close()
    
    def close():
        self.name_session.close()
    
    def get_substance(self, substance_id):
        return self.name_session.query(self.name_store.Substance).get(substance_id)
        
    def get_name(self, name_id):
        return self.name_session.query(self.name_store.Name).get(name_id)
    
    def gen_query_index(self):
        def get_substance(substance_id):
            if not substance_id in substance_cache:
                substance_cache[substance_id] = name_session.query(self.name_store.Substance).get(substance_id)
            return substance_cache[substance_id]
        
        def get_name(name_id):
            if not name_id in name_cache:
                name_cache[name_id] = name_session.query(self.name_store.Name).get(name_id)
            return name_cache[name_id]
        
        if self.name_store is None:
            self.connect_to_namestore()

        name_session = self.name_store.SessionMaker()
        
        name_cache = {}
        substance_cache = {}
            
        self.name2substances = defaultdict(set)
        self.substance2names = defaultdict(set)
        for substance2name in tqdm(name_session.query(self.name_store.Substance2Name).all(), leave=True, position=0, desc='Generating query indexes'):
            name = get_name(substance2name.name_id)
            substance = get_substance(substance2name.substance_id)
            if name.lower:
                self.name2substances[name.lower].add((substance2name.name_id, substance2name.substance_id))
                ci = self.chembl_index.get_chembl_ident(drugbase_id=substance.drugbase_id, molregno=substance.molregno, chembl_id=substance.chembl_id)
                self.substance2names[ci].add((substance2name.name_id, substance2name.substance_id))
        self.name2substances = dict(self.name2substances)
        self.substance2names = dict(self.substance2names)

        with open(f'{self.data_dir}/name2substances.pkl', 'wb') as f:
             pickle.dump(self.name2substances, f)
        with open(f'{self.data_dir}/substance2names.pkl', 'wb') as f:
             pickle.dump({k.__tuple__():vs for k,vs in self.substance2names.items()}, f)
                
        name_session.close()
        
        self.filtered_name2substances = self.gen_filtered_query_index()
        
        with open(f'{self.data_dir}/filtered_name2substances.pkl', 'wb') as f:
             pickle.dump(self.filtered_name2substances, f)

    def gen_filtered_query_index(self):
        self.filtered_name2substances = defaultdict(set)
        for name, d in tqdm(self.name2substances.items(), leave=True, position=0):
            self.filtered_name2substances[self.filter_name(name)].update(d)
        self.filtered_name2substances = dict(self.filtered_name2substances)
        return self.filtered_name2substances
        
    def load_query_index(self):
        try:
            with open(f'{self.data_dir}/name2substances.pkl', 'rb') as f:
                self.name2substances = pickle.load(f)
            with open(f'{self.data_dir}/substance2names.pkl', 'rb') as f:
                self.substance2names = {chembl_ident.ChemblIdent(*k):vs for k,vs in pickle.load(f).items()}
            with open(f'{self.data_dir}/filtered_name2substances.pkl', 'rb') as f:
                self.filtered_name2substances = pickle.load(f)
        except:
            self.gen_query_index()
    
    def get_query_index(self):
        if (self.name2substances is None) and (self.substance2names is None) and (self.filtered_name2substances is None):
            self.load_query_index()

        return self.name2substances, self.substance2names, self.filtered_name2substances

    def query_name(self, q, filter_name=True):
        if self.name2substances is None:
            self.load_query_index()
        
        if q is None:
            return None
        
        if filter_name:
            filtered_q = self.filter_name(q)
            if bool(filtered_q):
                q = filtered_q
            else:
                q = q.lower()
            index = self.filtered_name2substances
        else:
            q = q.lower()
            index = self.name2substances
        
        if q in index:
            r = []
            for name_id, substance_id in index[q]:
                substance = self.get_substance(substance_id)
                ci = self.chembl_index.get_chembl_ident(drugbase_id=substance.drugbase_id, molregno=substance.molregno, chembl_id=substance.chembl_id)
                
                name = self.get_name(name_id)
                name = {'name': name.name, 'type': name.type, 'table': name.table}
                
                r.append({'substance': ci, 'name': name})
            
            return r
        
    def query_substance_names(self, ci):
        if self.substance2names is None:
            self.load_query_index()
        
        if ci in self.substance2names:
            r = []
            for name_id, substance_id in self.substance2names[ci]:
                name = self.get_name(name_id)
                name = {'name': name.name, 'type': name.type, 'table': name.table}
                
                r.append(name)
            
            return r
        
    @staticmethod
    def filter_name(s):
        def normalise_whitespace(s):
            s = re.sub('\s+', ' ', s)
            s = re.sub('\s+$', '', s)
            s = re.sub('^\s+', '', s)
            return s
        
        if s is None:
            return ''
        
        s = s.lower()
        if re.match('.*(\s|^)\(.*\)(\s|$).*', s):
            s = re.sub('(\s|^)\(.*\)(\s|$)', ' ', s)

            # tidy up problems that occur due to removing brackets
            s = normalise_whitespace(s)
            
            if not bool(s):
                return normalise_whitespace(s)
            
            if s[-1] == ',':
                s = normalise_whitespace(s[:-1])
            if s[0] == ',':
                s = normalise_whitespace(s[1:])
            
            if not bool(s):
                return normalise_whitespace(s)
            
            return s

        else:
            return normalise_whitespace(s)
