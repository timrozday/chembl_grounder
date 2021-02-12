import pickle
import json
from collections import defaultdict
import requests

import chembl_ident
from . import gsrs_index as gi
from . import chembl_name_index as cni
from . import chembl_structure_index as csi

class ChemblGrounder():
    def __init__(self, data_dir='.', chembl_index=None, gsrs_index=None, structure_index=None, chembl_name_index=None):
        self.data_dir = data_dir
        
        if chembl_index is None:
            self.chembl_index = chembl_ident.ChemblIndexes(data_dir=self.data_dir)
        else:
            self.chembl_index = chembl_index
            
        if gsrs_index is None:
            self.gsrs_index = gi.GsrsIndex(data_dir=self.data_dir)
        else:
            self.gsrs_index = gsrs_index
        
        if structure_index is None:
            self.structure_index = csi.ChemblStructureIndex(data_dir=self.data_dir, chembl_index=self.chembl_index)
        else:
            self.structure_index = structure_index
        
        if chembl_name_index is None:
            self.chembl_name_index = cni.ChemblNameIndex(data_dir=self.data_dir, chembl_index=self.chembl_index)
        else:
            self.chembl_name_index = chembl_name_index
        self.chembl_name_index.connect_to_namestore(create=True)
        
        self.caches = {'unichem': {}, 'name': {}, 'structure': {}}
        
        
    def pick_top_evidence(self, evidence):
        type_rank = {
            'structure': 0,
            'code': 1, 
            'name': 2, 
        }
        code_rank = {
            'gsrs': 0, 
            'unichem': 1
        }
        structure_rank = {
            'drugbase': 0, 
            'chembl': 1
        }
        name_rank = {
            'drugbase pref': 0,
            'drugbase syn': 1,
            'chembl pref': 2,
            'chembl syn': 3,
            'drugbase dm': 4,
            'drugbase chem': 5,
            'chembl trade':6 ,
            'chembl compound': 7,
        }

        t = sorted(evidence.keys(), key=lambda x:type_rank[x])[0]
        es = evidence[t]

        if t == "code":
            data = [(x,(type_rank[t],
                        len(x),
                        code_rank[x[-1]['source']],
                        len(x[0]['link_data']))) for x in es]
            top, score = sorted(data, key=lambda x:x[1])[0]
            return top, score

        if t == "name":
            data = [(x,(type_rank[t],
                        len(x),
                        name_rank[x[-1]['source']],
                        len(x[0]['link_data']))) for x in es]
            top, score = sorted(data, key=lambda x:x[1])[0]
            return top, score

        if t == "structure":
            data = [(x,(type_rank[t],
                        len(x),
                        structure_rank[x[-1]['source']],  # drugbase > chembl, shortest inchi
                        len(x[0]['link_data']))) for x in es]
            top, score = sorted(data, key=lambda x:x[1])[0]
            return top, score
    
    
    def score_chembl_match(self, chembl_ident, evidence):
        inchi = self.chembl_index.get_structure(chembl_ident)
        top_evidence, (type_score, evidence_len_score, matched_entity_score, id_len_score) = self.score_evidence(evidence)
        return (type_score, 
                0 if inchi else 1, 
                evidence_len_score, 
                matched_entity_score, 
                id_len_score, 
                len(inchi) if inchi else 0)

    def pick_top_chembl_ident(self, chembl_idents):
        def score(chembl_ident):  # sort by: inchi, drugbase, inchi length
            inchi = self.structure_index.get_structure(chembl_ident)
            return (1 if chembl_ident.drugbase_id else 0, 
                    1 if inchi else 0, 
                    -len(inchi) if inchi else 0)

        return sorted(chembl_idents, key=lambda x:score(x))[-1]
    
    @staticmethod
    def rank_name_types(nt):
        ranks = {'spl': 0, 'cn': 1, 'bn': 2, 'sys': 3, 'of': 4, 'cd': 5}

        if nt in ranks:
            return ranks[nt]

        return 100

    def lookup_ingredient_unichem(self, code, src_id=14):
        response = requests.get(f'https://www.ebi.ac.uk/unichem/rest/src_compound_id/{code}/{src_id}')
        rs = json.loads(response.text)
        data = defaultdict(list)
        for r in rs:
            if not r == 'error': 
                data[int(r['src_id'])].append(r['src_compound_id'])

        if 1 in data:
            return data[1][0]

    def lookup_ingredient_structure(self, unii):
        results = []
        try:
            gsrs_inchi = self.gsrs_index.get_inchi(unii)['standardised']
        except:
            return None, None

        try:
            r = self.structure_index.query(gsrs_inchi, connectivity=True, strip=True, consistency=True, split=True)
        except Exception as e:
            print(gsrs_inchi)
            raise e
        if r:
            for chembl_ident,(consistency,evidence) in r.items():
                if consistency:
                    results.append((chembl_ident, evidence))

        return gsrs_inchi, results

    def lookup_ingredient_name(self, name, unii):
        results = defaultdict(lambda :defaultdict(set))

        names = {(name,'spl')}
        try:
            gsrs_data = self.gsrs_index.query(unii)
            if 'names' in gsrs_data:
                gsrs_names = {(name,gsrs_name_type) for name,gsrs_name_type in gsrs_data['names'] if gsrs_name_type in {'cn', 'bn', 'cd', 'of', 'sys'}}
                names.update(gsrs_names)
        except:
            pass

        for n, name_type in names:
            r = self.chembl_name_index.query_name(n, filter_name=True)
            if r:
                for v in r:
                    s = v['substance']
                    n = (v['name']['name'], v['name']['type'], v['name']['table'])
                    results[s][n].add(name_type)
        
        results = {k1:{k2:v2 for k2,v2 in v1.items()} for k1,v1 in results.items()}
        return results

    def query(self, name, unii):
        ingredient_matches_evidence = defaultdict(lambda :defaultdict(list))

        # structure matches
        if not unii in self.caches['structure']:
            self.caches['structure'][unii] = self.lookup_ingredient_structure(unii)
        inchi, structure_matches = self.caches['structure'][unii]

        if structure_matches:
            for chembl_ident, score in structure_matches:
                if chembl_ident.chembl_id:
                    t = {'source': 'chembl', 'link_type': 'chembl_id', 'link_data': chembl_ident.chembl_id}
                elif chembl_ident.molregno:
                    t = {'source': 'chembl', 'link_type': 'molregno', 'link_data': chembl_ident.molregno}
                elif chembl_ident.drugbase_id:
                    t = {'source': 'chembl', 'link_type': 'drugbase_id', 'link_data': chembl_ident.drugbase_id}
                evidence_dict = [
                    {'source': 'spl', 'link_type': 'unii', 'link_data': unii}, 
                    {'source': 'gsrs', 'link_type': 'inchi', 'link_data': inchi}, 
                    t
                ]
                ingredient_matches_evidence[chembl_ident]['structure'].append(evidence_dict)

        # UniChem code matches
        if not unii in self.caches['unichem']:
            self.caches['unichem'][unii] = self.lookup_ingredient_unichem(unii)
        unichem_match = self.caches['unichem'][unii]

        if unichem_match:
            chembl_id = unichem_match
            chembl_ident = self.chembl_index.get_chembl_ident(chembl_id=chembl_id)
            evidence_dict = [{'source': 'spl', 'link_type': 'unii', 'link_data': unii}, {'source': 'unichem', 'link_type': 'chembl_id', 'link_data': chembl_id}]
            ingredient_matches_evidence[chembl_ident]['code'].append(evidence_dict)

        # Name matches
        if not (name, unii) in self.caches['name']:
            self.caches['name'][(name, unii)] = self.lookup_ingredient_name(name, unii)
        name_matches = self.caches['name'][(name, unii)]

        if name_matches:
            for chembl_ident,matched_names in name_matches.items():

                # get the preferred ChEMBL identifier
                chembl_code_type, chembl_code = next(
                    ((t,i) for t,i in [
                        ("drugbase_id",chembl_ident.drugbase_id), 
                        ("chembl_id",chembl_ident.chembl_id), 
                        ("molregno",chembl_ident.molregno)
                    ] if not i is None), 
                    None
                )

                for name_match, name_types in matched_names.items():
                    match_name, match_name_type, match_name_table = name_match
                    name_type = sorted(list(name_types), key=self.rank_name_types)[0]

                    if name_type == 'spl':
                        evidence_dict = [
                            {'source': 'spl', 'link_type': 'name', 'link_data': match_name}, 
                            {'source': match_name_table, 'link_type': chembl_code_type, 'link_data': chembl_code}
                        ]
                    else:
                        evidence_dict = [
                            {'source': 'spl', 'link_type': 'unii', 'link_data': unii}, 
                            {'source': 'gsrs', 'link_type': name_type, 'link_data': name_match}, 
                            {'source': match_name_table, 'link_type': chembl_code_type, 'link_data': chembl_code}
                        ]

                    ingredient_matches_evidence[chembl_ident]['name'].append(evidence_dict)
        
        ingredient_matches_evidence = {k1:{k2:v2 for k2,v2 in v1.items()} for k1,v1 in ingredient_matches_evidence.items()}

        # pick top
        if ingredient_matches_evidence:
            chembl_ident_evidence = {c:self.pick_top_evidence(e) for c,e in ingredient_matches_evidence.items()}
            top_evidence_score = min({s for c,(t,s) in chembl_ident_evidence.items()})
            top_chembl_ident = self.pick_top_chembl_ident({c for c,(t,s) in chembl_ident_evidence.items() if s==top_evidence_score})

            return top_chembl_ident, ingredient_matches_evidence
        
        else:
            return None, None
