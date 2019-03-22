from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)

import json
import gzip
import hashlib
import pandas as pd
from rdkit import Chem
from joblib import Parallel, delayed
from time import time


import template_extractor

def can_parse(rsmi):
    react, spec, prod = rsmi.split('>')
    if Chem.MolFromSmiles(react) and Chem.MolFromSmiles(prod):
        return True
    else:
        return False
    
t0 = time()

uspto = pd.read_csv('data/1976_Sep2016_USPTOgrants_smiles.rsmi', sep='\t')

uspto['ReactionSmiles'] = uspto['ReactionSmiles'].str.split(' ', expand=True)[0]
split_smiles = uspto['ReactionSmiles'].str.split('>', expand=True)
uspto['reactants'] = split_smiles[0]
uspto['spectators'] = split_smiles[1]
uspto['products'] = split_smiles[2]

parsable = Parallel(n_jobs=-1, verbose=1)(delayed(can_parse)(rsmi) for rsmi in uspto['ReactionSmiles'].values)
# parsable = uspto['ReactionSmiles'].map(can_parse)

uspto = uspto[parsable]
print('{} parsable reactions'.format(len(uspto)))

hexhash = (uspto['ReactionSmiles']+uspto['PatentNumber']).apply(lambda x: hashlib.sha256(x.encode('utf-8')).hexdigest())

uspto['source'] = 'uspto'
uspto['source_id'] = hexhash

uspto = uspto.reset_index().rename(columns={'index': '_id'})

reactions = uspto[['_id', 'reactants', 'products', 'spectators', 'source', 'source_id']]

reactions.to_json('data/uspto.reactions.json.gz', orient='records', compression='gzip')

with gzip.open('data/uspto.reactions.json.gz') as f:
    reactions = json.load(f)

def extract(reaction):
    try:
        return template_extractor.extract_from_reaction(reaction)
    except KeyboardInterrupt:
        print('Interrupted')
        raise KeyboardInterrupt
    except Exception as e:
        print(e)
        return {'reaction_id': reaction['_id']}

templates = Parallel(n_jobs=-1, verbose=4)(delayed(extract)(reaction) for reaction in reactions)

with gzip.open('data/uspto.templates.json.gz', 'w') as f:
    json.dump(templates, f)
    
print('elapsed seconds: {}'.format(int(time()-t0)))