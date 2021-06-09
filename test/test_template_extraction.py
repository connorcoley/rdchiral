from rdkit import Chem
from rdkit.Chem import AllChem
import os, sys, json
import re
sys.path = [os.path.dirname(os.path.dirname((__file__)))] + sys.path 

from rdchiral.main import rdchiralReaction, rdchiralReactants, rdchiralRunText, rdchiralRun
from rdchiral.template_extractor import extract_from_reaction

#used stupid fix to avoid errors -- replace all "\" with "!"
with open(os.path.join(os.path.dirname(__file__), 'test_rdchiral_extraction.json'), 'r') as fid:
    test_cases = json.load(fid)


#helper functions
def union(list_of_lists):
    u = set()
    for ls in list_of_lists:
      u = u.union(set(ls))
    return u


def intersection(list_of_lists):
    i = set(list_of_lists[0])
    for ls in list_of_lists[1:]:
      i = i.intersection(set(ls))
    return i

def unmap(smiles):
  mol = Chem.MolFromSmiles(smiles)
  [a.SetAtomMapNum(0) for a in mol.GetAtoms()]
  return Chem.MolToSmiles(mol)

def smiles_processing(smiles):
    smiles = re.sub('[/\\\]C=N[/\\\]', 'C=N', smiles)
    return smiles

def simplify_stereo(template):
    """
    Given reactants, products, and a reaction template,
    Remove 
    """
    smarts = template["reaction_smarts"]
    p_smarts, _, r_smarts = smarts.split(">")
    
    reactants = Chem.MolFromSmiles(smiles_processing(template["reactants"]))
    products = Chem.MolFromSmiles(smiles_processing(template["products"]))
    
#     print(Chem.MolToSmiles(reactants), Chem.MolToSmiles(products))
    r_smarts = AllChem.MolFromSmarts(r_smarts)
    p_smarts = AllChem.MolFromSmarts(p_smarts)

    # keep chirality for these

    react_ctr_idxs = union(reactants.GetSubstructMatches(r_smarts, useChirality=True))
    prod_ctr_idxs = union(products.GetSubstructMatches(p_smarts, useChirality=True))

    react_ctr_maps = [
      reactants.GetAtomWithIdx(i).GetAtomMapNum() for i in react_ctr_idxs
    ]
    prod_ctr_maps = [products.GetAtomWithIdx(i).GetAtomMapNum() for i in prod_ctr_idxs]

    react_maps = intersection([react_ctr_maps, prod_ctr_maps])

    react_idxs = {
      at.GetAtomMapNum(): at.GetIdx()
      for at in reactants.GetAtoms()
      if at.GetAtomMapNum() not in react_maps
    }
    prod_idxs = {
      at.GetAtomMapNum(): at.GetIdx()
      for at in products.GetAtoms()
      if at.GetAtomMapNum() not in react_maps
    }

    mapnums = intersection((prod_idxs.keys(), react_idxs.keys()))
    for mapnum in mapnums:
        
      #remove @ @@ information from atoms not in reaction center
      distant_react_atom = reactants.GetAtomWithIdx(react_idxs[mapnum])
      distant_prod_atom = products.GetAtomWithIdx(prod_idxs[mapnum])
      distant_react_atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
      distant_prod_atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
      
      #remove E/Z information from bonds not in reaction center
      for bond in distant_react_atom.GetBonds():
        if bond.GetBeginAtomIdx() in mapnums and bond.GetEndAtomIdx() in mapnums:
            bond.SetStereo(Chem.rdchem.BondStereo.STEREONONE)
              
      for bond in distant_prod_atom.GetBonds():
        if bond.GetBeginAtomIdx() in mapnums and bond.GetEndAtomIdx() in mapnums:
            bond.SetStereo(Chem.rdchem.BondStereo.STEREONONE)

    return Chem.MolToSmiles(reactants), Chem.MolToSmiles(products)


all_passed = True
for i, test_case in enumerate(test_cases):

    print('\n# Test {:2d}/{} -- reaction id: {}'.format(i+1, len(test_cases), test_case['_id']))
    #restoring original syntax
    test_case['expected_template'] = test_case['expected_template']
    test_case['reactants'] = test_case['reactants']
    test_case['products'] = test_case['products']
    try:
        extracted = extract_from_reaction(test_case)['reaction_smarts']
        if extracted == test_case['expected_template']:
            print('    from text: passed template extraction')
        else:
            all_passed = False
            print('    from text: failed template extraction')
            print('    Expected : {}'.format(test_case['expected_template']))
            print('    Returned : {}'.format(extracted))
        
        #simplify the stereochemistry of the molecules based on the newly extracted template
        test_case['reaction_smarts'] = extracted
        simp_reactants, simp_products = simplify_stereo(test_case)
        pred_reactants = rdchiralRunText('('+extracted.replace('>>',')>>'), unmap(simp_products), combine_enantiomers = False)
        
        for pred_reactant in pred_reactants:
          no_matches = True
          if sorted(unmap(pred_reactant).split('.')) == sorted(unmap(simp_reactants).split('.')):
              print('    from text: applying template recovers reactants')
              no_matches = False
              break
        if no_matches:
          all_passed = False
          print('    from text: template does not recover reactants')
          print ('Expected:', unmap(simp_reactants) ,'\nReturned:', pred_reactants)

    except KeyError:
        all_passed = False
        print('    from text: Failed to extract template')
        
        
        

all_passed = 'All passed! -- Extracted templates matched the expected templates or the extracted template correctly recovers the reactants' if all_passed else 'Failed!'
print('\n# Final result: {}'.format(all_passed))