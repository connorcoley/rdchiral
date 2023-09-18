import os, sys, json
sys.path = [os.path.dirname(os.path.dirname((__file__)))] + sys.path 

from rdchiral.main import rdchiralReaction, rdchiralReactants, rdchiralRunText, rdchiralRun
from rdkit import Chem

with open(os.path.join(os.path.dirname(__file__), 'test_atom_mapping_cases.json'), 'r') as fid:
    test_cases = json.load(fid)

all_passed = True

def canonicalize_outcomes(outcomes):
    ''' Convert all SMILES in a list of outcomes to the canonical form '''
    return list(map(lambda x: Chem.CanonSmiles(x), outcomes))

for i, test_case in enumerate(test_cases):

    print('\n# Test {:2d}/{}'.format(i+1, len(test_cases)))

    # Initialize with test case SMILES/SMARTS
    reaction_smarts = test_case['smarts']
    reactant_smiles = test_case['smiles']
    reactants = rdchiralReactants(reactant_smiles, custom_reactant_mapping=True)
    expected = canonicalize_outcomes(test_case['expected'])

    # Test rdchiralRunText
    if canonicalize_outcomes(rdchiralRunText(reaction_smarts, reactant_smiles, custom_reactant_mapping=True, keep_mapnums=True)) == expected:
        print('    from text: passed')
    else:
        print('    from text: failed')
        all_passed = False

    # Pre-initialize & repeat with rdChiralRun
    rxn = rdchiralReaction(reaction_smarts)
    if all(canonicalize_outcomes(rdchiralRun(rxn, reactants, keep_mapnums=True)) == expected for j in range(3)):
        print('    from init: passed')
    else:
        print('    from init: failed')
        all_passed = False

all_passed = 'All passed!' if all_passed else 'Failed!'
print('\n# Final result: {}'.format(all_passed))