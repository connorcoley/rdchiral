from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdchiral.initialization import rdchiralReactants, rdchiralReaction
from rdchiral.main import rdchiralRun

# reactant = (Chem.MolFromSmiles('CCC/C=C\\CC'),)
reaction = '[C:1]/[CH:2]=[CH:3]\\[C:4]>>[C:1][CH0:2]#[CH0:3][C:4]'
reactant = '[CH3:10][CH2:11][CH2:12]/[CH:13]=[CH:14]\\[CH2:15][CH3:16]'
# products = reaction.RunReactants(reactant)
# print(Chem.CanonSmiles(Chem.MolToSmiles(products[0][0])))
# mol = Chem.MolFromSmiles(react)
# breakpoint()
reactants = rdchiralReactants(reactant)
outcomes = rdchiralRun(rdchiralReaction(reaction), reactants, keep_mapnums=True)
print(outcomes)