from __future__ import print_function
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem.rdchem import ChiralType, BondType, BondDir
from itertools import chain

PLEVEL = 2 

def vprint(level, txt):
    if PLEVEL >= level:
        print(txt)

def initialize_rxn_from_smarts(reaction_smarts):
    # Initialize reaction
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    rxn.Initialize()
    if rxn.Validate()[1] != 0:
        raise ValueError('Could not validate reaction')
    vprint(2, 'Validated rxn without errors')

    unmapped = 800
    for rct in rxn.GetReactants():
        rct.UpdatePropertyCache()
        Chem.AssignStereochemistry(rct)
        # Fill in atom map numbers
        for a in rct.GetAtoms():
            if not a.HasProp('molAtomMapNumber'):
                a.SetIntProp('molAtomMapNumber', unmapped)
                unmapped += 1
    vprint(2, 'Added {} map nums to unmapped reactants'.format(unmapped-800))
    if unmapped > 900:
        raise ValueError('Why do you have so many unmapped atoms in the template reactants?')

            ### mark chiral closures?


    return rxn

def initialize_reactants_from_smiles(reactant_smiles):
    # Initialize reactants
    reactants = Chem.MolFromSmiles(reactant_smiles)
    Chem.AssignStereochemistry(reactants)
    reactants.UpdatePropertyCache()
    vprint(2, 'Initialized reactants')
    return reactants

def run_from_text(reaction_smarts, reactant_smiles, **kwargs):
    rxn = initialize_rxn_from_smarts(reaction_smarts)
    reactants = initialize_reactants_from_smiles(reactant_smiles)

    return run(rxn, reactants, **kwargs)

def parity4(data):
    '''
    Thanks to http://www.dalkescientific.com/writings/diary/archive/2016/08/15/fragment_parity_calculation.html'''
    if data[0] < data[1]:
        if data[2] < data[3]:
            if data[0] < data[2]:
                if data[1] < data[2]:
                    return 0 # (0, 1, 2, 3) 
                else:
                    if data[1] < data[3]:
                        return 1 # (0, 2, 1, 3) 
                    else:
                        return 0 # (0, 3, 1, 2) 
            else:
                if data[0] < data[3]:
                    if data[1] < data[3]:
                        return 0 # (1, 2, 0, 3) 
                    else:
                        return 1 # (1, 3, 0, 2) 
                else:
                    return 0 # (2, 3, 0, 1) 
        else:
            if data[0] < data[3]:
                if data[1] < data[2]:
                    if data[1] < data[3]:
                        return 1 # (0, 1, 3, 2) 
                    else:
                        return 0 # (0, 2, 3, 1) 
                else:
                    return 1 # (0, 3, 2, 1) 
            else:
                if data[0] < data[2]:
                    if data[1] < data[2]:
                        return 1 # (1, 2, 3, 0) 
                    else:
                        return 0 # (1, 3, 2, 0) 
                else:
                    return 1 # (2, 3, 1, 0) 
    else:
        if data[2] < data[3]:
            if data[0] < data[3]:
                if data[0] < data[2]:
                    return 1 # (1, 0, 2, 3) 
                else:
                    if data[1] < data[2]:
                        return 0 # (2, 0, 1, 3) 
                    else:
                        return 1 # (2, 1, 0, 3) 
            else:
                if data[1] < data[2]:
                    return 1 # (3, 0, 1, 2) 
                else:
                    if data[1] < data[3]:
                        return 0 # (3, 1, 0, 2) 
                    else:
                        return 1 # (3, 2, 0, 1) 
        else:
            if data[0] < data[2]:
                if data[0] < data[3]:
                    return 0 # (1, 0, 3, 2) 
                else:
                    if data[1] < data[3]:
                        return 1 # (2, 0, 3, 1) 
                    else:
                        return 0 # (2, 1, 3, 0) 
            else:
                if data[1] < data[2]:
                    if data[1] < data[3]:
                        return 0 # (3, 0, 2, 1) 
                    else:
                        return 1 # (3, 1, 2, 0) 
                else:
                    return 0 # (3, 2, 1, 0) 

def copy_chirality(a_src, a_new):

    # Not possible to be a tetrahedral center anymore?
    if a_new.GetDegree() < 3:
        return 
    if a_new.GetDegree() == 3 and \
            any(b.GetBondType() != BondType.SINGLE for b in a_new.GetBonds()):
        return

    a_new.SetChiralTag(a_src.GetChiralTag())
    
    if not atom_chirality_matches(a_src, a_new):
        a_new.InvertChirality()

def atom_chirality_matches(a_tmp, a_mol):
    '''
    Checks for consistency in chirality between a template atom and a molecule atom.

    Also checks to see if chirality needs to be inverted in copy_chirality
    '''
    if a_mol.GetChiralTag() == ChiralType.UNSPECIFIED:
        return True # chiral template, achiral molecule -> should match, why not
    if a_tmp.GetChiralTag() == ChiralType.UNSPECIFIED:
        return False # achiral template, chiral molecule -> should not match

    isotopes_tmp = [a.GetIsotope() for a in a_tmp.GetNeighbors()]
    isotopes_new = [a.GetIsotope() for a in a_mol.GetNeighbors()]
    if len(isotopes_tmp) < 4:
        isotopes_tmp.append(0)
    if len(isotopes_mol) < 4:
        isotopes_mol.append(0)

    only_in_src = set(isotopes_tmp) - set(isotopes_mol)
    if len(only_in_src) <= 1:
        parity_matches = parity4(isotopes_tmp) == \
            parity4([i if i in isotopes_tmp else only_in_src[0] for i in isotopes_mol])
        tag_matches = a_tmp.GetChiralTag() == a_mol.GetChiralTag()
        return parity_matches == tag_matches
    else:
        return True # ambiguous case, just return for now
        # TODO: fix this?



def run(rxn, reactants, keep_isotopes=False):
    '''
    rxn = RDKit reaction
    reactants = RDKit mol with all reactants together
    '''
    final_outcomes = set()

    # To have the product atoms match reactant atoms, we
    # need to populate the Isotope field, since this field
    # gets copied over during the reaction.
    [a.SetIsotope(i+1) for (i, a) in enumerate(reactants.GetAtoms())]

    # We need to keep track of what map numbers 
    # (i.e., isotopes) correspond to which atoms
    # note: all reactant atoms must be mapped, so this is safe
    atoms_r = {a.GetIsotope(): a for a in reactants.GetAtoms()}

    # Copy reaction template so we can play around with isotopes
    for i, rct in enumerate(rxn.GetReactants()):
        if i == 0:
            template_r = rct
        else:
            template_r = AllChem.CombineMols(template_r, rct)
    for i, prd in enumerate(rxn.GetProducts()):
        if i == 0:
            template_p = prd
        else:
            template_p = AllChem.CombineMols(template_p, prd)

    atoms_rt_map = {a.GetIntProp('molAtomMapNumber'): a for a in template_r.GetAtoms() if a.HasProp('molAtomMapNumber')}
    atoms_pt_map = {a.GetIntProp('molAtomMapNumber'): a for a in template_p.GetAtoms() if a.HasProp('molAtomMapNumber')}

    # Run naive RDKit
    outcomes = rxn.RunReactants((reactants,))
    vprint(1, 'Using naive RunReactants, {} outcomes'.format(len(outcomes)))
    if not outcomes:
        return []

    for outcome in outcomes:
        ###############################################################################
        # Look for new atoms in products that were not in 
        # reactants (e.g., LGs for a retro reaction)
        vprint(1, 'Processing {}'.format([Chem.MolToSmiles(x, True) for x in outcome]))
        unmapped = 900
        for m in outcome:
            for a in m.GetAtoms():
                # Assign "map" number via isotope
                if not a.GetIsotope():
                    a.SetIsotope(unmapped)
                    unmapped += 1
                else:
                    pass
                    #TODO: copy chirality
        vprint(2, 'Added {} map numbers to product'.format(unmapped-900))
        ###############################################################################


        ###############################################################################
        # Check to see if reactants should not have been matched (based on chirality)

        # Define isotope -> reactant template atom map
        atoms_rt =  {a.GetIsotope(): atoms_rt_map[a.GetIntProp('old_mapno')] \
            for m in outcome for a in m.GetAtoms() if a.HasProp('old_mapno')}

        # Make sure each atom matches
        if not all(atom_chirality_matches(atoms_rt[i], atoms_r[i]) for i in atoms_rt):
            vprint(1, 'Chirality violated! Should not have gotten this match')
            continue
        vprint(2, 'Chirality matches! Just checked with atom_chirality_matches')

        # Check bond chirality
        # - add implicit cis to "reactant" bonds in rings swith double bond
        for b in reactants:
            if b.IsInRing() and b.GetBondType() == BondType.DOUBLE:
                pass
        #TODO: add bond chirality?

        ###############################################################################

        ###############################################################################
        # Convert product(s) to single product so that all 
        # reactions can be treated as pseudo-intramolecular
        # But! check for ring openings mistakenly split into multiple
        # This can be diagnosed by duplicate map numbers (i.e., SMILES)
        isotopes = [a.GetIsotope() for m in outcome \
            for a in m.GetAtoms() if a.GetIsotope()]
        if len(isotopes) != len(set(isotopes)): # duplicate?
            vprint(1, 'Found duplicate isotopes in product - need to stitch')
            # need to do a fancy merge
            merged_mol = Chem.RWMol(outcome[0])
            merged_iso_to_id = {a.GetIsotope(): a.GetIdx() for a in outcome[0].GetAtoms() if a.GetIsotope()}
            for j in range(1, len(outcome)):
                new_mol = outcome[j]
                for a in new_mol.GetAtoms():
                    if a.GetIsotope() not in merged_iso_to_id:
                        merged_iso_to_id[a.GetIsotope()] = merged_mol.AddAtom(a)
                for b in new_mol.GetBonds():
                    try:
                        merged_mol.AddBond(
                            merged_iso_to_id[b.GetBeginAtom().GetIsotope()],
                            merged_iso_to_id[b.GetEndAtom().GetIsotope()],
                            b.GetBondType(),
                        )
                    except: # bond already exists, will throw error
                        pass 
            outcome = merged_mol.GetMol()
            vprint(1, 'Merged editable mol, converted back to real mol')
        else:
            new_outcome = outcome[0]
            for j in range(1, len(outcome)):
                new_outcome = AllChem.CombineMols(new_outcome, outcome[j])
            outcome = new_outcome
        vprint(2, 'Converted all outcomes to single molecules')
        ###############################################################################




        ###############################################################################
        # Figure out which atoms were matched in the templates
        # note: at least one outcome at this point, so outcomes[0] is safe
        # atoms_rt and atoms_p will be outcome-specific.
        atoms_pt = {a.GetIsotope(): atoms_pt_map[a.GetIntProp('old_mapno')] \
            for a in outcome.GetAtoms() if a.HasProp('old_mapno')}
        atoms_p = {a.GetIsotope(): a for a in outcome.GetAtoms() if a.GetIsotope()}

        # Set isotopes of templates
        # note: this is okay to do within the loop, because ALL atoms must be matched
        # in the templates, so the isotopes will get overwritten every time
        # This makes it easier to check parity changes
        [a.SetIsotope(i) for (i, a) in atoms_rt.iteritems()]
        [a.SetIsotope(i) for (i, a) in atoms_pt.iteritems()]
        ###############################################################################




        ###############################################################################
        # Check for missing bonds. These are bonds that are presentt in the reactants,
        # not specified in the reactant template, and not in the product. Accidental
        # fragmentation can occur for intramolecular ring openings
        missing_bonds = []
        for b in reactants.GetBonds():
            i = b.GetBeginAtom().GetIsotope()
            j = b.GetEndAtom().GetIsotope()
            if not outcome.GetBondBetweenAtoms(atoms_p[i].GetIdx(), atoms_p[j].GetIdx()) and \
                    not template_r.GetBondBetweenAtoms(atoms_rt[i].GetIdx(), atoms_rt[j].GetIdx()):
                missing_bonds.append(i, j, b)
        if missing_bonds:
            outcome = Chem.RWMol(outcome)
            rwmol_iso_to_id = {a.GetIsotope(): a.GetIdx() for a in outcome.GetAtoms() if a.GetIsotope()}
            for (i, j, b) in missing_bonds:
                outcome.AddBond(rwmol_iso_to_id(i), rwmol_iso_to_id(j))
                new_b = outcome.GetBondBetweenAtoms(rwmol_iso_to_id(i), rwmol_iso_to_id(j))
                new_b.SetBondType(b.GetBondType())
                new_b.SetBondDir(b.GetBondDir())
                new_b.SetIsAromatic(b.GetIsAromatic())
            outcome = outcome.GetMol()
        ###############################################################################



        ###############################################################################
        # Check for chirality
        for a in outcome.GetAtoms():
            # Participants in reaction core (from reactants) will have old_mapno
            # Spectators present in reactants will have react_atom_idx
            # ...so new atoms will have neither!
            if not a.HasProp('old_mapno'):
                # Not part of the reactants template
                if not a.HasProp('react_atom_idx'):
                    # Atoms only appear in product template - their chirality
                    # should be properly instantiated by RDKit...hopefully...
                    pass
                else:
                    # Copy from reactants
                    copy_chirality(atoms_r[a.GetIsotope()], a)
            else:
                # Part of reactants and reaction core -> copy from product template
                copy_chirality(atoms_pt[a.GetIsotope()], a)
        ###############################################################################


        # Ring change?

        # Clear isotope
        if not keep_isotopes:
            [a.SetIsotope(0) for a in outcome.GetAtoms()]

        # Update property cache
        outcome.UpdatePropertyCache()

        # Uniquify via SMILES string - a little sloppy
        # Need a full SMILES->MOL->SMILES cycle to get a true canonical string
        # also, split by '.' and sort when outcome contains multiple molecules
        smiles = '.'.join(sorted(
            Chem.MolToSmiles(
                Chem.MolFromSmiles(
                    Chem.MolToSmiles(outcome, True)
                ), True 
            ).split('.')
        ))
        final_outcomes.add(smiles)

    return list(final_outcomes)


if __name__ == '__main__':
    reaction_smarts = '[C:1][OH:2]>>[C:1][O:2][C]'
    reactant_smiles = 'CC(=O)OCCCO'
    outcomes = run_from_text(reaction_smarts, reactant_smiles)
    print(outcomes)