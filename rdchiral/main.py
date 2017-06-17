from __future__ import print_function
import sys 
sys.path.append('../../')

import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem.rdchem import ChiralType, BondType, BondDir
from itertools import chain
from rdchiral.utils.parity4 import parity4

PLEVEL = 10

def vprint(level, txt, *args):
    if PLEVEL >= level:
        print(txt.format(*args))

def initialize_rxn_from_smarts(reaction_smarts):
    # Initialize reaction
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    rxn.Initialize()
    if rxn.Validate()[1] != 0:
        raise ValueError('Could not validate reaction')
    vprint(2, 'Validated rxn without errors')

    unmapped = 700
    for rct in rxn.GetReactants():
        rct.UpdatePropertyCache()
        Chem.AssignStereochemistry(rct)
        # Fill in atom map numbers
        for a in rct.GetAtoms():
            if not a.HasProp('molAtomMapNumber'):
                a.SetIntProp('molAtomMapNumber', unmapped)
                unmapped += 1
    vprint(2, 'Added {} map nums to unmapped reactants', unmapped-700)
    if unmapped > 800:
        raise ValueError('Why do you have so many unmapped atoms in the template reactants?')

    return rxn

def initialize_reactants_from_smiles(reactant_smiles):
    # Initialize reactants
    reactants = Chem.MolFromSmiles(reactant_smiles)
    Chem.AssignStereochemistry(reactants, flagPossibleStereoCenters=True)
    reactants.UpdatePropertyCache()
    vprint(2, 'Initialized reactants, assigned stereochem with flagpossiblestereocenters')
    return reactants

def run_from_text(reaction_smarts, reactant_smiles, **kwargs):
    rxn = initialize_rxn_from_smarts(reaction_smarts)
    reactants = initialize_reactants_from_smiles(reactant_smiles)

    return run(rxn, reactants, **kwargs)

def copy_chirality(a_src, a_new):

    # Not possible to be a tetrahedral center anymore?
    if a_new.GetDegree() < 3:
        return 
    if a_new.GetDegree() == 3 and \
            any(b.GetBondType() != BondType.SINGLE for b in a_new.GetBonds()):
        return

    vprint(3, 'For isotope {}, copying src {} chirality tag to new',
        a_src.GetIsotope(), a_src.GetChiralTag())
    a_new.SetChiralTag(a_src.GetChiralTag())
    
    if not atom_chirality_matches(a_src, a_new):
        vprint(3, 'For isotope {}, inverting chirality', a_new.GetIsotope())
        a_new.InvertChirality()

def template_atom_could_have_been_tetra(a):
    '''
    Could this atom have been a tetrahedral center?
    If yes, template atom is considered achiral and will not match a chiral rct
    If no, the tempalte atom is auxilliary and should not rule out a product set
    '''

    if a.HasProp('tetra_possible'):
        return a.GetBoolProp('tetra_possible')
    if a.GetDegree() < 3 or (a.GetDegree() == 3 and 'H' not in a.GetSmarts()):
        a.SetBoolProp('tetra_possible', False)
        return False 
    a.SetBoolProp('tetra_possible', True)
    return True 

def atom_chirality_matches(a_tmp, a_mol):
    '''
    Checks for consistency in chirality between a template atom and a molecule atom.

    Also checks to see if chirality needs to be inverted in copy_chirality
    '''
    if a_mol.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
        if a_tmp.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
            vprint(3, 'atom {} is achiral & achiral -> match', a_mol.GetIsotope())
            return True # achiral template, achiral molecule -> match
        # What if the template was chiral, but the reactant isn't just due to symmetry?
        if not a_mol.HasProp('_ChiralityPossible'):
            # It's okay to make a match, as long as the product is achiral (even
            # though the product template will try to impose chirality)
            vprint(3, 'atom {} is specified in template, but cant possibly be chiral in mol', a_mol.GetIsotope())
            return True

        # TODO: figure out if we want this behavior - should a chiral template
        # be applied to an achiral molecule? For the retro case, if we have
        # a retro reaction that requires a specific stereochem, return False;
        # however, there will be many cases where the reaction would probably work
        vprint(3, 'atom {} is achiral in mol, but specified in template', a_mol.GetIsotope())
        return False
    if a_tmp.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
        vprint(3, 'Reactant {} atom chiral, rtemplate achiral...', a_tmp.GetIsotope())
        if template_atom_could_have_been_tetra(a_tmp):
            vprint(3, '...and that atom could have had its chirality specified! no_match')
            return False
        vprint(3, '...but the rtemplate atom could not have had chirality specified, match anyway')
        return True

    isotopes_tmp = [a.GetIsotope() for a in a_tmp.GetNeighbors()]
    isotopes_mol = [a.GetIsotope() for a in a_mol.GetNeighbors()]
    if len(isotopes_tmp) < 4:
        isotopes_tmp.append(-1) # H
    if len(isotopes_mol) < 4:
        isotopes_mol.append(-1) # H

    try:
        vprint(10, str(isotopes_tmp))
        vprint(10, str(isotopes_mol))
        vprint(10, str(a_tmp.GetChiralTag()))
        vprint(10, str(a_mol.GetChiralTag()))
        only_in_src = set(isotopes_tmp) - set(isotopes_mol)
        if len(only_in_src) <= 1:
            tmp_parity = parity4(isotopes_tmp)
            mol_parity = parity4([i if i in isotopes_tmp else only_in_src.pop() for i in isotopes_mol])
            vprint(10, str(tmp_parity))
            vprint(10, str(mol_parity))
            parity_matches = tmp_parity == mol_parity
            tag_matches = a_tmp.GetChiralTag() == a_mol.GetChiralTag()
            chirality_matches = parity_matches == tag_matches
            vprint(2, 'Isotope {} chiral match? {}', a_tmp.GetIsotope(), chirality_matches)
            return chirality_matches
        else:
            return True # ambiguous case, just return for now
            # TODO: fix this?
    except KeyError:
        print(isotopes_tmp)
        print(isotopes_mol)
        print(only_in_src)
        raise KeyError('Pop from empty set')



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

    ###############################################################################
    # Run naive RDKit

    # Strip chirality first
    chi_tags = [a.GetChiralTag() for a in reactants.GetAtoms()]
    [a.SetChiralTag(ChiralType.CHI_UNSPECIFIED) for a in reactants.GetAtoms()]
    vprint(3, 'Stripped chiral tags from reactants')
    
    outcomes = rxn.RunReactants((reactants,))
    vprint(2, 'Using naive RunReactants on {}, {} outcomes', 
        Chem.MolToSmiles(reactants, True), len(outcomes))

    # Restore chiral tags
    [a.SetChiralTag(chi_tags[i]) for (i, a) in enumerate(reactants.GetAtoms())]
    vprint(3, 'Restored chiral tags to reactants')
    if not outcomes:
        return []


    ###############################################################################

    for outcome in outcomes:
        ###############################################################################
        # Look for new atoms in products that were not in 
        # reactants (e.g., LGs for a retro reaction)
        vprint(2, 'Processing {}', str([Chem.MolToSmiles(x, True) for x in outcome]))
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
        vprint(2, 'Added {} map numbers to product', unmapped-900)
        ###############################################################################


        ###############################################################################
        # Check to see if reactants should not have been matched (based on chirality)

        # Define isotope -> reactant template atom map
        vprint(3, 'Isotopes with old_mapno field: {}'.format(
            [a.GetIsotope() for m in outcome for a in m.GetAtoms() if a.HasProp('old_mapno')]
        ))
        atoms_rt =  {a.GetIsotope(): atoms_rt_map[a.GetIntProp('old_mapno')] \
            for m in outcome for a in m.GetAtoms() if a.HasProp('old_mapno')}

        # Set isotopes of reactant template
        # note: this is okay to do within the loop, because ALL atoms must be matched
        # in the templates, so the isotopes will get overwritten every time
        [a.SetIsotope(i) for (i, a) in atoms_rt.iteritems()]

        # Make sure each atom matches
        if not all(atom_chirality_matches(atoms_rt[i], atoms_r[i]) for i in atoms_rt):
            vprint(2, 'Chirality violated! Should not have gotten this match')
            continue
        vprint(2, 'Chirality matches! Just checked with atom_chirality_matches')

        # Check bond chirality
        # - add implicit cis to "reactant" bonds in rings swith double bond
        for b in reactants.GetBonds():
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
                    except RuntimeError: # bond already exists, will throw error
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

        # Set isotopes of product template
        # note: this is okay to do within the loop, because ALL atoms must be matched
        # in the templates, so the isotopes will get overwritten every time
        # This makes it easier to check parity changes
        [a.SetIsotope(i) for (i, a) in atoms_pt.iteritems()]
        ###############################################################################



        ###############################################################################
        # Check for missing bonds. These are bonds that are present in the reactants,
        # not specified in the reactant template, and not in the product. Accidental
        # fragmentation can occur for intramolecular ring openings
        missing_bonds = []
        for b in reactants.GetBonds():
            i = b.GetBeginAtom().GetIsotope()
            j = b.GetEndAtom().GetIsotope()
            if i in atoms_p and j in atoms_p:
                # atoms from reactant bond show up in product
                if not outcome.GetBondBetweenAtoms(atoms_p[i].GetIdx(), atoms_p[j].GetIdx()):
                    #...but there is not a bond in the product between those atoms
                    if i not in atoms_rt or j not in atoms_rt or not template_r.GetBondBetweenAtoms(atoms_rt[i].GetIdx(), atoms_rt[j].GetIdx()):
                        # the reactant template did not specify a bond between those atoms (e.g., intentionally destroy)
                        missing_bonds.append((i, j, b))
        if missing_bonds:
            vprint(1, 'Product is missing non-reacted bonds that were present in reactants!')
            outcome = Chem.RWMol(outcome)
            rwmol_iso_to_id = {a.GetIsotope(): a.GetIdx() for a in outcome.GetAtoms() if a.GetIsotope()}
            for (i, j, b) in missing_bonds:
                outcome.AddBond(rwmol_iso_to_id[i], rwmol_iso_to_id[j])
                new_b = outcome.GetBondBetweenAtoms(rwmol_iso_to_id[i], rwmol_iso_to_id[j])
                new_b.SetBondType(b.GetBondType())
                new_b.SetBondDir(b.GetBondDir())
                new_b.SetIsAromatic(b.GetIsAromatic())
            outcome = outcome.GetMol()
        else:
            vprint(3, 'No missing bonds')
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
                    vprint(4, 'Atom {} created by product template, should have right chirality', a.GetIsotope())
                else:
                    vprint(4, 'Atom {} outside of template, copy chirality from reactants', a.GetIsotope())
                    copy_chirality(atoms_r[a.GetIsotope()], a)
            else:
                # Part of reactants and reaction core
                if template_atom_could_have_been_tetra(atoms_rt[a.GetIsotope()]):
                    vprint(3, 'Atom {} was in rct template (could have been tetra)', a.GetIsotope())
                    if template_atom_could_have_been_tetra(atoms_pt[a.GetIsotope()]):
                        vprint(3, 'Atom {} in product template could have been tetra, too', a.GetIsotope())
                        # Was the product template specified?
                        if atoms_pt[a.GetIsotope()].GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
                            # No, leave unspecified in product
                            vprint(3, '...but it is not specified in product, so destroy chirality')
                            a.SetChiralTag(ChiralType.CHI_UNSPECIFIED)
                        else:
                            # Yes
                            vprint(3, '...and product is specified')
                            # Was the reactant template specified?
                            if atoms_rt[a.GetIsotope()].GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
                                # No, so the reaction introduced chirality
                                vprint(3, '...but reactant template was not, so copy from product template')
                                copy_chirality(atoms_pt[a.GetIsotope()], a)
                            else:
                                # Yes, so we need to check if chirality should be preserved or inverted
                                vprint(3, '...and reactant template was, too! copy from reactants')
                                copy_chirality(atoms_r[a.GetIsotope()], a)
                                if not atom_chirality_matches(atoms_pt[a.GetIsotope()], atoms_rt[a.GetIsotope()]):
                                    vprint(3, 'but! reactant template and product template have opposite stereochem, so invert')
                                    a.InvertChirality()
                    else:
                        vprint(3, 'If reactant template could have been ' +
                            'chiral, but the product template could not, then we dont need ' +
                            'to worry about specifying product atom chirality')

                else:
                    vprint(3, 'Atom {} could not have been chiral in reactant template', a.GetIsotope())
                    if not template_atom_could_have_been_tetra(atoms_pt[a.GetIsotope()]):
                        vprint(3, 'Atom {} also could not have been chiral in product template', a.GetIsotope())
                        vprint(3, '...so, copy chirality from reactant instead')
                        copy_chirality(atoms_r[a.GetIsotope()], a)
                    else:
                        vprint(3, 'Atom could/does have product template chirality!', a.GetIsotope())
                        vprint(3, '...so, copy chirality from product template')
                        copy_chirality(atoms_pt[a.GetIsotope()], a)
                    print('New chiral tag {}'.format(a.GetChiralTag()))
                    print('new mol: {}'.format(Chem.MolToSmiles(outcome, True)))
        vprint(2, 'After attempting to re-introduce chirality, outcome = {}',
            Chem.MolToSmiles(outcome, True))
        ###############################################################################

        # Clear isotope
        if not keep_isotopes:
            [a.SetIsotope(0) for a in outcome.GetAtoms()]

        # Update property cache
        outcome.UpdatePropertyCache()

        # Uniquify via SMILES string - a little sloppy
        # Need a full SMILES->MOL->SMILES cycle to get a true canonical string
        # also, split by '.' and sort when outcome contains multiple molecules
        smiles = Chem.MolToSmiles(outcome, True)
        outcome = Chem.MolFromSmiles(smiles)
        if outcome is None:
            vprint(1, '~~ could not parse self?')
            vprint(1, 'Attempted SMILES: {}', smiles)
            continue
        smiles = '.'.join(sorted(Chem.MolToSmiles(outcome, True).split('.')))
        final_outcomes.add(smiles)

    return list(final_outcomes)


if __name__ == '__main__':
    reaction_smarts = '[C:1][OH:2]>>[C:1][O:2][C]'
    reactant_smiles = 'CC(=O)OCCCO'
    outcomes = run_from_text(reaction_smarts, reactant_smiles)
    print(outcomes)