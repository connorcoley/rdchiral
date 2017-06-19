from __future__ import print_function
import sys 
import os

import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem.rdchem import ChiralType, BondType, BondDir

from rdchiral.utils import vprint
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdchiral.chiral import template_atom_could_have_been_tetra, copy_chirality, atom_chirality_matches
from rdchiral.clean import canonicalize_outcome_smiles, combine_enantiomers_into_racemic

def rdchiralRunText(reaction_smarts, reactant_smiles, **kwargs):
    '''Run from SMARTS string and SMILES string. This is NOT recommended
    for library application, since initialization is pretty slow. You should
    separately initialize the template and molecules and call run()'''
    rxn = rdchiralReaction(reaction_smarts)
    reactants = rdchiralReactants(reactant_smiles)
    return rdchiralRun(rxn, reactants, **kwargs)

def rdchiralRun(rxn, reactants, keep_isotopes=False, combine_enantiomers=True):
    '''
    rxn = rdchiralReaction (rdkit reaction + auxilliary information)
    reactants = rdchiralReactants (rdkit mol + auxilliary information)

    note: there is a fair amount of initialization (assigning stereochem), most
    importantly assigning isotope numbers to the reactant atoms. It is 
    HIGHLY recommended to use the custom classes for initialization.
    '''

    final_outcomes = set()

    # We need to keep track of what map numbers 
    # (i.e., isotopes) correspond to which atoms
    # note: all reactant atoms must be mapped, so this is safe
    atoms_r = reactants.atoms_r

    # Copy reaction template so we can play around with isotopes
    template_r, template_p = rxn.template_r, rxn.template_p

    # Get molAtomMapNum->atom dictionary for tempalte reactants and products
    atoms_rt_map = rxn.atoms_rt_map
    atoms_pt_map = rxn.atoms_pt_map

    ###############################################################################
    # Run naive RDKit on ACHIRAL version of molecules

    outcomes = rxn.rxn.RunReactants((reactants.reactants_achiral,))
    vprint(2, 'Using naive RunReactants, {} outcomes', len(outcomes))
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
        vprint(2, 'Added {} map numbers to product', unmapped-900)
        ###############################################################################


        ###############################################################################
        # Check to see if reactants should not have been matched (based on chirality)

        # Define isotope -> reactant template atom map
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
        #TODO: add bond chirality considerations to exclude improper matches

        ###############################################################################



        ###############################################################################
        # Convert product(s) to single product so that all 
        # reactions can be treated as pseudo-intramolecular
        # But! check for ring openings mistakenly split into multiple
        # This can be diagnosed by duplicate map numbers (i.e., SMILES)

        isotopes = [a.GetIsotope() for m in outcome for a in m.GetAtoms() if a.GetIsotope()]
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
                    bi = b.GetBeginAtom().GetIsotope()
                    bj = b.GetEndAtom().GetIsotope()
                    vprint(10, 'stitching bond between {} and {} in stich has chirality {}, {}'.format(
                        bi, bj, b.GetStereo(), b.GetBondDir()
                    ))
                    if not merged_mol.GetBondBetweenAtoms(
                            merged_iso_to_id[bi], merged_iso_to_id[bj]):
                        merged_mol.AddBond(merged_iso_to_id[bi],
                            merged_iso_to_id[bj], b.GetBondType())
                        merged_mol.GetBondBetweenAtoms(
                            merged_iso_to_id[bi], merged_iso_to_id[bj]
                        ).SetStereo(b.GetStereo())
                        merged_mol.GetBondBetweenAtoms(
                            merged_iso_to_id[bi], merged_iso_to_id[bj]
                        ).SetBondDir(b.GetBondDir())
            outcome = merged_mol.GetMol()
            vprint(1, 'Merged editable mol, converted back to real mol, {}', Chem.MolToSmiles(outcome, True))
        else:
            new_outcome = outcome[0]
            for j in range(1, len(outcome)):
                new_outcome = AllChem.CombineMols(new_outcome, outcome[j])
            outcome = new_outcome
        vprint(2, 'Converted all outcomes to single molecules')
        ###############################################################################




        ###############################################################################
        # Figure out which atoms were matched in the templates
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
        for (i, j, b) in reactants.bonds_by_isotope:
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


        # Now that we've fixed any bonds, connectivity is set. This is a good time
        # to udpate the property cache, since all that is left is fixing atom/bond
        # stereochemistry.
        try:
            outcome.UpdatePropertyCache()
        except ValueError as e: 
            vprint(1, '{}, {}'.format(Chem.MolToSmiles(outcome, True), e))
            continue


        ###############################################################################
        # Correct tetra chirality in the outcome

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
                        # Reactant template chiral, product template not - the
                        # reaction is supposed to destroy chirality, so leave
                        # unspecified
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
                    
            vprint(3, 'New chiral tag {}', a.GetChiralTag())
        vprint(2, 'After attempting to re-introduce chirality, outcome = {}',
            Chem.MolToSmiles(outcome, True))
        ###############################################################################


        ###############################################################################
        # Correct bond directionality in the outcome
        # TODO


        # Clear isotope
        if not keep_isotopes:
            [a.SetIsotope(0) for a in outcome.GetAtoms()]

        # Canonicalize
        smiles = canonicalize_outcome_smiles(outcome)
        if smiles is not None:
            final_outcomes.add(smiles)

    ###############################################################################
    # One last fix for consolidating multiple stereospecified products...
    if combine_enantiomers:
        final_outcomes = combine_enantiomers_into_racemic(final_outcomes)
    ###############################################################################

    return list(final_outcomes)


if __name__ == '__main__':
    reaction_smarts = '[C:1][OH:2]>>[C:1][O:2][C]'
    reactant_smiles = 'CC(=O)OCCCO'
    outcomes = rdchiralRunText(reaction_smarts, reactant_smiles)
    print(outcomes)