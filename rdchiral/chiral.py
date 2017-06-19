from __future__ import print_function
from rdkit.Chem.rdchem import ChiralType, BondType, BondDir

from rdchiral.utils import vprint, parity4

def template_atom_could_have_been_tetra(a, strip_if_spec=False):
    '''
    Could this atom have been a tetrahedral center?
    If yes, template atom is considered achiral and will not match a chiral rct
    If no, the tempalte atom is auxilliary and we should not use it to remove
    a matched reaction. For example, a fully-generalized terminal [C:1] 
    '''

    if a.HasProp('tetra_possible'):
        return a.GetBoolProp('tetra_possible')
    if a.GetDegree() < 3 or (a.GetDegree() == 3 and 'H' not in a.GetSmarts()):
        a.SetBoolProp('tetra_possible', False)
        if strip_if_spec: # Clear chiral tag in case improperly set
            a.SetChiralTag(ChiralType.CHI_UNSPECIFIED)
        return False 
    a.SetBoolProp('tetra_possible', True)
    return True 



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

    # When there are fewer than 3 heavy neighbors, chirality is ambiguous...
    if len(isotopes_tmp) < 3 or len(isotopes_mol) < 3:
        return True

    # Degree of 3 -> remaining atom is a hydrogen, add to list
    if len(isotopes_tmp) < 4:
        isotopes_tmp.append(-1) # H
    if len(isotopes_mol) < 4:
        isotopes_mol.append(-1) # H

    try:
        vprint(10, str(isotopes_tmp))
        vprint(10, str(isotopes_mol))
        vprint(10, str(a_tmp.GetChiralTag()))
        vprint(10, str(a_mol.GetChiralTag()))
        only_in_src = [i for i in isotopes_tmp if i not in isotopes_mol][::-1] # reverse for popping
        only_in_mol = [i for i in isotopes_mol if i not in isotopes_tmp]
        if len(only_in_src) <= 1 and len(only_in_mol) <= 1:
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
            vprint(2, 'Isotope {} chiral match? Based on isotope lists, ambiguous -> True', a_tmp.GetIsotope())
            return True # ambiguous case, just return for now
            # TODO: fix this?

    except IndexError as e:
        print(a_tmp.GetPropsAsDict())
        print(a_mol.GetPropsAsDict())
        print(a_tmp.GetChiralTag())
        print(a_mol.GetChiralTag())
        print(str(e))
        print(str(isotopes_tmp))
        print(str(isotopes_mol))
        raise KeyError('Pop from empty set - this should not happen!')
