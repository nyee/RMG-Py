#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module provides classes and methods for working with molecular substructure
groups. These enable molecules to be searched for common motifs (e.g.
reaction sites).
"""

import cython

from .graph import Vertex, Edge, Graph
from .atomtype import atomTypes, allElements, nonSpecifics, getFeatures
import rmgpy.molecule.molecule as mol
from copy import deepcopy, copy

################################################################################

class ActionError(Exception):
    """
    An exception class for errors that occur while applying reaction recipe
    actions. Pass a string describing the circumstances that caused the
    exceptional behavior.
    """
    pass

################################################################################

class GroupAtom(Vertex):
    """
    An atom group. This class is based on the :class:`Atom` class, except that
    it uses :ref:`atom types <atom-types>` instead of elements, and all
    attributes are lists rather than individual values. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `atomType`          ``list``            The allowed atom types (as :class:`AtomType` objects)
    `radicalElectrons`  ``list``            The allowed numbers of radical electrons (as short integers)
    `charge`            ``list``            The allowed formal charges (as short integers)
    `label`             ``str``             A string label that can be used to tag individual atoms
    `lonePairs`         ``list``            The number of lone electron pairs
    =================== =================== ====================================

    Each list represents a logical OR construct, i.e. an atom will match the
    group if it matches *any* item in the list. However, the
    `radicalElectrons`, and `charge` attributes are linked
    such that an atom must match values from the same index in each of these in
    order to match.
    """

    def __init__(self, atomType=None, radicalElectrons=None, charge=None, label='', lonePairs=None):
        Vertex.__init__(self)
        self.atomType = atomType or []
        for index in range(len(self.atomType)):
            if isinstance(self.atomType[index], str):
                self.atomType[index] = atomTypes[self.atomType[index]]
        self.radicalElectrons = radicalElectrons or []
        self.charge = charge or []
        self.label = label
        self.lonePairs = lonePairs or []

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        d = {
            'edges': self.edges,
            'connectivity1': self.connectivity1,
            'connectivity2': self.connectivity2,
            'connectivity3': self.connectivity3,
            'sortingLabel': self.sortingLabel,
        }
        atomType = self.atomType
        if atomType is not None:
            atomType = [a.label for a in atomType]
        return (GroupAtom, (atomType, self.radicalElectrons, self.charge, self.label, self.lonePairs), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling an object.
        """
        self.edges = d['edges']
        self.connectivity1 = d['connectivity1']
        self.connectivity2 = d['connectivity2']
        self.connectivity3 = d['connectivity3']
        self.sortingLabel = d['sortingLabel']

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return '[{0}]'.format(','.join([repr(a.label) for a in self.atomType]))

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "<GroupAtom {0!s}>".format(self)

    @property
    def bonds(self): return self.edges

    def copy(self):
        """
        Return a deep copy of the :class:`GroupAtom` object. Modifying the
        attributes of the copy will not affect the original.
        """
        return GroupAtom(self.atomType[:], self.radicalElectrons[:], self.charge[:], self.label)

    def __changeBond(self, order):
        """
        Update the atom group as a result of applying a CHANGE_BOND action,
        where `order` specifies whether the bond is incremented or decremented
        in bond order, and should be 1 or -1.
        """
        atomType = []
        for atom in self.atomType:
            if order == 1:
                atomType.extend(atom.incrementBond)
            elif order == -1:
                atomType.extend(atom.decrementBond)
            else:
                raise ActionError('Unable to update GroupAtom due to CHANGE_BOND action: Invalid order "{0}".'.format(order))
        if len(atomType) == 0:
            raise ActionError('Unable to update GroupAtom due to CHANGE_BOND action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def __formBond(self, order):
        """
        Update the atom group as a result of applying a FORM_BOND action,
        where `order` specifies the order of the forming bond, and should be
        'S' (since we only allow forming of single bonds).
        """
        if order != 'S':
            raise ActionError('Unable to update GroupAtom due to FORM_BOND action: Invalid order "{0}".'.format(order))
        atomType = []
        for atom in self.atomType:
            atomType.extend(atom.formBond)
        if len(atomType) == 0:
            raise ActionError('Unable to update GroupAtom due to FORM_BOND action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def __breakBond(self, order):
        """
        Update the atom group as a result of applying a BREAK_BOND action,
        where `order` specifies the order of the breaking bond, and should be
        'S' (since we only allow breaking of single bonds).
        """
        if order != 'S':
            raise ActionError('Unable to update GroupAtom due to BREAK_BOND action: Invalid order "{0}".'.format(order))
        atomType = []
        for atom in self.atomType:
            atomType.extend(atom.breakBond)
        if len(atomType) == 0:
            raise ActionError('Unable to update GroupAtom due to BREAK_BOND action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def __gainRadical(self, radical):
        """
        Update the atom group as a result of applying a GAIN_RADICAL action,
        where `radical` specifies the number of radical electrons to add.
        """
        radicalElectrons = []
        if any([len(atomType.incrementRadical) == 0 for atomType in self.atomType]):
            raise ActionError('Unable to update GroupAtom due to GAIN_RADICAL action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        for electron in self.radicalElectrons:
            radicalElectrons.append(electron + radical)
        # Set the new radical electron counts
        self.radicalElectrons = radicalElectrons

    def __loseRadical(self, radical):
        """
        Update the atom group as a result of applying a LOSE_RADICAL action,
        where `radical` specifies the number of radical electrons to remove.
        """
        radicalElectrons = []
        pairs = set()
        if any([len(atomType.decrementRadical) == 0 for atomType in self.atomType]):
            raise ActionError('Unable to update GroupAtom due to LOSE_RADICAL action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        for electron in self.radicalElectrons:
            electron = electron - radical
            if electron < 0:
                raise ActionError('Unable to update GroupAtom due to LOSE_RADICAL action: Invalid radical electron set "{0}".'.format(self.radicalElectrons))    
            radicalElectrons.append(electron)
            
        # Set the new radical electron counts
        self.radicalElectrons = radicalElectrons
        
    def __gainPair(self, pair):
        """
        Update the atom group as a result of applying a GAIN_PAIR action,
        where `pair` specifies the number of lone electron pairs to add.
        """
        lonePairs = []
        if any([len(atomType.incrementLonePair) == 0 for atomType in self.atomType]):
            raise ActionError('Unable to update GroupAtom due to GAIN_PAIR action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        for lonePairs in zip(self.lonePairs):
            lonePairs.append(lonePairs + pair)
        # Set the new lone electron pair count
        self.lonePairs = lonePairs
        
    def __losePair(self, pair):
        """
        Update the atom group as a result of applying a LOSE_PAIR action,
        where `pair` specifies the number of lone electron pairs to remove.
        """
        lonePairs = []
        if any([len(atomType.decrementLonePair) == 0 for atomType in self.atomType]):
            raise ActionError('Unable to update GroupAtom due to LOSE_PAIR action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        for lonePairs in zip(self.lonePairs):
            if lonePairs - pair < 0:
                raise ActionError('Unable to update GroupAtom due to LOSE_PAIR action: Invalid lone electron pairs set "{0}".'.format(self.lonePairs))
            lonePairs.append(lonePairs - pair)
        # Set the new lone electron pair count
        self.lonePairs = lonePairs

    def applyAction(self, action):
        """
        Update the atom group as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        act = action[0].upper()
        if act == 'CHANGE_BOND':
            self.__changeBond(action[2])
        elif act == 'FORM_BOND':
            self.__formBond(action[2])
        elif act == 'BREAK_BOND':
            self.__breakBond(action[2])
        elif act == 'GAIN_RADICAL':
            self.__gainRadical(action[2])
        elif act == 'LOSE_RADICAL':
            self.__loseRadical(action[2])
        elif action[0].upper() == 'GAIN_PAIR':
            self.__gainPair(action[2])
        elif action[0].upper() == 'LOSE_PAIR':
            self.__losePair(action[2])
        else:
            raise ActionError('Unable to update GroupAtom: Invalid action {0}".'.format(action))

    def equivalent(self, other):
        """
        Returns ``True`` if `other` is equivalent to `self` or ``False`` if not,
        where `other` can be either an :class:`Atom` or an :class:`GroupAtom`
        object. When comparing two :class:`GroupAtom` objects, this function
        respects wildcards, e.g. ``R!H`` is equivalent to ``C``.
        
        """
        cython.declare(group=GroupAtom)
        if not isinstance(other, GroupAtom):
            # Let the equivalent method of other handle it
            # We expect self to be an Atom object, but can't test for it here
            # because that would create an import cycle
            return other.equivalent(self)
        group=other
        
        cython.declare(atomType1=AtomType, atomtype2=AtomType, radical1=cython.short, radical2=cython.short,
                       lp1=cython.short, lp2=cython.short, charge1=cython.short, charge2=cython.short)
        # Compare two atom groups for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for atomType1 in self.atomType:
            for atomType2 in group.atomType:
                if atomType1.equivalent(atomType2): break
            else:
                return False
        for atomType1 in group.atomType:
            for atomType2 in self.atomType:
                if atomType1.equivalent(atomType2): break
            else:
                return False
        # Each free radical electron state in self must have an equivalent in other (and vice versa)
        for radical1 in self.radicalElectrons:
            if group.radicalElectrons:  # Only check if the list is non-empty.  An empty list indicates a wildcard.
                for radical2  in group.radicalElectrons:
                    if radical1 == radical2: break
                else:
                    return False
        for radical1 in group.radicalElectrons:
            if self.radicalElectrons:
                for radical2 in self.radicalElectrons:
                    if radical1 == radical2: break
                else:
                    return False
        for lp1 in self.lonePairs:
            if group.lonePairs:
                for lp2 in group.lonePairs:
                    if lp1 == lp2: break
                else:
                    return False
        #Each charge in self must have an equivalent in other (and vice versa)
        for charge1 in self.charge:
            if group.charge:
                for charge2 in group.charge:
                    if charge1 == charge2: break
                else:
                    return False
        for charge1 in group.charge:
            if self.charge:
                for charge2 in self.charge:
                    if charge1 == charge2: break
                else:
                    return False
        # Otherwise the two atom groups are equivalent
        return True

    def isSpecificCaseOf(self, other):
        """
        Returns ``True`` if `other` is the same as `self` or is a more
        specific case of `self`. Returns ``False`` if some of `self` is not
        included in `other` or they are mutually exclusive. 
        """
        cython.declare(group=GroupAtom)
        if not isinstance(other, GroupAtom):
            # Let the isSpecificCaseOf method of other handle it
            # We expect self to be an Atom object, but can't test for it here
            # because that would create an import cycle
            return other.isSpecificCaseOf(self)
        group=other
        
        cython.declare(atomType1=AtomType, atomtype2=AtomType, radical1=cython.short, radical2=cython.short, 
                       lp1=cython.short, lp2=cython.short, charge1=cython.short, charge2=cython.short)
        # Compare two atom groups for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for atomType1 in self.atomType: # all these must match
            for atomType2 in group.atomType: # can match any of these
                if atomType1.isSpecificCaseOf(atomType2): break
            else:
                return False
        # Each free radical electron state in self must have an equivalent in other (and vice versa)
        if self.radicalElectrons:
            for radical1 in self.radicalElectrons:
                if group.radicalElectrons:
                    for radical2 in group.radicalElectrons:
                        if radical1 == radical2: break
                    else:
                        return False
        else:
            if group.radicalElectrons: return False
        if self.lonePairs:
            for lp1 in self.lonePairs:
                if group.lonePairs:
                    for lp2 in group.lonePairs:
                        if lp1 == lp2: break
                    else:
                        return False
        else:
            if group.lonePairs: return False
        #Each charge in self must have an equivalent in other
        if self.charge:
            for charge1 in self.charge:
                if group.charge:
                    for charge2 in group.charge:
                        if charge1 == charge2: break
                    else:
                        return False
        else:
            if group.charge: return False
        # Otherwise self is in fact a specific case of other
        return True

    def isOxygen(self):
        """
        Return ``True`` if the atom represents an oxygen atom or ``False`` if
        not.
        """
        allOxygens = [atomTypes['O']] + atomTypes['O'].specific
        checkList=[x in allOxygens for x in self.atomType]

        return not False in checkList

    def isSulfur(self):
        """
        Return ``True`` if the atom represents an sulfur atom or ``False`` if
        not.
        """
        allSulfur = [atomTypes['S']] + atomTypes['S'].specific
        checkList=[x in allSulfur for x in self.atomType]

        return not False in checkList

    def hasWildcards(self):
        """
        Return ``True`` if the atom has wildcards in any of the attributes:
        atomtype, electronpairs, lone pairs, charge, and bond order. Returns
        ''False'' if no attribute has wildcards.
        """

        if len(self.atomType) > 1:
            return True
        elif len(self.radicalElectrons) > 1 or len(self.radicalElectrons) == 0:
            return True
        elif len(self.lonePairs) > 1:
            return True
        for bond in self.bonds.values():
            if len(bond.order) > 1:
                return True

        return False

    def makeSampleAtom(self):
        """

        Returns: a class :Atom: object analagous to the GroupAtom

        This makes a sample, so it takes the first element when there are multiple options inside of
        self.atomtype, self.radicalElectrons, self.lonePairs, and self.charge

        """

        #Use the first atomtype to determine element, even if there is more than one atomtype
        atomtype = self.atomType[0]
        element = None

        defaultLonePairs={'H': 0,
                          'D': 0,
                          'T': 0,
                          'He':1,
                          'C': 0,
                          'O': 2,
                          'N': 1,
                          'Si':0,
                          'S': 1,
                          'Ne':4,
                          'Cl':3,
                          'Ar':4,
        }

        for elementLabel in allElements:
            if atomtype is atomTypes[elementLabel] or atomtype in atomTypes[elementLabel].specific:
                element = elementLabel
                break
        else:
            #For types that correspond to more than one type of element, pick the first that appears in specific
            for subtype in atomtype.specific:
                if subtype.label in allElements:
                    element = subtype.label
                    break

        #dummy defaultAtom to get default values
        defaultAtom = mol.Atom()

        newAtom = mol.Atom(element = element,
                           radicalElectrons = self.radicalElectrons[0] if self.radicalElectrons else defaultAtom.radicalElectrons,
                           charge = self.charge[0] if self.charge else defaultAtom.charge,
                           lonePairs = self.lonePairs[0] if self.lonePairs else defaultAtom.lonePairs)

        #For some reason the default when no lone pairs is set to -100,
        #Based on git history, it is probably because RDKit requires a number instead of None
        #Instead we will set it to 0 here
        if newAtom.lonePairs == -100:
            newAtom.lonePairs = defaultLonePairs[newAtom.symbol]

        return newAtom
################################################################################

class GroupBond(Edge):
    """
    A bond group. This class is based on the :class:`Bond` class, except that
    all attributes are lists rather than individual values. The allowed bond
    types are given :ref:`here <bond-types>`. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `order`             ``list``            The allowed bond orders (as character strings)
    =================== =================== ====================================

    Each list represents a logical OR construct, i.e. a bond will match the
    group if it matches *any* item in the list.
    """

    def __init__(self, atom1, atom2, order=None):
        Edge.__init__(self, atom1, atom2)
        self.order = order or []

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return str(self.order)

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "<GroupBond {0!r}>".format(self.order)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (GroupBond, (self.vertex1, self.vertex2, self.order))

    def copy(self):
        """
        Return a deep copy of the :class:`GroupBond` object. Modifying the
        attributes of the copy will not affect the original.
        """
        return GroupBond(self.vertex1, self.vertex2, self.order[:])

    def isSingle(self):
        """
        Return ``True`` if the bond represents a single bond or ``False`` if
        not.
        """
        return self.order[0] == 'S' and len(self.order) == 1

    def isDouble(self):
        """
        Return ``True`` if the bond represents a double bond or ``False`` if
        not.
        """
        return self.order[0] == 'D' and len(self.order) == 1

    def isTriple(self):
        """
        Return ``True`` if the bond represents a triple bond or ``False`` if
        not.
        """
        return self.order[0] == 'T' and len(self.order) == 1

    def isBenzene(self):
        """
        Return ``True`` if the bond represents a benzene bond or ``False`` if
        not.
        """
        return self.order[0] == 'B' and len(self.order) == 1

    def __changeBond(self, order):
        """
        Update the bond group as a result of applying a CHANGE_BOND action,
        where `order` specifies whether the bond is incremented or decremented
        in bond order, and should be 1 or -1.
        """
        newOrder = []
        for bond in self.order:
            if order == 1:
                if bond == 'S':         newOrder.append('D')
                elif bond == 'D':       newOrder.append('T')
                else:
                    raise ActionError('Unable to update GroupBond due to CHANGE_BOND action: Invalid bond order "{0}" in set {1}".'.format(bond, self.order))
            elif order == -1:
                if bond == 'D':         newOrder.append('S')
                elif bond == 'T':       newOrder.append('D')
                else:
                    raise ActionError('Unable to update GroupBond due to CHANGE_BOND action: Invalid bond order "{0}" in set {1}".'.format(bond, self.order))
            else:
                raise ActionError('Unable to update GroupBond due to CHANGE_BOND action: Invalid order "{0}".'.format(order))
        # Set the new bond orders, removing any duplicates
        self.order = list(set(newOrder))

    def applyAction(self, action):
        """
        Update the bond group as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        if action[0].upper() == 'CHANGE_BOND':
            self.__changeBond(action[2])
        else:
            raise ActionError('Unable to update GroupBond: Invalid action {0}".'.format(action))

    def equivalent(self, other):
        """
        Returns ``True`` if `other` is equivalent to `self` or ``False`` if not,
        where `other` can be either an :class:`Bond` or an :class:`GroupBond`
        object.
        """
        cython.declare(gb=GroupBond)
        if not isinstance(other, GroupBond):
            # Let the equivalent method of other handle it
            # We expect self to be a Bond object, but can't test for it here
            # because that would create an import cycle
            return other.equivalent(self)
        gb = other
        
        cython.declare(order1=str, order2=str)
        # Compare two bond groups for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for order1 in self.order:
            for order2 in gb.order:
                if order1 == order2: break
            else:
                return False
        for order1 in gb.order:
            for order2 in self.order:
                if order1 == order2: break
            else:
                return False
        # Otherwise the two bond groups are equivalent
        return True

    def isSpecificCaseOf(self, other):
        """
        Returns ``True`` if `other` is the same as `self` or is a more
        specific case of `self`. Returns ``False`` if some of `self` is not
        included in `other` or they are mutually exclusive.
        """
        cython.declare(gb=GroupBond)
        if not isinstance(other, GroupBond):
            # Let the isSpecificCaseOf method of other handle it
            # We expect self to be a Bond object, but can't test for it here
            # because that would create an import cycle
            return other.isSpecificCaseOf(self)
        gb = other
        
        cython.declare(order1=str, order2=str)
        # Compare two bond groups for equivalence
        # Each atom type in self must have an equivalent in other
        for order1 in self.order: # all these must match
            for order2 in gb.order: # can match any of these
                if order1 == order2: break
            else:
                return False
        # Otherwise self is in fact a specific case of other
        return True

    def makeBond(self, molecule, atom1, atom2):
        """
        Creates a :class: Bond between atom1 and atom2 analogous to self

        The intended input arguments should be class :Atom: not class :GroupAtom:
        Args:
            atom1: First :class: Atom the bond connects
            atom2: Second :class: Atom the bond connects

        """
        newBond = mol.Bond(atom1, atom2, order = self.order[0])
        molecule.addBond(newBond)
################################################################################

class Group(Graph):
    """
    A representation of a molecular substructure group using a graph data
    type, extending the :class:`Graph` class. The `atoms` and `bonds` attributes
    are aliases for the `vertices` and `edges` attributes, and store 
    :class:`GroupAtom` and :class:`GroupBond` objects, respectively.
    Corresponding alias methods have also been provided.
    """

    def __init__(self, atoms=None, multiplicity=None):
        Graph.__init__(self, atoms)
        self.multiplicity = multiplicity if multiplicity else []
        self.update()

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Group, (self.vertices,))

    def __getAtoms(self): return self.vertices
    def __setAtoms(self, atoms): self.vertices = atoms
    atoms = property(__getAtoms, __setAtoms)

    def addAtom(self, atom):
        """
        Add an `atom` to the graph. The atom is initialized with no bonds.
        """
        return self.addVertex(atom)

    def addBond(self, bond):
        """
        Add a `bond` to the graph as an edge connecting the two atoms `atom1`
        and `atom2`.
        """
        return self.addEdge(bond)

    def getBonds(self, atom):
        """
        Return a list of the bonds involving the specified `atom`.
        """
        return self.getEdges(atom)

    def getBond(self, atom1, atom2):
        """
        Returns the bond connecting atoms `atom1` and `atom2`.
        """
        return self.getEdge(atom1, atom2)

    def hasAtom(self, atom):
        """
        Returns ``True`` if `atom` is an atom in the graph, or ``False`` if
        not.
        """
        return self.hasVertex(atom)

    def hasBond(self, atom1, atom2):
        """
        Returns ``True`` if atoms `atom1` and `atom2` are connected
        by an bond, or ``False`` if not.
        """
        return self.hasEdge(atom1, atom2)

    def removeAtom(self, atom):
        """
        Remove `atom` and all bonds associated with it from the graph. Does
        not remove atoms that no longer have any bonds as a result of this
        removal.
        """
        return self.removeVertex(atom)

    def removeBond(self, bond):
        """
        Remove the bond between atoms `atom1` and `atom2` from the graph.
        Does not remove atoms that no longer have any bonds as a result of
        this removal.
        """
        return self.removeEdge(bond)

    def sortAtoms(self):
        """
        Sort the atoms in the graph. This can make certain operations, e.g.
        the isomorphism functions, much more efficient.
        """
        return self.sortVertices()

    def copy(self, deep=False):
        """
        Create a copy of the current graph. If `deep` is ``True``, a deep copy
        is made: copies of the vertices and edges are used in the new graph.
        If `deep` is ``False`` or not specified, a shallow copy is made: the
        original vertices and edges are used in the new graph.
        """
        other = cython.declare(Group)
        g = Graph.copy(self, deep)
        other = Group(g.vertices)
        return other

    def update(self):

        self.updateConnectivityValues()
        self.updateFingerprint()


    def merge(self, other):
        """
        Merge two groups so as to store them in a single
        :class:`Group` object. The merged :class:`Group`
        object is returned.
        """
        g = Graph.merge(self, other)
        molecule = Group(atoms=g.vertices)
        return molecule

    def split(self):
        """
        Convert a single :class:`Group` object containing two or more
        unconnected groups into separate class:`Group` objects.
        """
        graphs = Graph.split(self)
        molecules = []
        for g in graphs:
            molecule = Group(atoms=g.vertices)
            molecules.append(molecule)
        return molecules

    def clearLabeledAtoms(self):
        """
        Remove the labels from all atoms in the molecular group.
        """
        for atom in self.vertices:
            atom.label = ''

    def containsLabeledAtom(self, label):
        """
        Return ``True`` if the group contains an atom with the label
        `label` and ``False`` otherwise.
        """
        for atom in self.vertices:
            if atom.label == label: return True
        return False

    def getLabeledAtom(self, label):
        """
        Return the atom in the group that is labeled with the given `label`.
        Raises :class:`ValueError` if no atom in the group has that label.
        """
        for atom in self.vertices:
            if atom.label == label: return atom
        raise ValueError('No atom in the functional group has the label "{0}".'.format(label))

    def getLabeledAtoms(self):
        """
        Return the labeled atoms as a ``dict`` with the keys being the labels
        and the values the atoms themselves. If two or more atoms have the
        same label, the value is converted to a list of these atoms.
        """
        labeled = {}
        for atom in self.vertices:
            if atom.label != '':
                if atom.label in labeled:
                    if isinstance(labeled[atom.label],list):
                        labeled[atom.label].append(atom)
                    else:
                        labeled[atom.label] = [labeled[atom.label]]
                        labeled[atom.label].append(atom)
                else:
                    labeled[atom.label] = atom
        return labeled

    def fromAdjacencyList(self, adjlist):
        """
        Convert a string adjacency list `adjlist` to a molecular structure.
        Skips the first line (assuming it's a label) unless `withLabel` is
        ``False``.
        """
        from .adjlist import fromAdjacencyList
        self.vertices, multiplicity = fromAdjacencyList(adjlist, group=True)
        if multiplicity is not None:
            self.multiplicity = multiplicity
        self.update()
        return self

    def toAdjacencyList(self, label=''):
        """
        Convert the molecular structure to a string adjacency list.
        """
        from .adjlist import toAdjacencyList
        return toAdjacencyList(self.vertices, multiplicity=self.multiplicity, label='', group=True)


    def updateFingerprint(self):
        """
        Update the molecular fingerprint used to accelerate the subgraph
        isomorphism checks.
        """
        cython.declare(atom=GroupAtom, atomType=AtomType)
        cython.declare(carbon=AtomType, nitrogen=AtomType, oxygen=AtomType, sulfur=AtomType)
        cython.declare(isCarbon=cython.bint, isNitrogen=cython.bint, isOxygen=cython.bint, isSulfur=cython.bint, radical=cython.int)
        
        carbon   = atomTypes['C']
        nitrogen = atomTypes['N']
        oxygen   = atomTypes['O']
        sulfur   = atomTypes['S']
        
        self.carbonCount   = 0
        self.nitrogenCount = 0
        self.oxygenCount   = 0
        self.sulfurCount   = 0
        self.radicalCount  = 0
        for atom in self.vertices:
            if len(atom.atomType) == 1:
                atomType   = atom.atomType[0]
                isCarbon   = atomType.equivalent(carbon)
                isNitrogen = atomType.equivalent(nitrogen)
                isOxygen   = atomType.equivalent(oxygen)
                isSulfur   = atomType.equivalent(sulfur)
                if isCarbon and not isNitrogen and not isOxygen and not isSulfur:
                    self.carbonCount += 1
                elif isNitrogen and not isCarbon and not isOxygen and not isSulfur:
                    self.nitrogenCount += 1
                elif isOxygen and not isCarbon and not isNitrogen and not isSulfur:
                    self.oxygenCount += 1
                elif isSulfur and not isCarbon and not isNitrogen and not isOxygen:
                    self.sulfurCount += 1
            if len(atom.radicalElectrons) == 1:
                radical = atom.radicalElectrons[0]
                self.radicalCount += radical

    def isIsomorphic(self, other, initialMap=None):
        """
        Returns ``True`` if two graphs are isomorphic and ``False``
        otherwise. The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Group` object, or a :class:`TypeError` is raised.
        """
        # It only makes sense to compare a Group to a Group for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        # Do the isomorphism comparison
        return Graph.isIsomorphic(self, other, initialMap)

    def findIsomorphism(self, other, initialMap=None):
        """
        Returns ``True`` if `other` is isomorphic and ``False``
        otherwise, and the matching mapping. The `initialMap` attribute can be
        used to specify a required mapping from `self` to `other` (i.e. the
        atoms of `self` are the keys, while the atoms of `other` are the
        values). The returned mapping also uses the atoms of `self` for the keys
        and the atoms of `other` for the values. The `other` parameter must
        be a :class:`Group` object, or a :class:`TypeError` is raised.
        """
        # It only makes sense to compare a Group to a Group for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        # Do the isomorphism comparison
        return Graph.findIsomorphism(self, other, initialMap)

    def isSubgraphIsomorphic(self, other, initialMap=None):
        """
        Returns ``True`` if `other` is subgraph isomorphic and ``False``
        otherwise. In other words, return ``True`` if self is more specific than other.
        The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Group` object, or a :class:`TypeError` is raised.
        """        
        cython.declare(group=Group)
        cython.declare(mult1=cython.short, mult2=cython.short)
        # It only makes sense to compare a Group to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        group = other
        
        if self.multiplicity:
            for mult1 in self.multiplicity:
                if group.multiplicity:
                    for mult2 in group.multiplicity:
                        if mult1 == mult2: break
                    else:
                        return False
        else:
            if group.multiplicity: return False
        # Do the isomorphism comparison
        return Graph.isSubgraphIsomorphic(self, other, initialMap)

    def findSubgraphIsomorphisms(self, other, initialMap=None):
        """
        Returns ``True`` if `other` is subgraph isomorphic and ``False``
        otherwise. In other words, return ``True`` is self is more specific than other.
        Also returns the lists all of valid mappings. The
        `initialMap` attribute can be used to specify a required mapping from
        `self` to `other` (i.e. the atoms of `self` are the keys, while the
        atoms of `other` are the values). The returned mappings also use the
        atoms of `self` for the keys and the atoms of `other` for the values.
        The `other` parameter must be a :class:`Group` object, or a
        :class:`TypeError` is raised.
        """
        cython.declare(group=Group)
        cython.declare(mult1=cython.short, mult2=cython.short)

        # It only makes sense to compare a Group to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        group = other
        
        if self.multiplicity:
            for mult1 in self.multiplicity:
                if group.multiplicity:
                    for mult2 in group.multiplicity:
                        if mult1 == mult2: break
                    else:
                        return []
        else:
            if group.multiplicity: return []
                
        # Do the isomorphism comparison
        return Graph.findSubgraphIsomorphisms(self, other, initialMap)
    
    def isIdentical(self, other):
        """
        Returns ``True`` if `other` is identical and ``False`` otherwise.
        The function `isIsomorphic` respects wildcards, while this function
        does not, make it more useful for checking groups to groups (as
        opposed to molecules to groups)
        """
        # It only makes sense to compare a Group to a Group for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        # An identical group is always a child of itself and 
        # is the only case where that is true. Therefore
        # if we do both directions of isSubgraphIsmorphic, we need
        # to get True twice for it to be identical
        if not self.isSubgraphIsomorphic(other):
            return False
        elif not other.isSubgraphIsomorphic(self):
            return False
        else:
            return True

    def isAromaticRing(self):
        """
        This method returns a boolean telling if the group has a 5 or 6 cyclic with
        benzene bonds exclusively
        """

        ring_size = len(self.atoms)
        if ring_size not in [5, 6]:
            return False
        for ringAtom in self.atoms:
            for bondedAtom, bond in ringAtom.edges.iteritems():
                if bondedAtom in self.atoms:
                    if not bond.isBenzene():
                        return False
        return True

    def standardizeAtomType(self):
        """
        This function changes the atomTypes in a group if the atom must
        be a specific atomType based on its bonds and valency.

        Currently only standardizes oxygen, carbon and sulfur atomTypes

        We also only check when there is exactly one atomType,
        one bondType, one radical setting.
        For any group where there are wildcards or multiple attributes,
        we cannot apply this check.

        In the case where the atomType is ambigious based on bonds
        and valency, this function will not change the type.

        Returns a 'True' if the group was modified otherwise returns 'False'
        """

        modified = False

        #If this atom or any of its ligands has wild cards, then don't try to standardize
        if self.hasWildCards: return modified
        for bond12, atom2 in self.bonds.iteritems():
            if atom2.hasWildCards: return modified

        #dictionary of element to expected valency
        valency = {atomTypes['C'] : 4,
                   atomTypes['O'] : 2,
                   atomTypes ['S']: 2,
                   atomTypes['Si']: 4
                   }

        #list of :class:AtomType which are elements with more sub-divided atomtypes beneath them
        specifics= [elementLabel for elementLabel in allElements if not elementLabel in nonSpecifics]
        for index, atom in enumerate(self.atoms):
            claimedAtomType = atom.atomType[0]
            newAtomType = None
            element = None
            #Ignore elements that do not have more than one atomtype
            if claimedAtomType.label in nonSpecifics: continue
            for elementLabel in specifics:
                if claimedAtomType.label == elementLabel or atomTypes[claimedAtomType.label] in atomTypes[elementLabel].specific:
                    element = atomTypes[elementLabel]
                    break

            #claimedAtomType is not in one of the specified elements
            if not element: continue
            #Don't standardize atomtypes for nitrogen for now. My feeling is that
            # the work on the nitrogen atomtypes is still incomplete
            elif element == atomTypes['N']: continue

            groupFeatures = getFeatures(atom, atom.bonds)
            # print groupFeatures

            single = groupFeatures[0]
            allDouble = groupFeatures[1]
            triple = groupFeatures[5]
            benzene = groupFeatures[6]

            bondValency = single + 2 * allDouble + 3 * triple + 4.0/3.0 * benzene
            filledValency =  atom.radicalElectrons[0] + bondValency
            # print index, atom, filledValency

            #For an atomtype to be known for certain, the valency must be filled
            #within 1 of the total valency available
            if filledValency >= valency[element] - 1:
                for specificAtomType in element.specific:
                    atomtypeFeatureList = specificAtomType.getFeatures()
                    for molFeature, atomtypeFeature in zip(groupFeatures, atomtypeFeatureList):
                        if atomtypeFeature == []:
                            continue
                        elif not molFeature in atomtypeFeature:
                            break
                    else:
                        if specificAtomType is atomTypes['Oa'] or specificAtomType is atomTypes['Sa']:
                            if atom.lonePairs == 3 or atom.radicalElectrons == 2:
                                newAtomType = specificAtomType
                                break
                        else:
                            newAtomType = specificAtomType
                            break

            #set the new atom type if the algorithm found one
            if newAtomType and not newAtomType is claimedAtomType:
                atom.atomType[0] = newAtomType
                modified = True

        return modified

    def addExplicitLigands(self):
        """
        This function Od/Sd ligand to CO or CS atomtypes if they are not already there.

        Returns a 'True' if the group was modified otherwise returns 'False'
        """

        modified = False

        atomsToAddTo=[]

        for index, atom in enumerate(self.atoms):
            claimedAtomType = atom.atomType[0]
            #Do not perform is this atom has wildCards
            if atom.hasWildCards: continue
            elif claimedAtomType is atomTypes['CO'] or claimedAtomType is atomTypes['CS']:
                for atom2, bond12 in atom.bonds.iteritems():
                    if bond12.isDouble():
                        break
                else: atomsToAddTo.append(index)

        for atomIndex in atomsToAddTo:
            modified = True
            if self.atoms[atomIndex].atomType[0] is atomTypes['CO']:
                newAtom = GroupAtom(atomType=[atomTypes['Od']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
            elif self.atoms[atomIndex].atomType[0] is atomTypes['CS']:
                newAtom = GroupAtom(atomType=[atomTypes['Sd']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
            self.addAtom(newAtom)
            newBond = GroupBond(self.atoms[atomIndex], newAtom, order=['D'])
            self.addBond(newBond)

        return modified

    def standardizeGroup(self):
        """
        This function modifies groups to make them have a standard AdjList form.

        Currently it makes atomtypes as specific as possible and makes CO/CS atomtypes
        have explicit Od/Sd ligands. Other functions can be added as necessary

        Returns a 'True' if the group was modified otherwise returns 'False'
        """

        modified = False

        #If viable then we apply current conventions:
        checkList=[]
        checkList.append(self.standardizeAtomType())
        checkList.append(self.addExplicitLigands())

        return True in checkList

    def addImplicitAtomsFromAtomType(self):
        """

        Returns: a modified group with implicit atoms added
        Add implicit double/triple bonded atoms O, S or R, for which we will use a C

        Not designed to work with wildcards
        """

        #dictionary of implicit atoms and their bonds
        implicitAtoms = {}
        lonePairsRequired = {}

        copyGroup = deepcopy(self)

        for atom1 in copyGroup.atoms:
            atomtypeFeatureList = atom1.atomType[0].getFeatures()
            lonePairsRequired[atom1]=atomtypeFeatureList[7]

            #set to 0 required if empty list
            atomtypeFeatureList = [featureList if featureList else [0] for featureList in atomtypeFeatureList]
            allDoubleRequired = atomtypeFeatureList[1]
            rDoubleRequired = atomtypeFeatureList[2]
            oDoubleRequired = atomtypeFeatureList[3]
            sDoubleRequired = atomtypeFeatureList[4]
            tripleRequired = atomtypeFeatureList[5]

            #count up number of bonds
            single = 0; rDouble = 0; oDouble = 0; sDouble = 0; triple = 0; benzene = 0
            for atom2, bond12 in atom1.bonds.iteritems():
                # Count numbers of each higher-order bond type
                if bond12.isSingle():
                    single += 1
                elif bond12.isDouble():
                    if atom2.isOxygen():
                        oDouble += 1
                    elif atom2.isSulfur():
                        sDouble += 1
                    else:
                        # rDouble is for double bonds NOT to oxygen or Sulfur
                        rDouble += 1
                elif bond12.isTriple(): triple += 1
                elif bond12.isBenzene(): benzene += 1


            while oDouble < oDoubleRequired[0]:
                oDouble +=1
                newAtom = GroupAtom(atomType=[atomTypes['O']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
                newBond = GroupBond(atom1, newAtom, order=['D'])
                implicitAtoms[newAtom] = newBond
            while sDouble < sDoubleRequired[0]:
                sDouble +=1
                newAtom = GroupAtom(atomType=[atomTypes['S']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
                newBond = GroupBond(atom1, newAtom, order=['D'])
                implicitAtoms[newAtom] = newBond
            while rDouble < rDoubleRequired[0] or rDouble + oDouble + sDouble < allDoubleRequired[0]:
                rDouble +=1
                newAtom = GroupAtom(atomType=[atomTypes['C']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
                newBond = GroupBond(atom1, newAtom, order=['D'])
                implicitAtoms[newAtom] = newBond
            while triple < tripleRequired[0]:
                triple +=1
                newAtom = GroupAtom(atomType=[atomTypes['C']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
                newBond = GroupBond(atom1, newAtom, order=['T'])
                implicitAtoms[newAtom] = newBond

        for atom, bond in implicitAtoms.iteritems():
            copyGroup.addAtom(atom)
            copyGroup.addBond(bond)

        for atom, lonePair in lonePairsRequired.iteritems():
            if lonePair: atom.lonePairs = lonePair

        return copyGroup

    def addImplicitBenzene(self):
        """
        Returns: A modified group with any implicit benzene rings added

        This method currently does not if there are wildcards in atomtypes or bond orders
        The current algorithm also requires that all Cb and Cbf are atomtyped

        There are other cases where the algorithm doesn't work. For example whenever there
        are many dangling Cb or Cbf atoms not in a ring, it is likely fail. In the database test
        (the only use thus far), we will require that any group with more than 3 Cbfs have
        complete rings. This is much stricter than this method can handle, but right now
        this method cannot handle very general cases, so it is better to be conservative. 
        """

        def classifyBenzeneCarbons(group):
            """
            Args:
                group: :class:Group with atoms to classify

            Returns: tuple with lists of each atom classification
            """
            cbAtomList = []
            cbfAtomList = [] #All Cbf Atoms
            cbfAtomList1 = [] #Cbf Atoms that are bonded to exactly one other Cbf (part of 2 rings)
            cbfAtomList2 = [] #Cbf that are sandwiched between two other Cbf (part of 2 rings)
            cbfAtomList3 = [] #Cbf atoms that are sandwiched between three other Cbf atoms (part of 3 rings)
            connectedCbfs={} #dictionary of connections to other cbfAtoms

            #Only want to work with benzene bonds on carbon
            labelsOfCarbonAtomTypes = [x.label for x in atomTypes['C'].specific] + ['C']

            for atom in group.atoms:
                if not atom.atomType[0].label in labelsOfCarbonAtomTypes: continue
                elif atom.atomType[0].label == 'Cb':
                    cbAtomList.append(atom)
                elif atom.atomType[0].label == 'Cbf':
                    cbfAtomList.append(atom)
                else:
                    benzeneBonds = 0
                    for atom2, bond12 in atom.bonds.iteritems():
                        if bond12.isBenzene(): benzeneBonds+=1
                    if benzeneBonds > 2: cbfAtomList.append(atom)
                    elif benzeneBonds >0: cbAtomList.append(atom)

            #further sort the cbf atoms
            for cbfAtom in cbfAtomList:
                fbBonds = 0
                connectedCbfs[cbfAtom] = []
                for atom2, bond in cbfAtom.bonds.iteritems():
                    if bond.order[0] == 'B' and atom2 in cbfAtomList:
                        fbBonds +=1
                        connectedCbfs[cbfAtom].append(atom2)
                if fbBonds < 2: cbfAtomList1.append(cbfAtom)
                elif fbBonds == 2: cbfAtomList2.append(cbfAtom)
                elif fbBonds == 3: cbfAtomList3.append(cbfAtom)

            #check that cbfAtoms only have benzene bonds
            for cbfAtom in cbfAtomList:
                for atom2, bond12 in cbfAtom.bonds.iteritems():
                    assert bond12.isBenzene(), "Cbf atom in {0} has a bond with an order other than 'B'".format(group)

            return (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, cbfAtomList3, connectedCbfs)

        def addBenzeneAtomToGroup(group, connectingAtom, cbf=False):
            """
            This function adds a benzene atom to the group

            Args:
                group: :class:Group that we want to add a benzene atom to
                connectingAtom: :class:GroupAtom that is connected to the new benzene atom
                cbf: boolean indicating if the new atom should be a cbf
            Returns: a tuple containing the modified group and the new atom

            """
            newAtom = None
            if cbf:
                newAtom = GroupAtom(atomType=[atomTypes['Cbf']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
            else:
                newAtom = GroupAtom(atomType=[atomTypes['Cb']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
            newBond = GroupBond(connectingAtom, newAtom, order=['B'])
            group.addAtom(newAtom)
            group.addBond(newBond)

            return (group, newAtom)

        def sortByConnectivity(atomList):
            """
            Args:
                atomList: input list of atoms

            Returns: a sorted list of atoms where each atom is connected to a previous
            atom in the list if possible
            """
            sortedAtomList=[]
            # print sortedAtomList, type(sortedAtomList)
            while atomList:
                sortedAtomList.append(atomList.pop(0))
                previousLength = len(sortedAtomList)
                for atom1 in atomList:
                    added = False
                    for atom2, bond12 in atom1.bonds.iteritems():
                        if atom2 in sortedAtomList:
                            sortedAtomList.append(atom1)
                            added = True
                            break
                    if added : continue
                for atom1 in sortedAtomList[previousLength:]:
                    atomList.remove(atom1)

            return sortedAtomList

        def checkSet(superList, subList):
            """
            Args:
                superList: list to check if superset of partList
                subList:  list to check if subset of superList

            Returns: Boolean to see if superList is a superset of subList

            """
            superSet = set(superList)
            subSet = set(subList)
            return superSet.issuperset(subSet)

        def mergeOverlappingBenzeneRings(ring1, ring2, od):
            """
            The input arguements of rings are always in the order that the atoms appear
            inside the ring. That is, each atom is connected to the ones adjacent on the
            list.

            Args:
                ring1: list of :class:GroupAtoms representing first partial ring to merge
                ring2: list :class:GroupAtoms representing second partial ring to merge
                od: in for overlap distance

            This function tries to see if the beginning or ends of each list have the
            same atom objects, i.e the two part rings should be merged together.

            Returns: If rings are mergable, returns a new list of the merged ring, otherwise
            an empty list

            """
            newRing = []
            #ring already complete
            if len(ring1) ==6 or len(ring2) == 6: return newRing

            #start of ring1 matches end of ring2
            matchList1 = [x1 is x2 for x1,x2 in zip(ring1[-od:],ring2[:od])]
            #end of ring1 matches end of ring2
            matchList2 = [x1 is x2 for x1,x2 in zip(ring1[-od:],ring2[:od-1:-1])]
            #start of ring1 matches end of ring2
            matchList3 = [x1 is x2 for x1,x2 in zip(ring1[:od],ring2[-od:])]
            #start of ring1 matches start of ring2
            matchList4 = [x1 is x2 for x1,x2 in zip(ring1[:od],ring2[od::-1])]
            if not False in matchList1:
                newRing = ring1 +ring2[od:]
            elif not False in matchList2:
                newRing = ring1 + ring2[-od-1::-1]
            elif not False in matchList3:
                newRing = ring2[:-od] + ring1
            elif not False in matchList4:
                newRing = ring2[:od-1:-1] + ring1

            return newRing

        def addCbAtomToRing(ring, cbAtom):
            """
            Every 'Cb' atom belongs in exactly one benzene ring. This function checks
            adds the cbAtom to the ring (in connectivity order) if the cbAtom is connected
            to any the last or first atom in the partial ring.

            Args:
                ring: list of :class:GroupAtoms representing a partial ring to merge
                cbAtom: :class:GroupAtom with atomtype 'Cb'

            Returns: If cbAtom connects to the beginning or end of ring, returns a
            new list of the merged ring, otherwise an empty list

            """

            mergedRing = []
            #ring already complete
            if len(ring) == 6 : return mergedRing
            for atom2, bond12 in cbAtom.bonds.iteritems():
                if bond12.isBenzene():
                    if atom2 is ring[-1]: mergedRing = ring+[cbAtom]
                    elif atom2 is ring[0]: mergedRing = [cbAtom] +ring

            return mergedRing

        #######################################################################################
        #start of main algorithm
        copyGroup = deepcopy(self)
        """
        Step 1. Classify all atoms as Cb, Cbf1, Cbf2, Cbf3, ignoring all non-benzene carbons

        Every carbon atom in a benzene ring can be defined as one of the following:
        Cb - benzene carbon in exclusively one ring. Can have one single bond
        Cbf - general classification for any benzene carbon that connects two fused benzene rings
        Cbf1 - Cbf that is bonded to exactly one other Cbf, exclusively in two different benzene rings
        Cbf2 - Cbf that is bonded to exactly two other Cbfs, exclusively in two different benzene ring
        Cbf3 - Cbf that is bonded to exactly three other Cbfs, exclusively in three different benzene rings

        The dictionary connectedCbfs has a cbf atom as key and the other Cbf atoms it is connected to as values
        """
        (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, cbfAtomList3, connectedCbfs) = classifyBenzeneCarbons(copyGroup)

        print "cbf1", cbfAtomList1
        print "cbf2", cbfAtomList2
        print "cbf3", cbfAtomList3

        """
        #Step 2. Partner up each Cbf1 and Cbf2 atom

        For any fused benzene rings, there will always be exactly two Cbf atoms that join the benzne
        rings. Therefore, we can say that every Cbf1 atom has one 'partner' Cbf atom in which it
        share lies in two benzene rings with. If you try to draw a couple example PAHs, you'll
        find that Cbf2 atoms also have one exclusive 'partner', while Cbf3 atoms are actually
        'partnered' with every atom it is bonded to.

        In this step, we attempt to find the partner for each Cbf1 and Cbf2 atom to make creating
        the rings easier later. We find (or create) partners for every Cbf1 atom first because
        they only have one option for a partner, by definition. Cbf2 atoms can be trickier to
        partner up correctly. However because some are partnered with Cbf1 atoms, we can
        eventually figure out correct Cbf2/Cbf2 pairs by eliminating the ones partnered
        with Cbf1 atoms.
        """
        partners = {} #dictionary of exclusive partners, has 1:2 and 2:1
        for cbfAtom in cbfAtomList1:
            if cbfAtom in partners: continue
            #if cbfAtom has a connected cbf it must be the partner
            elif connectedCbfs[cbfAtom] and connectedCbfs[cbfAtom][0] not in partners:
                partners[cbfAtom] = connectedCbfs[cbfAtom][0]
                if connectedCbfs[cbfAtom][0] not in cbfAtomList3:
                    partners[connectedCbfs[cbfAtom][0]] = cbfAtom
            else:
                potentialPartner = None
                for atom2, bond12 in cbfAtom.bonds.iteritems():
                    if atom2 in partners: continue
                    elif bond12.isBenzene():
                        hasSingle = [True if bond23.isSingle() else False for bond23 in atom2.bonds.values()]
                        # print "potentialPartner", atom2, hasSingle
                        if not True in hasSingle:
                            potentialPartner= atom2
                #Make a Cb atom the partner, now marking it as a Cbfatom
                if potentialPartner:
                    partners[cbfAtom] = potentialPartner
                    partners[potentialPartner]= cbfAtom
                #otherwise create a new atom to be the partner
                else:
                    print "make new group at cbf1 partner step"
                    (copyGroup, newAtom) = addBenzeneAtomToGroup(copyGroup, cbfAtom, True)
                    partners[cbfAtom] = newAtom
                    partners[newAtom] = cbfAtom

        #reclassify all atoms since we may have added new ones
        (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, cbfAtomList3, connectedCbfs) = classifyBenzeneCarbons(copyGroup)

        print "group after partners cb1", copyGroup.atoms
        print "partners after cb1", partners
        print "cbf after partners cb1", cbfAtomList
        print "cbf1 after partners cb1", cbfAtomList1
        print "cbf2 after partners cb1", cbfAtomList2
        print "cbf3 after partners cb1", cbfAtomList3

        #Partner up leftover cbfAtom2
        cbf2CheckList = [True if cbfAtom in partners else False for cbfAtom in cbfAtomList2]
        while False in cbf2CheckList:
            for cbfAtom in cbfAtomList2:
                if cbfAtom in partners: continue
                cbf2Connections = cbfAtom.bonds.keys()
                for atom2, bond12 in cbfAtom.bonds.iteritems():
                    #if connected to a cbf3, it must be the partner
                    if cbf2Connections[0] in cbfAtomList3:
                        partners[cbfAtom] = cbf2Connections[0]
                        break
                    #otherwise we have found the partner if atom2 is not in partners
                    #and the other atom is in partners
                    elif not atom2 not in partners:
                        #get the other connected atom
                        other = None
                        for atom3 in connectedCbfs[cbfAtom]:
                            if atom3 is not atom2: other = atom3
                        if other in partners:
                            partners[cbfAtom] = atom2
                            partners[atom2] = cbfAtom
                            break
                        #if we are not sure which atom is the partner, we can skip it for now
                        #and we will eventually figure it out as we knock out other options
            cbf2CheckList = [True if cbfAtom in partners else False for cbfAtom in cbfAtomList2]

        print "partners after cbf2", partners
        print "cbf after partners cb2", cbfAtomList
        print "cbf1 after partners cbf2", cbfAtomList1
        print "cbf2 after partners cbf2", cbfAtomList2
        print "cbf3 after partners cbf2", cbfAtomList3

        #debug lines everything should have a partner now
        for cbfAtom in cbfAtomList1+cbfAtomList2:
            # print cbfAtom
            assert(cbfAtom in partners)
            if not partners[cbfAtom] in cbfAtomList3:
                assert cbfAtom is partners[partners[cbfAtom]]

        #I think this is now redundat....but not sure
        #It is impossible to have more an odd number of cbf atoms
        if len(cbfAtomList)%2 == 1: raise Exception("Should be even number of Cbf atoms at this point")

        """
        Step 3. Sort all lists by connectivity

        In the coming steps, we will sort Cb/Cbf atom into their benzene rings. If we cannot
        find a ring to sort an atom into, we will create a new ring containing that atom.
        It is important that we always check atoms that are already connected to existing rings
        before completely disconnected atoms. Otherwise, we will erroneously create new rings.
        """
        cbAtomList = sortByConnectivity(cbAtomList)
        cbfAtomList1 = sortByConnectivity(cbfAtomList1)
        cbfAtomList2 = sortByConnectivity(cbfAtomList2)
        cbfAtomList3 = sortByConnectivity(cbfAtomList3)

        """
        Step 4. Initalize the list of rings with any benzene rings that are already explicitly stated

        The variable rings is a list of lists. Each list in rings represents one full benzene rings,
        so it will eventually have six benzene carbons in it. Each ring's list will have the atoms
        sorted by connectivity, such that any atom is bonded to the atoms preceding and following it
        in the list. The first and last atom of the list will also be bonded together.
        """
        rings=[cycle for cycle in copyGroup.getAllCyclesOfSize(6) if Group(atoms = cycle).isAromaticRing()]

        print "rings when first initialized", len(rings), rings

        """
        Step 5. Add Cbf3 atoms to the correct rings

        Every Cbf3 atom is part of three different rings. These three rings can be defined by
        three-membered combinations of ligand1-Cbf3atom-ligand2, where the ligands are any
        unique atom connected to the Cbf3 atom. In this step, we start with 'ring seeds' of
        the three-membered combinations and try to merge them into existing rings. If no
        suitable ring is found, we create a new ring from the ring seed.
        """
        for cbfAtom in cbfAtomList3:
            ringsToAdd = []
            ligands = cbfAtom.bonds.keys() #should always have 3 ligands
            #these cannot be duplicate rings since the three in newRingSeeds we test are unique
            newRingSeeds = [[ligands[0], cbfAtom, ligands[1]],
                            [ligands[0], cbfAtom, ligands[2]],
                            [ligands[1], cbfAtom, ligands[2]]]
            #Check for duplicates, merge or create new rings
            # rings = addRingseedsToRings(rings, newRingSeeds, 2)
            for index1, ring1 in enumerate(newRingSeeds):
                mergeRingDict={}
                for index, ring2 in enumerate(rings):
                    #check if a duplicate of a fully created ring
                    if checkSet(ring2, ring1): break
                    #Next try to merge the ringseed into rings
                    mergeRing = mergeOverlappingBenzeneRings(ring2, ring1, 2)
                    if mergeRing:
                        print index1, index
                        mergeRingDict[index] = mergeRing
                        break
                #otherwise add this ringSeed because it represents a completely new ring
                else: rings.append(ring1)
                #if we merged a ring, we need to remove the old ring from rings and add the merged ring
                rings = [rings[index] if not index in mergeRingDict else mergeRingDict[index] for index in range(len(rings))]
                print "rings after iteration of cbf add ringSeeds", rings


        print "rings after initial cbf3 handle", len(rings), rings
        """
        Step 6. Add Cbf2 atoms to the correct rings

        Every Cbf2 atom is in two rings with unique ring seeds defined by
        partneredCbf-Cbf2atom-otherCbf and partneredCbf-Cbf2Atom-CbAtom. This step is
        almost analogous to the last step, except we may need to create the Cbatom
        in the last ring seed if it is not available.
        """
        for cbfAtom in cbfAtomList2:
            if connectedCbfs[cbfAtom][0] is partners[cbfAtom]: otherCbf = connectedCbfs[cbfAtom][1]
            else: otherCbf = connectedCbfs[cbfAtom][0]
            #These two ring seeds represent the two unique rings
            newRingSeeds = [[partners[cbfAtom], cbfAtom, otherCbf],
                            [partners[cbfAtom], cbfAtom]]
            allLigands = cbfAtom.bonds.keys()
            #add a new cb atom to the second newRing seed
            if len(allLigands) == 2:
                (copyGroup, newAtom) = addBenzeneAtomToGroup(copyGroup, cbfAtom)
                newRingSeeds[1].append(newAtom)
            #join the existing atom to the ringSeed
            elif len(allLigands) == 3:
                for atom2 in allLigands:
                    if atom2 not in connectedCbfs[cbfAtom]:
                        newRingSeeds[1].append(atom2)
                        break
            print "new ringSeeds in cbf2", newRingSeeds
            #Check for duplicates, merge or create new rings
            # rings = addRingseedsToRings(rings, newRingSeeds, 2)
            for index1, ring1 in enumerate(newRingSeeds):
                mergeRingDict={}
                for index, ring2 in enumerate(rings):
                    #check if a duplicate of a fully created ring
                    if checkSet(ring2, ring1): break
                    #Next try to merge the ringseed into rings
                    mergeRing = mergeOverlappingBenzeneRings(ring2, ring1, 2)
                    if mergeRing:
                        print index1, index
                        mergeRingDict[index] = mergeRing
                        break
                #otherwise add this ringSeed because it represents a completely new ring
                else: rings.append(ring1)
                #if we merged a ring, we need to remove the old ring from rings and add the merged ring
                rings = [rings[index] if not index in mergeRingDict else mergeRingDict[index] for index in range(len(rings))]
                print "rings after iteration of cbf add ringSeeds", rings

        print "rings after cbf2 handle", len(rings), rings
        """
        Step 7. Add Cbf1 atoms to the correct rings

        Every Cbf1 atom is in two rings with its partner. In this step, we add this ring seed
        twice to the rings.
        """
        print "new cbf1 atoms", cbfAtomList1
        for cbfAtom in cbfAtomList1:
            newRingSeed = [partners[cbfAtom], cbfAtom]
            inRing = 0
            #check to see if duplicate of an existing ring
            for ring in rings:
                print ring
                if checkSet(ring, newRingSeed):
                    inRing +=1
                    print "already in Ring", inRing
            #move on to next cbfAtom if we found two rings
            if inRing ==2: continue
            #try to merge into existing rings, if cbf1 is connected
            for index, ring in enumerate(rings):
                mergeRingDict={}
                mergeRing = mergeOverlappingBenzeneRings(ring, newRingSeed, 1)
                if mergeRing:
                    inRing+=1
                    print "merged into Ring", inRing
                    mergeRingDict[index]=mergeRing
                    #move on to next cbfAtom if we found two rings
                    if inRing ==2: break
            #if we merged a ring, we need to remove the old ring from rings and add the merged ring
            rings = [rings[index] if not index in mergeRingDict else mergeRingDict[index] for index in range(len(rings))]
            #if we still dont have two ring, we create a completely new ring
            if inRing < 2:
                for x in range(2-inRing):
                    print "add new Ring", inRing
                    rings.append(copy(newRingSeed))

        print "rings after cbf1 handle", len(rings), rings

        """
        Step 8. Add Cb atoms to the correct rings

        Each Cb atom is part of exactly one benzene ring. In this step we merge or make
        a new ring for each Cb atom.
        """
        for cbAtom in cbAtomList:
            inRing = 0
            #check to see if already in a ring
            for ring in rings:
                if checkSet(ring, [cbAtom]): inRing +=1
            #move on to next ring cbAtom if in a ring
            if inRing == 1 : continue
            #check to see if can be merged to an existing ring
            for index, ring in enumerate(rings):
                mergeRingDict={}
                mergeRing = addCbAtomToRing(ring, cbAtom)
                if mergeRing:
                    inRing+=1
                    mergeRingDict[index]=mergeRing
                    break
            #if we merged a ring, we need to remove the old ring from rings and add the merged ring
            rings = [rings[index] if not index in mergeRingDict else mergeRingDict[index] for index in range(len(rings))]
            #Start completely new ring if not any of above true
            if inRing == 0:
                rings.append([cbAtom])

        print "rings after cb handle", len(rings), rings

        """
        Step 9. Grow each partial ring up to six carbon atoms

        In this step we create new Cb atoms and add them to any rings which do not have 6 atoms.
        """
        mergedRingDict={}
        for index, ring in enumerate(rings):
            carbonsToGrow = 6-len(ring)
            mergedRingDict[index] = []
            for x in range(carbonsToGrow):
                if x ==0: lastAtom = ring[-1]
                else: lastAtom = mergedRingDict[index][-1]
                #add a new atom to the ring and the group
                (copyGroup, newAtom) = addBenzeneAtomToGroup(copyGroup, lastAtom)
                mergedRingDict[index].append(newAtom)
                #At the end attach to the other endpoint
                if x == carbonsToGrow -1:
                    newBond = GroupBond(ring[0], newAtom, order=['B'])
                    copyGroup.addBond(newBond)

        #Replace elements in rings for debugging
        for index in range(len(rings)):
            rings[index].extend(mergedRingDict[index])
        print "grown rings", len(rings), rings

        return copyGroup

    def makeSampleMolecule(self):
        """
        Returns: A sample class :Molecule: from the group
        """

        # print self
        # print self.atoms

        #Add implicit atoms
        modifiedGroup = self.addImplicitAtomsFromAtomType()
        #Make dictionary of :GroupAtoms: to :Atoms:
        atomDict = {}
        for atom in modifiedGroup.atoms:
            atomDict[atom] = atom.makeSampleAtom()

        #create the molecule
        newMolecule = mol.Molecule(atoms = atomDict.values())
        #Add explicit bonds to :Atoms:
        for atom1 in modifiedGroup.atoms:
            for atom2, bond12 in atom1.bonds.iteritems():
                bond12.makeBond(newMolecule, atomDict[atom1], atomDict[atom2])


        #Saturate up to expected valency
        for atom in newMolecule.atoms:
            statedCharge = atom.charge
            atom.updateCharge()
            if atom.charge - statedCharge:
                hydrogenNeeded = atom.charge - statedCharge
                for x in range(hydrogenNeeded):
                    newH = mol.Atom('H', radicalElectrons=0, lonePairs=0, charge=0)
                    newBond = mol.Bond(atom, newH, 'S')
                    newMolecule.addAtom(newH)
                    newMolecule.addBond(newBond)
                atom.updateCharge()
            # print statedCharge, atom.charge, type(atom.charge)

        newMolecule.update()
        # print newMolecule.atoms

        return newMolecule