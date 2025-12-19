#!/usr/bin/env python2.3
#
#
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/AutoDockTools/MoleculePreparation.py,v 1.89.2.2 2016/02/11 09:24:07 annao Exp $ 
#
#
#
import os, string
import numpy, math 

from scripts.mgltools.MolKit import Read
from scripts.mgltools.MolKit.torTree import TorTree
from scripts.mgltools.MolKit.pdbWriter import PdbqWriter, PdbqsWriter, PdbqtWriter
import scripts.mgltools.MolKit.molecule
import scripts.mgltools.MolKit.protein
from scripts.mgltools.MolKit.protein import ResidueSet
from scripts.mgltools.MolKit.mol2Parser import Mol2Parser

from scripts.mgltools.MolKit.molecule import AtomSet, Bond, BondSet
from scripts.mgltools.MolKit.chargeCalculator import GasteigerChargeCalculator
from scripts.mgltools.MolKit.chargeCalculator import KollmanChargeCalculator
from scripts.mgltools.MolKit.hydrogenBuilder import HydrogenBuilder
#import the Kollman charges dictionary to get keys it recognizes...
from scripts.mgltools.MolKit.qkollua import q
from scripts.mgltools.MolKit.bondSelector import RotatableBondSelector, AmideBondSelector, LeafBondSelector
from scripts.mgltools.MolKit.bondSelector import GuanidiniumBondSelector

from scripts.mgltools.AutoDockTools.atomTypeTools import NonpolarHydrogenMerger
from scripts.mgltools.AutoDockTools.atomTypeTools import LonepairMerger
from scripts.mgltools.AutoDockTools.atomTypeTools import SolvationParameterizer
from scripts.mgltools.AutoDockTools.atomTypeTools import AromaticCarbonManager
from scripts.mgltools.AutoDockTools.atomTypeTools import AutoDock4_AtomTyper

from scripts.mgltools.AutoDockTools.AutoDockBondClassifier import AutoDockBondClassifier

# for LigandRandomizer
from scripts.mgltools.mglutil.math.statetocoords import StateToCoords
from scripts.mgltools.mglutil.math.transformation import Transformation
from scripts.mgltools.AutoDockTools.Conformation import Conformation, State
import random


class AutoDockMoleculePreparation:

    std_types = ['CYS', 'ILE', 'SER', 'VAL', 'GLN', 'LYS', 'ASN', 'PRO', 'THR', 'PHE', 'ALA', 'HIS', 'GLY', 'ASP', 'LEU', 'ARG', 'TRP', 'GLU', 'TYR', 'MET', 'HID', 'HSP', 'HIE', 'HIP', 'CYX', 'CSS']

    def __init__(self, molecule, mode='', repairs='',
                    charges_to_add=None, cleanup='',
                    outputfilename=None, debug=False,
                    version=3, delete_single_nonstd_residues=False):
        #debug = True
        self.debug = debug
        self.molecule = mol = molecule
        self.mode = mode
        self.lenNPHS = 0
        #fix atom elements if necessary:
        self.autodock_element_dict = {'Cl':'c', 'Fe': 'f', 'Br':'b'}
        self.reverse_autodock_elements = {'c':'Cl', 'f':'Fe','b':'Br'}
        self.version = version
        self.delete_single_nonstd_residues = delete_single_nonstd_residues
        if version==3:
            revlist = ['c','f','b']
            reversedAts = AtomSet(list(filter(lambda x: x.element in revlist, mol.allAtoms)))
            for at in reversedAts:
                if debug: print('restoring ', at.element, ' to ',)
                at.autodock_element = at.element
                at.element = self.reverse_autodock_elements[at.element]
                if debug: print(at.element, ' for processing...')
        # build bonds if necessary      #ALWAYS DO THIS!!!
        if not len(mol.allAtoms.bonds[0]): #what if not enough bonds
            mol.buildBondsByDistance()

        #PROCESS molecule
        self.cleanup_type_list = cleanup.split('_')
        #if mode=='automatic':  #DO this ANYWAY???
        #remove waters + chains composed of nonstd residues
        self.cleanUpResidues(mol, self.cleanup_type_list)

        # REPAIR: possibly add bonds to non-bonded atoms + add hydrogens
        self.repair_type_list = repairs.split('_')
        #if mode=='automatic':  #SHOULD THIS BE DONE ANYWAY???
        self.repairMol(mol, self.repair_type_list)
        ##3/25/2008: have to remove lps first because there is no atomic
        #number for Lp or LP in PyBabel/babelElements:
        merge_lonepairs = 'lps' in self.cleanup_type_list
        if merge_lonepairs:
            LPM = self.LPM = LonepairMerger()
            #lps = filter(lambda x: x.element=='Xx' and \
                  #(x.name[0]=='L' or x.name[1]=='L'), mol.allAtoms)
            lenLPS = self.lenLPS = LPM.mergeLPS(mol.allAtoms)
            if debug: print('merged ', lenLPS, ' lonepairs')
            self.cleanup_type_list.remove('lps')
        # CHARGES: ALWAYS add charges last after atoms are set but
        # before merging nphs and lps
        self.chargeType = charges_to_add
        self.chargeError = 'ERROR'
        if debug: print("charges_to_add=", charges_to_add)
        #if mode=='automatic':  #SHOULD THIS BE DONE ANYWAY???
        self.addCharges(mol, charges_to_add)
        if self.chargeType!='ERROR':
            self.chargeError = self.calcChargeError()
        #need a two step cleanup to do merges after adding bonds/hydrogens
        #and AFTER adding charges
        # CLEANUP: possibly merge nonpolarhydrogens and/or lonepairs
        self.cleanUpAtomTypes(mol, self.cleanup_type_list)
        #set autodock_element for all atoms
        mol.allAtoms.autodock_element = mol.allAtoms.element
        self.outputfilename = outputfilename
        if debug: print("end of base class init")


    def repairMol(self, mol, repair_list):
        #@@ what about metals?
        # possibly check for non-bonded atoms and add bonds to them
        # possibly add hydrogens 
        if 'bonds' in repair_list:
            non_bonded_atoms = mol.allAtoms.get(lambda x: len(x.bonds)==0)
            if non_bonded_atoms:
                if self.debug: print("!!adding bonds to ", non_bonded_atoms.name, "!!")
                for a in non_bonded_atoms:
                    nearestAt = self.findNearest(a, mol.allAtoms)
                    Bond(a, nearestAt, bondOrder=1, origin='AddedBond', check=0,
                                addBond=1)
                    if self.debug: print('added bond:', a.bonds[0])
            repair_list.remove('bonds')
        # possibly add hydrogens 
        #MODE SWITCH 1: adding hydrogens  ||IN USE||
        self.newHs = 0
        hs = mol.allAtoms.get(lambda x: x.element=='H')
        if 'hydrogens' in repair_list:
            if self.debug: print("addinghydrogens!")
            self.newHs = self.addHydrogens(mol)
            if self.debug: print("added ", self.newHs, " hydrogens")
            repair_list.remove('hydrogens')
        elif not hs and 'checkhydrogens' in self.repair_type_list:
            if self.debug: print("addinghydrogens from check hydrogens")
            self.newHs = self.addHydrogens(mol)
            if self.debug: print("added ", self.newHs, " hydrogens")
            repair_list.remove('checkhydrogens')

    def addHydrogens(self, mol):
        beforeLen = len(mol.allAtoms)
        #NB: this adds all hydrogens (?!?)
        HB = HydrogenBuilder()
        HB.addHydrogens(mol)
        afterLen = len(mol.allAtoms)
        newHs = afterLen - beforeLen
        #NB this adds the HYDROGENS as "ATOM" not HETATM
        return newHs


    def addCharges(self, mol, charges_to_add):
        #detect whether isPeptide first
        debug = self.debug
        if debug: print("in addCharges:mol.name=", mol.name, " and charges_to_add=", charges_to_add)
        isPeptide = mol.isPeptide = self.detectIsPeptide()
        if charges_to_add=='gasteiger':
            chargeCalculator = self.chargeCalculator = GasteigerChargeCalculator()
            if isPeptide:
                print("adding gasteiger charges to peptide")
                chargeCalculator = self.chargeCalculator = GasteigerChargeCalculator()
        elif charges_to_add=='Kollman':
            if not mol.isPeptide:
                print("adding Kollman charges to non-peptide")
            chargeCalculator = self.chargeCalculator = KollmanChargeCalculator()
        #KEEP charges from file, if there are any
        #MODE SWITCH 2: adding charges????  ||NOT IN USE||
        if debug: print( "self.chargeType=", self.chargeType)
        if self.chargeType is None:
            # Python3: filter() returns an iterator; use generator expressions for counts
            atoms = mol.chains.residues.atoms
            len_zero_charges = sum(1 for x in atoms if x.charge == 0)
            if len_zero_charges == len(mol.allAtoms):
                print("WARNING: all atoms in '%s' had zero charges! Adding gasteiger charges..." % mol.name)
                chargeCalculator = self.chargeCalculator = GasteigerChargeCalculator()
                chargeCalculator.addCharges(mol.allAtoms)
            else:
                len_charged = sum(1 for x in atoms if x.chargeSet is not None)
                if len_charged!=len(mol.chains.residues.atoms):
                    #charges could come from file or ?
                    print("WARNING: some atoms in '%s' had no charge! Adding gasteiger charges to all..."%mol.name)
                if debug: print("no change to charges!")
            self.chargeType = mol.allAtoms[0].chargeSet
        else: 
            chargeCalculator.addCharges(mol.allAtoms)
            if debug: print('added ' + charges_to_add + 'charges to ', mol.name)


    def cleanUpAtomTypes(self, mol, cleanup_list):
        debug = self.debug
        merge_nonpolar_hydrogens = 'nphs' in cleanup_list
        if merge_nonpolar_hydrogens:
            NPHM = self.NPHM = NonpolarHydrogenMerger()
            lenNPHS = self.lenNPHS = NPHM.mergeNPHS(mol.allAtoms)
            cleanup_list.remove('nphs')
        merge_lonepairs = 'lps' in cleanup_list
        if merge_lonepairs:
            LPM = self.LPM = LonepairMerger()
            #lps = filter(lambda x: x.element=='Xx' and \
                  #(x.name[0]=='L' or x.name[1]=='L'), mol.allAtoms)
            lenLPS = self.lenLPS = LPM.mergeLPS(mol.allAtoms)
            cleanup_list.remove('lps')
        if self.delete_single_nonstd_residues:
            #delete non-standard residues from within a single chain
            non_std = []
            for res in mol.chains.residues:
                if res.type not in self.std_types:
                    non_std.append(res)
            #delete non_std
            names = ""
            #remove all bonds to other atoms:
            for r in non_std:
                for a in r.atoms:
                    for b in a.bonds:
                        at2 = b.atom1
                        if at2==a: at2 = b.atom2
                        at2.bonds.remove(b)
                        a.bonds.remove(b)
                names = names + r.parent.id +r.name + '_'
            for res in non_std:
                res.parent.remove(res)
                del res
            #fix allAtoms short cut
            mol.allAtoms = mol.chains.residues.atoms
        delete_alternateB = 'deleteAltB' in cleanup_list        
        if delete_alternateB:
            alt_ats = mol.allAtoms.get(lambda x: '@' in x.name)
            if not len(alt_ats):
                return
            dict_alts = {}
            for a in alt_ats: 
                dict_alts[a.name[-1]] = 1
            alt_types = dict_alts.keys()

            alt_types.sort()
            Astr = "@" + alt_types[0]
            altA_ats = mol.allAtoms.get(lambda x: Astr in x.name)
            for t in alt_types[1:]:
                alt_s = '@' + t
                t_alt_ats = mol.allAtoms.get(lambda x: alt_s in x.name)
                ct = 0
                for a in t_alt_ats:
                    res = a.parent
                    for b in a.bonds:
                        at2 = b.atom1
                        if at2==a: at2 = b.atom2
                        at2.bonds.remove(b)
                        a.bonds.remove(b)
                    res.remove(a)
                    del a
                    ct = ct+1
            count = 0
            #rename any remaining @
            for a in altA_ats:
                ind = a.name.index('@')
                a.name = a.name[:ind]
                count = count+1
            mol.allAtoms = mol.chains.residues.atoms
            mol.allAtoms.number = list(range(1, len(mol.allAtoms)+1))


    def cleanUpResidues(self, mol, cleanup_list):
        self.lenHOHS = 0
        remove_waters = 'waters' in cleanup_list
        if remove_waters:
            hohs = mol.allAtoms.get(lambda x: x.parent.type=='HOH')
            #this line added at Alex Perryman's suggestion:
            hohs2 = mol.allAtoms.get(lambda x: x.parent.type=='WAT')
            hohs = hohs + hohs2 
            if hohs:
                #remove(hohs)
                lenHOHS = self.lenHOHS = len(hohs)
                for h in hohs:
                    for b in h.bonds:
                        c = b.atom1
                        if c==h:
                            c = b.atom2
                        c.bonds.remove(b)
                    h.bonds = BondSet()
                    res = h.parent
                    h.parent.remove(h)
                    if len(h.parent.children)==0:
                        res = h.parent
                        chain = res.parent
                        chain.remove(res)
                        if len(chain.children)==0:
                            mol = chain.parent
                            mol.remove(chain)
                            del chain
                        del res
                    del h 
                #fix allAtoms short cut
                mol.allAtoms = mol.chains.residues.atoms    
            cleanup_list.remove('waters')

        remove_nonstdres = 'nonstdres' in cleanup_list
        self.lenChainsRemoved = 0
        if remove_nonstdres: 
            if len(mol.chains)>1:
                chains_to_delete = []
                for c in mol.chains:
                    num_res = len(c.residues)
                    non_std = []
                    for res in c.residues:
                        if res.type not in self.std_types:
                            non_std.append(res)
                    if len(c.residues)==len(non_std):
                        chains_to_delete.append(c)
                if len(chains_to_delete):
                    self.lenChainsRemoved = len(chains_to_delete)
                    for c in chains_to_delete:
                        #remove(c)     #NB this would remove HOH and ligand from 1HSG
                        c.parent.remove(c)
                        del c
                    #fix allAtoms short cut
                    mol.allAtoms = mol.chains.residues.atoms    
            #if self.delete_single_nonstd_residues:
            #    #delete non-standard residues from within a single chain
            #    non_std = []
            #    for res in mol.chains.residues:
            #        if res.type not in self.std_types:
            #            non_std.append(res)
            #    #delete non_std
            #    names = ""
            #    #remove all bonds to other atoms:
            #    for r in non_std:
            #        for a in r.atoms:
            #            for b in a.bonds:
            #                at2 = b.atom1
            #                if at2==a: at2 = b.atom2
            #                at2.bonds.remove(b)
            #                a.bonds.remove(b)
            #        names = names + r.parent.id +r.name + '_'
            #    print "\'Deleting non-standard residues:" + names + " from " + mol.name
            #    for res in non_std:
            #        res.parent.remove(res)
            #        del res
            #    #fix allAtoms short cut
            #    mol.allAtoms = mol.chains.residues.atoms
            cleanup_list.remove('nonstdres')


    def calcChargeError(self):
        s = numpy.add.reduce(self.molecule.allAtoms.charge)
        return min(math.ceil(s)-s, s-math.floor(s))


    def detectIsPeptide(self):
        isPeptide = True
        d = {}
        for r in self.molecule.chains.residues:
            d[r.type] = 0
        for t in d.keys():
            if t not in q.keys():
                isPeptide = False
                break
        return isPeptide


    def findNearest(self, atom, bonded_atoms):
        lenK = len(bonded_atoms)
        lenC = 1
        nonbonded_atoms = AtomSet([atom])
        c = numpy.array(nonbonded_atoms.coords, 'f')
        k = numpy.array(bonded_atoms.coords, 'f')
        bigK = numpy.resize(k, (lenK, lenK, 3))
        c.shape = (1, 1, 3)
        bigM = bigK[:1]
        d = bigM - c
        dSQ = d*d
        dSQMAT = numpy.sum(dSQ, 2)
        mm = dSQMAT[0][0]
        mindex = 0
        for i in range(lenK):
            if dSQMAT[0][i]<mm:
                mm = dSQMAT[0][i]
                mindex = i
        #print "found closest atom %d at dist %f" %( mindex, mm)
        return bonded_atoms[mindex]



class ReceptorPreparation(AutoDockMoleculePreparation):
    """
    Facade for preparing a scripts.mgltools.MolKit molecule to be used as the receptor in an
    AutoDock experiment derived from AutoDockMoleculePreparation class.


    Receptor preparation involves:
        1. adding solvation parameters
        2. writing a pdbqs outputfile

    Receptor preparation Collaborators extend base class collaborators:
        in AutoDockTools:
            in atomTypeManager classes for conforming to AD atomtypes:
                SolvationParameterizer
       in this file:
            ReceptorWriter

"""

    def __init__(self, molecule, mode='automatic', repairs='checkhydrogens',
                    charges_to_add='Kollman',
                    cleanup='nphs_lps_waters_nonstdres',
                    outputfilename=None, debug=False,
                    preserved={},
                    delete_single_nonstd_residues=False, dict=None):

        AutoDockMoleculePreparation.__init__(self, molecule, mode, repairs,
                                            charges_to_add, cleanup,
                                            debug=debug,
                                            delete_single_nonstd_residues=False)

        self.dict = dict
        if len(preserved):
            for atom, chargeList in preserved.items():
                atom._charges[chargeList[0]] = chargeList[1]
                atom.chargeSet = chargeList[0]

        # add solvation parameters here:
        SP = self.SP = SolvationParameterizer()
        unknownatoms = SP.addParameters(molecule.chains.residues.atoms)


        self.writer = ReceptorWriter()
        #MODE SWITCH 5: write outputfile     ||IN USE||
        if mode=='automatic':
            #write outputfile now without waiting ...
            #fix this: what records should be written?
            self.write(outputfilename)
            self.outputfilename = outputfilename


    def summarize(self):
        mol = self.molecule
        msg = "setup "+ mol.name+ ":\n"
        if self.newHs!=0:
            msg = msg + " added %d new hydrogen(s)\n" %self.newHs
        if self.chargeType!=None:
            if self.chargeType in ['pdbq', 'mol2', 'pdbqt', 'pdbqs']:
                msg = msg + " kept charges from " + self.chargeType + " file\n"
            else:
                msg = msg + " added %s charges\n" %self.chargeType
        if self.lenNPHS:
            msg = msg + " merged %d non-polar hydrogens\n" %self.lenNPHS
        if self.lenLPS:
            msg = msg + " merged %d lone pairs\n" %self.lenLPS
        totalCh = numpy.add.reduce(mol.allAtoms.charge)
        if hasattr(self, 'chargeError') and abs(self.chargeError)>0.000005:
            msg = msg + " total charge error = %6.4f\n"%self.chargeError
        return msg


    def write(self, outputfilename=None):
        #write to outputfile
        mol = self.molecule
        if not outputfilename and not self.outputfilename:
            outputfilename = self.outputfilename = mol.name + ".pdbqs"
        #should this be done at Molecule level or at Atom level?
        #Atom level would allow multiple molecules in Receptor, at no price
        self.writer.write(mol, outputfilename)



class AD4ReceptorPreparation(AutoDockMoleculePreparation):
    """
    Facade for preparing a scripts.mgltools.MolKit molecule to be used as the receptor in an
    AutoDock4 experiment derived from AutoDockMoleculePreparation class.

    Receptor preparation involves:
        1. adding autodock_element types
        2. writing a pdbqt outputfile

    Receptor preparation Collaborators extend base class collaborators:
        in AutoDockTools:
            in atomTypeManager classes for conforming to AD atomtypes:
                AutoDock4_AtomTyper
       in this file:
            AD4ReceptorWriter

"""

    def __init__(self, molecule, mode='automatic', repairs='checkhydrogens',
                    charges_to_add='gasteiger',
                    cleanup='nphs_lps_waters_nonstdres',
                    outputfilename=None, debug=False,
                    version=4, preserved={},
                    delete_single_nonstd_residues=False,
                    dict=None):


        self.dict = dict
        AutoDockMoleculePreparation.__init__(self, molecule,
                        mode=mode, repairs=repairs,
                        charges_to_add=charges_to_add, cleanup=cleanup,
                        outputfilename=outputfilename, debug=debug,
                        version=version, delete_single_nonstd_residues=delete_single_nonstd_residues)

        #NEED TO TYPE ATOMS!!
        try:
            delattr(mol.allAtoms, 'autodock_element')
        except:
            pass
        atomTyper = AutoDock4_AtomTyper()
        atomTyper.setAutoDockElements(molecule, reassign=True) #4/15/05 catch here?
        if len(preserved):
            for atom, chargeList in preserved.items():
                atom._charges[chargeList[0]] = chargeList[1]
                atom.chargeSet = chargeList[0]

        self.writer = AD4ReceptorWriter()
        #MODE SWITCH 5: write outputfile     ||IN USE||
        if mode=='automatic':
            #write outputfile now without waiting ...
            #fix this: what records should be written?
            self.write(outputfilename)
            self.outputfilename = outputfilename


    def write(self, outputfilename=None):
        #write to outputfile
        mol = self.molecule
        if not outputfilename and not self.outputfilename:
            outputfilename = self.outputfilename = mol.name + ".pdbqt"
        #should this be done at Molecule level or at Atom level?
        #Atom level would allow multiple molecules in Receptor, at no price
        self.writer.write(mol, outputfilename)
        if self.dict is not None:
            if not os.path.exists(self.dict):
                fptr = open(self.dict, 'a')
                ostr = "summary = d = {}\n"
                fptr.write(ostr)
            else:
                fptr = open(self.dict, 'a')
            # setup dict
            ostr = "d['" + self.molecule.name +"'] = {"
            # autodock_element
            ostr = ostr + "'atom_types': ["
            atom_types = list(set(self.molecule.chains.residues.atoms.autodock_element))
            for t in atom_types[:-1]:
                ostr = ostr + "'%s', "%t
            ostr = ostr + "'%s' "%atom_types[-1]
            # number of atoms
            ostr = ostr + "],\n\t\t\t'number of atoms' : "
            num_atoms = len(mol.chains.residues.atoms)
            ostr = ostr + "%d"%num_atoms
            ostr = ostr + ",\n\t\t\t'res_types' : ["
            # residue types
            res_types = list(set(self.molecule.chains.residues.type))
            res_types.sort()
            res_ct = 0
            for t in res_types[:-1]:
                ostr = ostr + "'%s', "%t
                res_ct += 1
                if res_ct==7:
                    ostr  = ostr + '\n\t\t\t\t\t\t   '
                    res_ct = 0
            ostr = ostr + "'%s' "%res_types[-1]
            ostr = ostr + "],\n\t\t\t'total_charge' : "
            # net charge
            total_charge = numpy.add.reduce(self.molecule.allAtoms.charge)
            ostr = ostr + "%6.2f"%total_charge
            ostr = ostr + ",\n\t\t\t'zero_charge' : ["
            # atoms with zero charge
            zc_ats = self.molecule.allAtoms.get(lambda x: x.charge==0)
            for at in zc_ats[:-1]:
                ostr = ostr + "'%s', "%at.name
            if len(zc_ats):
                ostr = ostr + "'%s'"%zc_ats[-1].name
            ostr = ostr + "],\n\t\t\t'zero_bonds' : ["
            # atoms with zero bonds
            zbond_ats = self.molecule.allAtoms.get(lambda x: len(x.bonds)==0)
            for at in zbond_ats:
                ostr = ostr + "'%s', "%at.name
            ostr = ostr + "],\n\t\t\t}\n"
            fptr.write(ostr)
            fptr.close()
        


class ReceptorWriter:
    """
    ReceptorWriter class writes a receptor molecule which has a ReceptorPreparationObject to a file
"""

    def __init__(self, write_CONECT=False):
        self.writer = PdbqsWriter()
        #each atom must have charge and AtVol and AtSolPar
        self.conditions = [lambda x: hasattr(x, 'AtVol'), \
                           lambda x: hasattr(x, 'AtSolPar'), \
                           lambda x: x.chargeSet is not None]
        self.write_CONECT = write_CONECT


    def write(self, receptor, outputfile):
        #should this be done at Molecule level or at Atom level?
        receptor_atoms = receptor.chains.residues.atoms
        len_receptor_atoms = len(receptor_atoms)

        ##check each condition for receptor_atoms
        for cond in self.conditions:
            assert sum(1 for a in receptor_atoms if cond(a))==len_receptor_atoms

        outptr = open(outputfile, 'w')
        records = ['ATOM']
        #if self.write_CONECT:
        #    records.append('CONECT')
        #for at in receptor_atoms:
        #    self.writer.write_atom(outptr, at)
        for ch in receptor.chains:
            for at in ch.residues.atoms:
                self.writer.write_atom(outptr, at)
            rec = self.writer.defineTERRecord(at)
            outptr.write(rec)
        #put in the TER and END lines???
        #optional CONECT records???
        if self.write_CONECT:
            self.writeCONECTRecords(receptor, outptr)
        outptr.close()
        

    def writeCONECTRecords(self, fptr, mol):
        atms = mol.allAtoms
        for atm in atms:
            rec = 'CONECT%5i'%atm.number
            for b in atm.bonds:
                a2 = b.atom1
                if a2==atm: a2 = b.atom2
                if not a2 in atms: continue #don't write intermolecular bonds
                rec = rec + '%5i'%a2.number 
            rec = rec + '\n'
            fptr.write(rec)




class AD4ReceptorWriter(ReceptorWriter):
    """
    AD4ReceptorWriter class writes a receptor molecule which has a ReceptorPreparationObject to a pdbqt file
"""

    def __init__(self, write_CONECT=False):
        self.writer = PdbqtWriter()
        self.conditions = [lambda x: hasattr(x, 'autodock_element'), \
                           lambda x: x.chargeSet is not None]
        self.write_CONECT = write_CONECT


class LigandPreparation(AutoDockMoleculePreparation):
    def __init__(self, molecule, mode='automatic', repairs='hydrogens_bonds',
                    charges_to_add='gasteiger', cleanup='nphs_lps',
                    allowed_bonds='backbone',  #amide + guanidinium off @@ backbone on (5/2014)
                    root='auto', outputfilename=None, dict=None, debug=False,
                    check_for_fragments=False, bonds_to_inactivate=[],
                    inactivate_all_torsions=False,
                    version=3, limit_torsions=False,
                    delete_single_nonstd_residues=False,
                    detect_bonds_between_cycles=False):
        # why aren't repairs, cleanup and allowed_bonds lists??
        # to run tests: use allowed_bonds = 'guanidinium'

        self.detect_bonds_between_cycles = detect_bonds_between_cycles


        AutoDockMoleculePreparation.__init__(self, molecule, mode, repairs,
                                            charges_to_add, cleanup,
                                            debug=debug, version=version,
                                            delete_single_nonstd_residues=False)

        self.version = version
        self.charge_error_tolerance = 0.000005
        #FOR FLEXIBILITY MODEL: setup RotatableBondManager
        #process bonds
        #MODE SWITCH 4: disallow allowed_bonds ||NOT IN USE||
        self.RBM = RotatableBondManager(molecule,
                            allowed_bonds, root, debug=False,
                            check_for_fragments=check_for_fragments,
                            bonds_to_inactivate=bonds_to_inactivate,
                            detectAll=self.detect_bonds_between_cycles)
        if inactivate_all_torsions is True:
            #in this case, make all active torsions inactive
            self.RBM.set_all_torsions(False)
        elif limit_torsions is not False:
            self.RBM.limit_torsions(limit_torsions)

        #detect aromatic cycle carbons and rename them for AutoDock3 
        rename = self.version==3
        ACM = self.ACM = AromaticCarbonManager(rename=rename)
        self.aromCs = self.ACM.setAromaticCarbons(molecule) 

        #optional output summary filename to write types and number of torsions
        self.dict = dict

        self.outputfilename = outputfilename
        #if mode is 'automatic': write outputfile now 
        #without waiting to set torsions etc interactively
        #MODE SWITCH 5: write outputfile     ||IN USE||
        self.writer = LigandWriter()
        if mode=='automatic':
            self.write(outputfilename)
            self.outputfilename = outputfilename


    def summarize(self):
        mol = self.molecule
        msg = "setup "+ mol.name+ ":\n"
        if self.newHs!=0:   
            msg = msg + " added %d new hydrogen(s)\n" %self.newHs
        if self.chargeType!=None:   
            if self.chargeType in ['pdbq', 'mol2']:
                msg = msg + " kept charges from " + self.chargeType + " file\n"
            else:
                msg = msg + " added %s charges\n" %self.chargeType
        if self.lenNPHS:
            msg = msg + " merged %d non-polar hydrogens\n" %self.lenNPHS
        if self.lenLPS:
            msg = msg + " merged %d lone pairs\n" %self.lenLPS
        if len(self.aromCs):
            msg = msg + " found %d aromatic carbons\n" %len(self.aromCs)
        totalCh = numpy.add.reduce(mol.allAtoms.charge)
        msg = msg + " detected %d rotatable bonds\n" %len(mol.possible_tors_bnds)
        msg = msg + " set TORSDOF to %d\n"%mol.TORSDOF
        #if hasattr(self, 'chargeError') and abs(self.chargeError)>0.000005:
        if hasattr(self, 'chargeError') and abs(self.chargeError)>self.charge_error_tolerance:
            msg = msg + " total charge error = %6.4f\n"%self.chargeError
        return msg

    def set_autodock_element_types_for_writing(self, allAtoms):
        sublist = self.autodock_element_dict.keys()
        problem_ats = AtomSet(list(filter(lambda x: x.element in sublist, allAtoms)))
        for at in problem_ats:
            #there is a substitute 'name' for this atom

            newval = self.autodock_element_dict[at.element]
            at.name = newval
            #have to change element to get correct written output
            #because the stupid pdbWriter writes the element not the name
            at.element = newval
            at.autodock_element = at.element

        return problem_ats


    def restore_autodock_element_types_after_writing(self, problem_ats):
        for at in problem_ats:
            #there is a substitute 'name' for this atom

            newval = self.reverse_autodock_elements[at.element]
            #have to change element to get correct written output
            #because the stupid pdbWriter writes the element not the name
            at.element = newval



    def write(self, outputfilename=None, write_CONECT=False):
        #write ligand to outputfile
        mol = self.molecule
        if not hasattr(mol, 'ROOT'):

            return 'ERROR' 
        if hasattr(mol, 'torTree') and not hasattr(mol, 'torscount'):
            mol.torscount = len(mol.torTree.torsionMap)
        if not outputfilename:
            stem = os.path.splitext(os.path.basename(mol.parser.filename))[0]
            if self.version>=4:
                outputfilename = stem + ".pdbqt"
            else:
                outputfilename = stem + ".pdbq"

        problem_ats = self.set_autodock_element_types_for_writing(mol.allAtoms)

        self.writer.write(mol, outputfilename)
        self.outputfilename = outputfilename
        if self.dict is not None:
            if not os.path.exists(self.dict):
                fptr = open(self.dict, 'a')
                ostr = "summary = d = {}\n"
                fptr.write(ostr)
            else:
                fptr = open(self.dict, 'a')
            type_dict = {}
            for a in self.molecule.allAtoms:
                type_dict[a.autodock_element] = 1
            atom_types = type_dict.keys()
            atom_types.sort()
            ostr = "d['" + self.molecule.name +"'] = {"
            ostr = ostr + "'atom_types': [" 
            for t in atom_types[:-1]:
                ostr = ostr + "'%s', "%t
            ostr = ostr + "'%s' "%atom_types[-1]
            ostr = ostr + "],\n\t\t\t'rbonds':" + str(self.molecule.torscount)
            ostr = ostr + ",\n\t\t\t'zero_charge' : ["
            zc_ats = self.molecule.allAtoms.get(lambda x: x.charge==0)
            for at in zc_ats:
                ostr = ostr + "'%s', "%at.name
            ostr = ostr + "],\n\t\t\t}\n"
            fptr.write(ostr)
            fptr.close()
        self.restore_autodock_element_types_after_writing(problem_ats)



    def autoroot(self):
        self.RBM.autoroot()


    def setroot(self, index):
        self.RBM.setroot(index)


    def limit_torsions(self, val, type):
        self.RBM.limit_torsions(val, type)


    def set_all_torsions(self, flag):
        mol = self.molecule
        self.RBM.set_all_torsions(flag)


    def set_amide_torsions(self, flag):
        if not len(self.molecule.amidebnds):
            return
        self.RBM.set_amide_torsions(flag)
        
        
    def set_ppbb_torsions(self, flag):
        if not len(self.molecule.ppbbbnds):
            return
        self.RBM.set_peptidebackbone_torsions(flag)
        
        
    def set_guanidinium_torsions(self, flag):
        if not len(self.molecule.guanidiniumbnds):
            return
        self.RBM.set_guanidinium_torsions(flag)
        

    def toggle_torsion(self,ind1, ind2):
        self.RBM.toggle_torsion(ind1, ind2)
        

    def changePlanarityCriteria(self, cutoff):
        oldAromCs = self.aromCs
        self.aromCs = self.ACM.setAromaticCarbons(self.molecule, cutoff)


    def set_carbon_names(self, atoms, type):
        self.ACM.set_carbon_names(atoms, type)
        #have to reset this field here
        self.aromCs = AtomSet(list(filter(lambda x: (x.autodock_element=='A' and x.element=='C'),
                                          self.molecule.allAtoms)))


    def setAromaticCarbons(self):
        self.aromCs = self.ACM.setAromaticCarbons(self.molecule)
        return self.aromCs

#LigandPreparation

class AD4LigandPreparation(LigandPreparation):

    def __init__(self, molecule, mode='automatic', repairs='checkhydrogens',
                    charges_to_add='gasteiger', cleanup='nphs_lps_waters_nonstdres',
                    allowed_bonds='backbone',  #@@ amide, guanidinium off, backbone on (5/2014)
                    root='auto', outputfilename=None, dict=None, debug=False,
                    check_for_fragments=False, bonds_to_inactivate=[],
                    inactivate_all_torsions=False,
                    version=4, typeAtoms=True, limit_torsion_number=False, 
                    limit_torsions_type='fewest', limit_torsions=False,
                    delete_single_nonstd_residues=False,
                    attach_nonbonded_fragments=False,
                    detect_bonds_between_cycles=False, 
                    attach_singletons=False,
                    write_CONECT=False):


        if attach_nonbonded_fragments==True:
            molecule.attach_nonbonded_fragments(attach_singletons=attach_singletons)
        #FIX THIS: what if molecule already has autodock_element set???
        LigandPreparation.__init__(self, molecule, mode='interactive', 
                    repairs=repairs, charges_to_add=charges_to_add, 
                    cleanup=cleanup, allowed_bonds=allowed_bonds, 
                    root=root, outputfilename=outputfilename, 
                    dict=dict, debug=debug,
                    check_for_fragments=check_for_fragments,
                    bonds_to_inactivate=bonds_to_inactivate, 
                    inactivate_all_torsions=inactivate_all_torsions,
                    version=version,
                    limit_torsions=limit_torsions,
                    detect_bonds_between_cycles=detect_bonds_between_cycles)

        # AD4 force field uses total rotatable bonds for TORSDOF:
        # calculation of loss of entropy on ligand binding
        #molecule.TORSDOF = molecule.possible_tors
        self.charge_error_tolerance = 0.008
        if typeAtoms:
            #NEED TO TYPE ATOMS!!
            delattr(molecule.allAtoms, 'autodock_element')
            atomTyper = AutoDock4_AtomTyper()
            atomTyper.setAutoDockElements(molecule, reassign=True) #4/15/05:catch here?
        #if self.debug: print "allAtoms.autodock_element=", molecule.allAtoms.autodock_element
        molecule.TORSDOF = molecule.possible_tors - len(molecule.guanidiniumbnds+molecule.amidebnds)
        #could limit the number of torsions here
        if limit_torsion_number is not False:
            ndihe = limit_torsion_number
            type = limit_torsions_type
            self.RBM.limit_torsions(ndihe, type)

        self.writer = AD4LigandWriter(write_CONECT=write_CONECT)
        #MODE SWITCH 5: write outputfile     ||IN USE||
        if mode=='automatic':
            #write outputfile now without waiting ...
            #fix this: what records should be written?
            self.write(outputfilename)
            self.outputfilename = outputfilename


    def set_autodock_element_types_for_writing(self, allAtoms):
        pass 

    def restore_autodock_element_types_after_writing(self, problem_atoms):
        pass



class LigandWriter:

    """
    LigandWriter class writes a ligand molecule which has a LigandPreparationObject to a file
    """

    def __init__(self, debug=False, write_CONECT=False):
        self.writer = PdbqWriter()
        self.maxtors = 32
        self.debug = debug
        self.write_CONECT = write_CONECT


    def write(self, ligand, outputfile):
        #setup rnum field
        ligand.allAtoms.rnum = -1

        old_torscount=ligand.torscount
        ligand.torscount = len(ligand.allAtoms.bonds[0].get(lambda x:x.activeTors))


        ttc = ligand.torscount

        assert hasattr(ligand, 'ROOT') # is this necessary or redundant
        assert hasattr(ligand, 'TORSDOF')
        self.outfptr = open(outputfile, 'w')
        #FIRST write out remarks about torsions
        if (ttc>self.maxtors):
            self.outfptr.write("REMARK WARNING: %d MAX_TORS EXCEEDED!!!\n"%(self.maxtors))
        self.outfptr.write("REMARK  " +"%d" %(ttc) + " active torsions:\n")
        self.outfptr.write("REMARK  status: ('A' for Active; 'I' for Inactive)\n")
        self.outatom_counter = 1
        bondset = ligand.allAtoms.bonds[0]
        for i, at in enumerate(ligand.allAtoms):
            at.number = i+1
        rootlist = [ligand.ROOT]
        rootnum = 1
        for b in bondset:
            if b.activeTors==1:
                bTors = 'A'
            else:
                bTors = 'I'
            # at1str = b.atom1.name + '_' + str(b.atom1.number)
            # at2str = b.atom2.name + '_' + str(b.atom2.number)
            # Use stable integer indices (do not rely on AtomSet-wide `number` sequences)
            at1_idx = ligand.allAtoms.index(b.atom1) + 1
            at2_idx = ligand.allAtoms.index(b.atom2) + 1
            at1str = "%s_%d" % (b.atom1.name, at1_idx)
            at2str = "%s_%d" % (b.atom2.name, at2_idx)
            if b.possibleTors :
                if bTors=='A' :
                    outstring = "REMARK " +" %3d  %s    between atoms: %-3s  and  %-3s \n" %(self.outatom_counter,bTors, at1str, at2str)
                    b.outNumber = self.outatom_counter
                    self.outatom_counter = self.outatom_counter + 1
                else:
                    outstring = "REMARK " +"      %s    between atoms: %-2s  and  %-3s \n" %(bTors, at1str, at2str)
                    #outstring = "REMARK " +"      %s    between atoms: %-3s  and  %-3s \n" %(bTors,b.atom1.name,b.atom2.name)
                self.outfptr.write(outstring)
        #next write out  root
        self.outfptr.write("ROOT\n")
        self.outatom_counter = 1
        #reset used field to serve as visited flag
        ligand.allAtoms.used = 0
        #rootlist grows to include atoms up to first active tors in each subtree
        #this starts w/ 1 atom:
        self.writtenAtoms = []
        for at in rootlist:
            #to start, rootlist has just 1 atom
            at.used = 1
            for bond in at.bonds:
                #if bond.activeTors and bond.possibleTors:continue
                at2 = bond.atom1
                if at2==at:
                    at2 = bond.atom2
                if bond.activeTors and bond.possibleTors:
                    #if at2 not in branchAtoms:
                        #branchAtoms.append(at2)
                    continue
                if at2.used:
                    continue
                if at2 not in rootlist:
                    rootlist.append(at2)
                    at2.rnum0 = rootnum
                    rootnum = rootnum + 1
        #before writing, sort rootlist
        newrootlist = []
        oldrootlist = rootlist
        r_ctr = 0
        for at in ligand.chains.residues.atoms:
            if hasattr(at, 'rnum0'):
                newrootlist.append(at)
                at.rnum = r_ctr
                r_ctr = r_ctr+1
        for at in newrootlist:
            at.number = self.outatom_counter
            self.writer.write_atom(self.outfptr,at)
            at.newindex = self.outatom_counter
            at.used = 1
            self.writtenAtoms.append(at)
            self.outatom_counter = self.outatom_counter+1
        #at this point, add chain_roots if there are any
        #for ch in dict['chain_rootlist']:
        #    for at in ch.residues.atoms:
        #        at.number = self.outatom_counter
        #        at.rnum = r_ctr
        #        r_ctr = r_ctr+1
        #        self.writer.write_atom(self.outfptr,at)
        #        at.newindex = self.outatom_counter
        #        at.used = 1
        #        self.writtenAtoms.append(at)
        #        self.outatom_counter = self.outatom_counter+1
        self.outfptr.write("ENDROOT\n")
        #last write out the rest of the stuff, using WriteSubtree.....
        #for at in rootlist:
        for at in rootlist:
            for bond in at.bonds:
                at2 = bond.atom1
                if at2==at:
                    at2 = bond.atom2
                if at2.used:
                    continue
                self.process(at, at2)
###                marker = self.outatom_counter
###                outstring = "BRANCH %3d %3d\n"%(at.rnum+1, marker)
###                self.outfptr.write(outstring)
###                ###11/30: change #1-> write atom before calling WriteSubtree
###                at2.newindex = self.outatom_counter
###                at2.number = self.outatom_counter
###                self.writer.write_atom(self.outfptr,at2)
###                self.writtenAtoms.append(at2)
###                self.outatom_counter = self.outatom_counter+1
###                self.WriteSubtree(at,at2)
###                outstring = "ENDBRANCH %3d %3d\n"%(at.rnum +1,marker)
###                self.outfptr.write(outstring)
        outstring = "TORSDOF " + str(ligand.TORSDOF) + "\n"
        self.outfptr.write(outstring)
        #call write_CONECT here
        if self.write_CONECT:
            self.writeCONECTRecords(ligand, self.outfptr)
        self.outfptr.close()
        ligand.returnCode = 0

        if len(self.writtenAtoms)!=len(ligand.allAtoms):
            ligand.returnCode = 1
            allAtoms = len(ligand.allAtoms)
            notWritten = allAtoms - len(self.writtenAtoms)
            ligand.returnMsg = 'WARNING: %d atoms of %d in %s  were not written\n' %(notWritten, allAtoms, ligand.parser.filename)
        for a in ligand.allAtoms:
            if hasattr(a, 'rnum0'):
                delattr(a, 'rnum0')
            if hasattr(a, 'cycleout'):
                delattr(a, 'cycleout')
            if hasattr(a, 'newindex'):
                a._ni = a.newindex
                delattr(a, 'newindex')
        ligand.ROOT.rnum0 = 0
        ligand.outputfile = outputfile


    def process(self, fromAtom, nextAtom):
        startIndex = fromAtom.number
        endIndex = self.outatom_counter
        new_branch = 0
        if hasattr(nextAtom,'used') and nextAtom.used: return
        for bond in nextAtom.bonds:
            if bond.neighborAtom(nextAtom)!=fromAtom:
                nAt = bond.neighborAtom(nextAtom)
                if hasattr(nAt, 'used') and nAt.used: continue
                else: new_branch = 1
        if new_branch:
            outstring = "BRANCH %3d %3d\n"%(startIndex, endIndex)
            #outstring = "BRANCH %3d %3d\n"%(fromAtom.rnum+1, marker)
            self.outfptr.write(outstring)
        queue = self.writeBreadthFirst(self.outfptr, fromAtom, nextAtom)
        if len(queue):
            for fromAtom, nextAtom in queue:
                self.process(fromAtom, nextAtom)
        if new_branch:
            outstring = "ENDBRANCH %3d %3d\n"%(startIndex, endIndex)
            self.outfptr.write(outstring)


    def writeLevel(self, atom, outfptr):
        """
        write all atoms bonded to atoms bonded to this atom by non-rotatable
        bonds
        """

        queue = []
        nextAts = []
        for b in atom.bonds:
            if b.activeTors:
                at2 = b.atom1
                if at2==atom: at2=b.atom2
                queue.append((atom, at2))
                continue
            a2 = b.atom1
            if a2==atom:
                a2 = b.atom2
            if a2.used:
                continue
            if a2 not in self.writtenAtoms:
                a2.number = a2.newindex = self.outatom_counter
                self.writer.write_atom(outfptr, a2)
                self.writtenAtoms.append(a2)
                a2.used = 1
                self.outatom_counter+=1
                nextAts.append(a2)
        for a2 in nextAts:
            nq = self.writeLevel(a2, outfptr)
            if len(nq):
                queue.extend(nq)
        return queue


    def writeBreadthFirst(self, outfptr, fromAtom, startAtom):
        """
            None <-writeBreadthFirst(outfptr, fromAtom, startAtom)
            writeBreadthFirst visits all the atoms in the current level
            then the first level down etc in a Breadth First Order traversal.
                            1                <-1
                        5 6   7 8            <-3
                     9 10   11 12            <-4
            It is used to write out the molecule with the correct format
            for AutoDock. Specifically, appropriate BRANCH/ENDBRANCH
            statements are added.
        """

        queue = []
        if startAtom.used==0:
            startAtom.used = 1
            startAtom.number = startAtom.newindex = self.outatom_counter
            self.writer.write_atom(outfptr,startAtom)
            self.writtenAtoms.append(startAtom)
            self.outatom_counter += 1
            activeTors = []
            #outfptr.write(outstring)
            for bond in startAtom.bonds:
                at2 = bond.atom1
                if at2==startAtom: at2 = bond.atom2
                if at2==fromAtom: continue  #this is current bond
                elif not at2.used:
                    if bond.activeTors:
                        queue.append((startAtom,at2))
                    else:
                        at2.number = at2.newindex = self.outatom_counter
                        self.writer.write_atom(outfptr, at2)
                        at2.written = 1
                        self.writtenAtoms.append(at2)
                        at2.newindex = self.outatom_counter
                        self.outatom_counter = self.outatom_counter + 1
                        newQ = self.writeLevel(at2, outfptr)
                        at2.used = 1
                        if len(newQ):
                            queue.extend(newQ)
        return  queue


    def WriteSubtree(self,fromAtom, startAtom):
        """WriteSubtree recursively visits the atoms in the current 'subtree' of the molecule
        in a BREADTH First Order traversal. It is used to write out the molecule with
        the correct format for AutoDock. Specifically, appropriate BRANCH/ENDBRANCH and
        TORS/ENDTORS statements are added.
        """

        if startAtom.used==0:
            startAtom.used = 1
            at = startAtom
            for bond in startAtom.bonds:
                if bond.activeTors:
                    continue
                marker = self.outatom_counter
                nextAtom = bond.atom1
                if nextAtom==startAtom:
                    nextAtom = bond.atom2
                if nextAtom==fromAtom:
                    continue
                if not nextAtom.used:
                    if hasattr(bond,'incycle'):
                        if not hasattr(nextAtom, 'cycleout'):
                            nextAtom.cycleout = 1
                            nextAtom.newindex = self.outatom_counter
                            nextAtom.number = self.outatom_counter
                            self.writer.write_atom(self.outfptr,nextAtom)
                            self.writtenAtoms.append(nextAtom)
                            self.outatom_counter = self.outatom_counter+1
                    else:
                        nextAtom.newindex = self.outatom_counter
                        nextAtom.number = self.outatom_counter
                        self.writer.write_atom(self.outfptr,nextAtom)
                        self.writtenAtoms.append(nextAtom)
                        self.outatom_counter = self.outatom_counter+1
            for bond in startAtom.bonds:
                marker = self.outatom_counter
                nextAtom = bond.atom1
                if nextAtom==startAtom:
                    nextAtom = bond.atom2
                if nextAtom==fromAtom:
                    #if self.debug:print "nextAtom%s==fromAtom%s"%(nextAtom.full_name(),fromAtom.full_name())
                    continue
                if not nextAtom.used: #skip this bond because bonded_atom has already been written
                    isleaf = len(nextAtom.bonds)==1
                    testcond = len(nextAtom.bonds)>1
                    #testcond: isleaf = len(nextAtom.bonds[0].neighborAtom(nextAtom).bonds)==1
                    if bond.activeTors and bond.possibleTors and not isleaf:
                        if testcond >0:
                            outstring = "BRANCH %3d %3d\n"%(at.newindex,marker)
                            self.outfptr.write(outstring)
                        nextAtom.newindex = self.outatom_counter
                        nextAtom.number = self.outatom_counter
                        self.writer.write_atom(self.outfptr,nextAtom)
                        self.writtenAtoms.append(nextAtom)
                        self.outatom_counter = self.outatom_counter+1
                    self.WriteSubtree(startAtom, nextAtom)
                    if bond.activeTors and bond.possibleTors:
                        if testcond >0:
                            outstring = "ENDBRANCH %3d %3d\n"%(at.newindex,marker)
                            self.outfptr.write(outstring)
        return


    def writeCONECTRecords(self, mol, fptr):
        #print "in writeCONECTRecords"
        atms = mol.allAtoms
        for atm in atms:
            rec = 'CONECT%5i'%atm.number
            for b in atm.bonds:
                a2 = b.atom1
                if a2==atm: a2 = b.atom2
                if not a2 in atms: continue #don't write intermolecular bonds
                rec = rec + '%5i'%a2.number
            rec = rec + '\n'
            fptr.write(rec)



class AD4LigandWriter(LigandWriter):
    """
    AD4LigandWriter class writes a ligand molecule which has a LigandPreparationObject to a file
    """

    def __init__(self, debug=False, write_CONECT=False):
        LigandWriter.__init__(self, write_CONECT)
        self.writer = PdbqtWriter()
        self.maxtors = 32
        self.debug = debug



class LigandRandomizer:
    """
    LigandRandomizer class tries to find a new random conformation of a ligand molecule and write it to a file
    Currently, ligand must have been created from previously written pdbqt file:
        i.e. it must have a torTree which is built from parsing a previously written pdbqt file
    """
    def QuatToAxisAngle(self, q):
        #(x,y,z,w)
        assert len(q)==4
        abs_w = abs(q[3])
        if q[3] >1.: q[3]=1.
        if q[3] <-1.: q[3]=-1.
        if q[3]==1:
            return [1.,0.,0.,0.]
        #conversion:
        angle = 2.* math.acos(q[3])
        inv_sin_half_angle = 1./math.sin(angle/2.)
        quat = []
        quat.append(q[0]*inv_sin_half_angle)
        quat.append(q[1]*inv_sin_half_angle)
        quat.append(q[2]*inv_sin_half_angle)
        w = q[3]
        #retval.ang = WrpModRad( angle );#line233 qmultiply.cc
        # define WrpModRad(angle)  WrpRad(ModRad(angle))
        # define WrpRad(angle)     IF (((angle)> PI)? ((angle)-TWOPI)
        #                          else IF (((angle)< -PI)?  ((angle)+TWOPI))
        #                          else: angle
        #               elif angle<-PI: ModRad(angle)+TWOPI
        #               if angle>PI: ModRad(angle)-TWOPI
        # angle == q[3]
        #WrpModRad(a) WrpRad( ModRad((a)))
        #^
        #WrpRad(a) (((a)>PI)? (ModRad(a)-TWOPI):(((a)<-1pi) ?  (ModRad(a)+TWOPI):(a)))
        #ModRad(a)  fmod(a, TWOPI)
        #ang = WrpModRad(angle)
        a = q[3] % 2*math.pi #ModRad(angle)
        if a > math.pi: #(a)>PI? ModRad(a) - TWOPI
            val = a - 2* math.pi
        elif a < -1*math.pi:
            val = a + 2*math.pi
        else:
            val = a
        quat.append(val)
        return normAxisAngle(quat)


    def normAxisAngle(self, quat):
        mag3 = hypotenuse(quat[0], quat[1], quat[2])
        ret = [1.,0.,0.,0]
        if mag3>0.000001:
            inv_mag3 = 1./mag3
            ret = [quat[0]*inv_mag3, quat[1]*inv_mag3, quat[2]*inv_mag3, quat[3]]
        return ret


    def __init__(self, molecule, outputfilename=None, #TRANSLATION
                        info_lo = "-1,-1,-1",#info->lo[X,Y,Z]
                        info_hi = "1,1,1",   #info->hi[X,Y,Z]
                        #qrange = 2.,   # see mkRandomState for quaternion
                        #drange = "-180,180",  # NO support for torsion constraints
                        ntries = 100,
                        torsOnly = 0,
                        relTrans = 0,
                        coord_slot = 1,
                        verbose=True):
        #validate input
        self.ok = False
        self.molecule = molecule
        m = molecule
        if verbose:
            msg =  "LR: molecule=%s\ninfo_lo=%s\ninfo_hi=%s\ntorsOnly=%d\nrelTrans=%d\n"%(molecule.name, info_lo,info_hi,torsOnly,relTrans)
            print(msg)
        #CHECK FOR TORTREE HERE:
        assert hasattr(m, 'torTree'), "specified ligand does not have a torsion tree"
        if m.parser.keys.count("ROOT")>0:
            ndihe = m.parser.keys.count('BRANCH')
        else:
            #do something else here @@ ??FILTER OUT ALL INACTIVE @@#1
            ndihe = len(m.allAtoms.bonds[0].get(lambda x: x.activeTors==True))
        if outputfilename is None:
            outputfilename= m.name + "_rand.pdbqt"

        m.buildBondsByDistance()
        # first time only
        if not hasattr(m, 'orig_coords'):
            m.orig_coords = m.allAtoms.coords[:] #need to reset between tries
            if not hasattr(m, 'mol_cen'):
                m.mol_cen = m.getCenter()
            m.orig_coord_index = coord_slot-1 #preserve original coords in slot '0'
        len_at_coords = len(m.allAtoms[0]._coords)
        if len_at_coords==1: #do this only once
            m.allAtoms.addConformation(m.orig_coords)
            m.allAtoms.setConformation(coord_slot) #set to 1
        else:
            #reset both copies to original coords
            for ix in range(len(m.allAtoms)):
                m.allAtoms[ix]._coords[m.orig_coord_index]=m.orig_coords[ix]
                m.allAtoms[ix]._coords[coord_slot]=m.orig_coords[ix]
            m.allAtoms.setConformation(coord_slot) #set to 1
        # SETUP StateToCoords
        if not hasattr(m, 'stoc'):
            m.stoc = StateToCoords(m, m.mol_cen, coord_slot)
        m.allAtoms.setConformation(coord_slot)
        #build reference dictionary of original bonds
        orig = len(m.allAtoms.bonds[0])
        orig_d = {}
        for a in m.allAtoms: orig_d[a]=set(a.bonds.getAtoms())
        # setup random torsion angles:
        dihe = []
        for ind in range(ndihe):
            dihe.append(round(random.uniform(-180,180),3))
        if torsOnly:
            trans = [0,0,0]
            #APRIL4
            m.stoc.applyAngList(dihe,Transformation(trans=m.mol_cen).getMatrix(transpose=1))
            self.ok = True
            newCrds = m.allAtoms.coords
            m.parser.write_with_new_coords(newCrds, outputfilename)
        else:
            #try a new random state up to ntries times = 100
            # convert trans to space centered on m
            #info_loX, info_loY, info_loZ = map(float, info_lo.split(','))
            #info_hiX, info_hiY, info_hiZ = map(float, info_hi.split(','))
            info_lo_list = info_lo.split(',')
            info_hi_list = info_lo.split(',')
            info_loX, info_loY, info_loZ = map(float, info_lo_list)
            info_hiX, info_hiY, info_hiZ = map(float, info_hi_list)
            self.info_loX = info_loX - m.mol_cen[0]
            self.info_loY = info_loY - m.mol_cen[1]
            self.info_loZ = info_loZ - m.mol_cen[2]
            self.info_hiX = info_hiX - m.mol_cen[0]
            self.info_hiY = info_hiY - m.mol_cen[1]
            self.info_hiZ = info_hiZ - m.mol_cen[2]
            # choose translation from these intervals:
            xrange = info_hiX-info_loX
            yrange = info_hiY-info_loY
            zrange = info_hiZ-info_loZ
            m.state = {}
            for new_try in range(ntries):
                #reset coords in working slot to original coords for new try
                m.allAtoms.setConformation(coord_slot)
                m.allAtoms.updateCoords(m.orig_coords, ind=coord_slot)
                #from mkRandomState.cc: /* ** Translation */
                #now.T.x = random_range( info->lo[X], info->hi[X]);
                #now.T.y = random_range( info->lo[Y], info->hi[Y]);
                #now.T.z = random_range( info->lo[Z], info->hi[Z]);
                trans = []
                if relTrans:
                    trans.append(round(random.uniform(0,1)*xrange,3))
                    trans.append(round(random.uniform(0,1)*yrange,3))
                    trans.append(round(random.uniform(0,1)*zrange,3))
                    #print "in relTrans: trans=", trans
                else:
                    trans.append(round(random.uniform(0,1)*xrange,3) + info_loX)
                    trans.append(round(random.uniform(0,1)*yrange,3) + info_loY)
                    trans.append(round(random.uniform(0,1)*zrange,3) + info_loZ)
                    #print "in else: trans=", trans
                self.trans = trans
                m.state['origin'] = m.mol_cen
                m.state['trans'] = trans
                #try to build torTree here to use to update torsions:
                if not hasattr(m, 'torTree'):
                    m.LPO.write("JUNK.pdbqt")
                m.torTree.setTorsionAngles(dihe) #@@ doesn't change coords!!
                m.state['dihe'] = dihe
                # /* ** Quaternion angular displacement  adapted based on:
                # /*This should produce a uniformly distributed quaternion, according to
                #**  Shoemake, Graphics Gems III.6, pp.124-132, "Uniform Random Rotations",
                #**  published by Academic Press, Inc., (1992) */
                r1 = random.uniform(0,1) #constants.h:local_random: >=0. and <1.0 @@??1??
                random_sign = +1
                if r1<0.5: random_sign = -1 #from constants.h
                PI = math.pi
                TWOPI = 2. * PI #TWOPI = 6.28318530717958647692  from autocomm.h
                # randomQuat from qmultiply.cc
                #@@#t1 = random.uniform(0,TWOPI)  #t1 = genunf(0., TWOPI)<-line 298 qmultiply.cc
                quat = []
                #r1 = random_sign * math.sqrt(1. - x0)
                x0 = random.uniform(0,1) #constants.h:local_random: >=0. and <1.0
                r1 =  random_sign * math.sqrt(1. - x0)
                t1 = TWOPI * random.uniform(0, 1)  #t1 = genunf(0., TWOPI)<-line 298 qmultiply.cc
                #now.Q.x
                quat.append(math.sin(t1)*r1) #q.x 'strict Shoemake version'
                #now.Q.y
                quat.append(math.cos(t1)*r1) #q.y @@lookup python function equivalents
                r2 = random_sign * math.sqrt(x0)
                t2 = TWOPI * random.uniform(0.,1.)
                #now.Q.z
                quat.append(math.sin(t2)*r2) #Q.z @@lookup python function equivalents
                #now.Q.w
                quat.append(math.cos(t2)*r2) #Q.w @@lookup python function equivalents
                self.quat = quat
                m.state['quat'] = quat
                # CONVERT TO AXISANGLE in order to use stoc:
                retval = [] #to hold x,y,z,ang
                x,y,z,w  = quat
                # TODO handle big W! Singularities line219 qmultiply.cc
                #check for identity
                if abs(w)>1.001: raise 'too big W' #?@@?
                if w>1.: w = 1.
                if w<-1.: w = -1. #if abs(w)==1 return sidentityAxisAngle==[1.,0.,0.,0.]
                angle = 2. * math.acos(w)
                inv_sin_half_angle = 1./math.sin(angle/2.)
                nx = inv_sin_half_angle * x
                ny = inv_sin_half_angle * y
                nz = inv_sin_half_angle * z
                ##ang = WrpModRad(angle)
                if angle>math.pi:
                    angle = angle -(2*math.pi)
                elif angle< -1*math.pi:
                    angle = angle +(2*math.pi)
                mag3 = math.sqrt(x**2+y**2+z**2)
                axisangle =[]
                if mag3>0.0000001:
                    inv_mag3 = 1./mag3
                    axisangle.append(inv_mag3*nx)
                    axisangle.append(inv_mag3*ny)
                    axisangle.append(inv_mag3*nz)
                    axisangle.append(angle*360./(2*math.pi))
                m.state['axisangle'] = axisangle
                #two steps:
                #1. use a place-holder axis-angle dummy to create Conformation
                #with appropriate torsion-angles
                #2.then generate random orientation using bits from autodock
                newState = Conformation(m, m.mol_cen, trans, axisangle, dihe) #@@
                m.stoc.applyState(newState)
                # from qmultiply.cc: set orientation
                #/*15 adds, 9 multiplies*/
                tx  = x+x;
                ty  = y+y;
                tz  = z+z;

                twx = w*tx;
                omtxx = 1. - x*tx;
                txy = y*tx;
                txz = z*tx;

                twy = w*ty;
                tyy = y*ty;
                tyz = z*ty;

                twz = w*tz;
                tzz = z*tz;

                r11 = 1. - tyy - tzz;
                r12 =      txy - twz;
                r13 =      txz + twy;
                r21 =      txy + twz;
                r22 = omtxx    - tzz;
                r23 =      tyz - twx;
                r31 =      txz - twy;
                r32 =      tyz + twx;
                r33 = omtxx    - tyy;
                newCrds = []
                ats = m.allAtoms
                tcoord = m.allAtoms.coords[:]
                Tx,Ty,Tz = trans
                for ix,at  in enumerate(ats):
                    tcoord = at.coords
                    tmpx = (tcoord[0])*r11 + (tcoord[1])*r21 + (tcoord[2])*r31 + Tx;
                    tmpy = (tcoord[0])*r12 + (tcoord[1])*r22 + (tcoord[2])*r32 + Ty;
                    tmpz = (tcoord[0])*r13 + (tcoord[1])*r23 + (tcoord[2])*r33 + Tz;
                    newCrds.append((tmpx,tmpy,tmpz))
                m.allAtoms.updateCoords(newCrds, ind=m.orig_coord_index)
                #remove all original bonds and set hasBonds to 0 on all levels
                #@@ REMEMBER WHICH ARE activeTors and not
                aBD = activeBondDict = {}
                rbm = RotatableBondManager(m, allowed_bonds='backbone')  #5/2014
                for b in m.allAtoms.bonds[0]:
                    v1 = b.possibleTors
                    v2 = b.activeTors
                    try:
                        aBD[(b.atom1, b.atom2)] = (v1,v2)
                    except:
                        aBD[(b.atom1.name, b.atom2.name)] = (v1,v2)
                del(m.allAtoms.bonds) #@@ REMOVE ALL BONDS AFTER BUILDING DICT
                m.hasBonds=0
                for c in m.chains: c.hasBonds=0
                for r in m.chains.residues: r.hasBonds=0
                for a in m.allAtoms:
                    a.bonds = BondSet([])
                    a.hasBonds = 0
                m.buildBondsByDistance()
                #do it again here??
                rbm = RotatableBondManager(m, allowed_bonds='backbone') #5/2014
                newLen = len(m.allAtoms.bonds[0])
                new_d = {}
                for a in m.allAtoms: new_d[a]=set(a.bonds.getAtoms())
                ok = True
                for a in m.allAtoms:
                    if orig_d[a]!=new_d[a]:
                        ok = False
                if ok:

                    break
            if ok:
                #reset bond activity:
                for k, v in aBD.items(): # k = (b.atom1, b.atom2);
                    for b in k[0].bonds:
                        found = 0
                        b.possibleTors = 0
                        b.activeTors = 0
                        if b.neighborAtom(k[0])==k[1]:
                            found = 1
                            b.possibleTors, b.activeTors = v

                #might have been a pdb or ? file so this won't do:
                m.parser.write_with_new_coords(newCrds, outputfilename)
                self.ok = 1
            else:
                self.ok = 0



class RotatableBondManager:
    """
    flags bonds of ligand molecules for possible rotability and active
    rotability.  [for backwards compatibility] The flags are 'possibleTors'
    and 'activeTors'
    """


    def __init__(self, molecule, allowed_bonds=['backbone'],
                    root='auto',  debug=False,
                    check_for_fragments=False,
                    bonds_to_inactivate='',
                    detectAll=False):

        self.detectAll = detectAll
        allowed_bond_list = []
        if allowed_bonds:
            # allowed_bond_list = string.split(allowed_bonds,'_')
            allowed_bond_list = allowed_bonds.split('_')
        #set up flags:
        allow_amide_torsions = 'amide' in allowed_bond_list
        molecule.has_amide = allow_amide_torsions
        allow_backbone_torsions = 'backbone' in allowed_bond_list
        molecule.has_backbone = allow_backbone_torsions
        allow_guanidinium_torsions = 'guanidinium' in allowed_bond_list
        molecule.has_guanidinium = allow_guanidinium_torsions
        allow_all_torsions = 'all' in allowed_bond_list
        self.molecule = molecule
        self.debug = debug
        self.__classifyBonds(molecule.allAtoms, allow_guanidinium_torsions)
        self.set_peptidebackbone_torsions(allow_backbone_torsions)
        self.set_amide_torsions(allow_amide_torsions)
        #self.set_guanidinium_torsions(allow_guanidinium_torsions)
        #NEW 11/22/2004
        if check_for_fragments:
            self.detect_bonded_fragments()

        if root=='auto':
            self.autoroot()
        else:
            self.setroot(int(root))
        if len(bonds_to_inactivate):
            bnds_list = map(int,string.split(bonds_to_inactivate,'_'))
            #???
            #molecule.has_amide = 'amide' not in bnds_list
            #molecule.has_guanidinium = 'guanidinium' not in bnds_list
            #molecule.has_backbone = 'backbone' not in bnds_list
            #molecule.has_active = 'all' not in bnds_list
            for i in range(len(bnds_list)/2):
                ind1 = bnds_list[i*2]
                ind2 = bnds_list[i*2+1]
                self.toggle_torsion(ind1, ind2)


    def __classifyBonds(self, atoms, allow_guanidinium_torsions):
        mols = atoms.top.uniq()
        #check that all are in the same molecule, for now
        assert  len(mols)==1
        mol = mols[0]
        if hasattr(mol,'processed_bonds'):
            return
        #check whether these atoms already have bonds
        if not len(atoms.bonds[0]):
            mol.buildBondsByDistance()
        ADBC = self.ADBC = AutoDockBondClassifier(detectAll=self.detectAll)
        dict =self.dict = ADBC.classify(mol.allAtoms.bonds[0])
        #turn everything off
        mol.allAtoms.bonds[0].possibleTors = 0
        mol.allAtoms.bonds[0].activeTors = 0
        mol.allAtoms.bonds[0].leaf = 0
        mol.allAtoms.bonds[0].incycle = 0
        # restore appropriate categories:
        if len(dict['leaf']):
            dict['leaf'].leaf = 1
        aromBnds = dict['aromatic']
        if len(dict['cycle']):
            dict['cycle'].incycle = 1
        if len(dict['rotatable']):
            for b in dict['rotatable']:
                #check for mistake
                if b.incycle:
                    b.possibleTors = 0
                    b.activeTors = 0
                    dict['rotatable'].remove(b)
            dict['rotatable'].possibleTors = 1
            dict['rotatable'].activeTors = 1
            #NB: this should set amide and ppbb
        #mol.guanidiniumbnds = BondSet()
        #self.set_guanidinium_torsions(allow_guanidinium_torsions)
        for b in dict['guanidinium']:
            if not b.incycle and not b.leaf and b.bondOrder==1:
                b.possibleTors = 1
                b.activeTors = allow_guanidinium_torsions
        #if len(dict['guanidinium']): and not allow_guanidinium_bonds:
        #    #as per gmmm 4/14:
        #    mol.guanidiniumbnds = dict['guanidinium']
        #    for b in dict['guanidinium']:
        #        if self.debug: print b, ' ', b.bondOrder
        #        if not b.incycle and not b.leaf and b.bondOrder==1:
        #            if self.debug: print 'set ', b,'.possible/active to 1/0'
        #            b.possibleTors = 1
        #            b.activeTors = 0
        #set keywords for this molecule
        #NB: 'torscount' must be adjusted if specific torsions are inactivated
        rotatable_bnds = mol.allAtoms.bonds[0].get(lambda x: x.activeTors==1)
        possible_tors_bnds = mol.allAtoms.bonds[0].get(lambda x: x.possibleTors==1)
        possible_tors_ct = 0
        mol.possible_tors_bnds = BondSet()
        if possible_tors_bnds:
            possible_tors_ct = len(possible_tors_bnds)
            mol.possible_tors_bnds = possible_tors_bnds
        mol.torscount = 0
        mol.possible_tors = possible_tors_ct
        if rotatable_bnds:
            #THIS IS WRONG!!rotatable includes too many
            #mol.torscount = len(rotatable_bnds)
            #torscount changes, possible_tors_ct doesn't
            mol.torscount = possible_tors_ct

        #mol.TORSDOF = mol.torscount - len(dict['hydrogenRotators'])
        mol.TORSDOF = possible_tors_ct - len(dict['hydrogenRotators'])
        mol.possible_tors = possible_tors_ct

        mol.amidebnds = dict['amide']
        mol.ppbbbnds = dict['ppbb']
        mol.guanidiniumbnds = dict['guanidinium']
        mol.processed_bonds = 1


    def set_torsions(self, mol, bnds, flag):
        ct = 0
        for b in bnds:
            if not b.possibleTors:
                continue
            if flag:
                if not b.activeTors:
                    b.activeTors = 1
                    mol.torscount = mol.torscount + 1
                    ct = ct + 1
            else:
                if b.activeTors==1:
                    b.activeTors = 0
                    mol.torscount = mol.torscount - 1
                    ct = ct - 1
        return ct


    def set_amide_torsions(self, flag):
        """
        set rotability of all amide bonds to flag: True/False
        """
        assert flag in [True, False, 1, 0]
        mol = self.molecule
        #print "before call to set_torsions: mol.amidebnds.activeTors=", mol.amidebnds.activeTors
        ct = self.set_torsions(mol, mol.amidebnds, flag)


    def set_peptidebackbone_torsions(self, flag):
        """
        set rotability of all peptidebackbone bonds to flag: True/False
        """
        assert flag in [True, False, 1, 0]
        mol = self.molecule
        ct = self.set_torsions(mol, mol.ppbbbnds, flag)


    def set_guanidinium_torsions(self, flag):
        """
        set rotability of all guanidinium bonds to flag: True/False
        """

        assert flag in [True, False, 1, 0]
        mol = self.molecule
        ct = self.set_torsions(mol, mol.guanidiniumbnds, flag)



    def set_all_torsions(self, flag):
        """
        set rotability of all rotatable bonds to flag: True/False
        """
        assert flag in [True, False, 1, 0]
        mol = self.molecule
        rotatable = list(filter(lambda x: x.possibleTors==1, mol.allAtoms.bonds[0]))
        if len(rotatable):
            ct = self.set_torsions(mol, rotatable, flag)



    def toggle_torsion(self, ind1, ind2, use_tt_ind=0):
        """
        flip the rotability of the bond
        between atoms with indices ind1 and ind2
        use_tt_ind is set when limiting torsion using a built-on-the-fly
        rotabilityManager
        """
        mol = self.molecule
        at1 = mol.allAtoms[ind1]
        at2 = mol.allAtoms[ind2]
        if use_tt_ind:
            at1 = mol.allAtoms.get(lambda x: hasattr(x, 'tt_ind') and x.tt_ind==ind1)
            if not at1:
                return "ERROR"
            at1 = at1[0]
            at2 = mol.allAtoms.get(lambda x: hasattr(x, 'tt_ind') and x.tt_ind==ind2)
            if not at2:
                return "ERROR"
            at2 = at2[0]
        bnds = AtomSet([at1,at2]).bonds[0]
        if not len(bnds):
            return
        bnd = bnds[0]
        if bnd.possibleTors:
            #print "calling set_torsions with ", bnd, " and flag=", abs(bnd.activeTors-1)
            self.set_torsions(mol, [bnd], abs(bnd.activeTors-1))
        else:
            return


    def detect_bonded_fragments(self):
        atoms = self.molecule.allAtoms
        #first: build a list of bonded fragments:
        l = []   # the list to hold dicts of bonded fragments
        fb = []  # keep track of which bonds are found
        #start with unique list of bonds 'bnds'
        bond_dict = {}
        for b in atoms.bonds[0]:
            bond_dict[b] = 1
        bnds = bond_dict.keys()
        #first pass: try to connect pieces
        for b in bnds:
            ind1 = b.atom1
            ind2 = b.atom2
            found = b in fb
            if not found:
                for d in l:
                    d_keys = d.keys()
                    if ind1 in d_keys:
                        d[ind2] = 1
                        fb.append(b)
                    elif ind2 in d_keys:
                        d[ind1] = 1
                        fb.append(b)
                if not b in fb:
                    #start new fragment
                    l.append({ind1:1, ind2:1})
                    fb.append(b)
        #now to try to merge fragments, use lists of keys
        key_list = []
        for dict in l:
            key_list.append(dict.keys())
        #check each list of keys against following lists
        #if there are any duplications, merge current
        #into the following one..
        for i in range(len(key_list)):
            #check this list
            kl = key_list[i]
            found = 0
            #...against the each of the subsequent ones
            for j in range(i+1, len(key_list)):
                jl = key_list[j]
                #....check each entry in subsequent one
                #.........against this list 'kl'
                for entry in jl:
                    if entry in kl:
                        #if a match
                        #...merge this dict into jth one..
                        l[j].update(l[i])
                        #.....reset this dict to {}
                        l[i] = {}
                        #.......set found flag
                        found = 1
                        #..........update jth list of keys
                        key_list[j]= l[j].keys()
                        #............skip rest of jl
                        break
                if found:
                    #.................and skip rest of key_list
                    break
        #now find an index in the largest fragment
        max_ind = 0
        largest_fragment = {}
        ct = 0
        for dict in l:
            len_d = len(dict)
            if len_d >0:
                ct = ct + 1
            if len_d > max_ind:
                max_ind = len(dict)
                largest_fragment = dict
                #the keys are ATOMS!
                #print z_keys
        if ct>0:
            #build an atomset of the keys
            self.molecule.largest_fragment = AtomSet(largest_fragment.keys())


    def autoroot(self, verbose=0):
        """
        autoRoot
        """
        mol = self.molecule
        #clear old root
        if hasattr(mol, 'ROOT') and hasattr(mol.ROOT, 'rnum0'):
            delattr(mol.ROOT, 'rnum0')
        if hasattr(mol, 'autoRoot'):
            mol.ROOT = mol.autoRoot
            mol.ROOT.rnum0 = 0
            return
        if len(mol.chains)>1:
            return  "AutoRoot not implemented for molecules with >1 chain"
        if hasattr(mol, 'largest_fragment'):
            atoms_to_check = mol.largest_fragment
        else:
            atoms_to_check = mol.allAtoms
        #mol.bestBranch = len(mol.allAtoms)
        mol.bestBranch = len(atoms_to_check)
        bestList = []
        #for item in mol.allAtoms:
        for item in atoms_to_check:
            if not hasattr(item, 'bonds'):
                continue
            if len(item.bonds)==1 and item.bonds[0].leaf:
                continue
            if hasattr(item,'leaf'):
                continue
            item.maxbranch = 0
            for b in item.bonds:
                nxtAtom = b.atom1
                if nxtAtom==item:
                    nxtAtom = b.atom2
                if not b.leaf:
                    thistree = mol.subTree(item, nxtAtom, atoms_to_check)
                    #thistree = mol.subTree(item, nxtAtom, mol.allAtoms)
                    thisbranch = len(thistree)
                    if thisbranch>item.maxbranch:
                        item.maxbranch = thisbranch
            #bestList holds all current best choices for Root..
            if item.maxbranch < mol.bestBranch:
                bestList = []
                bestList.append(item)
                mol.bestBranch = item.maxbranch

            if item.maxbranch==mol.bestBranch and item not in bestList:
                bestList.append(item)
        if len(bestList)>1:
            foundCycle = 0
            for at in bestList:
                at.cycleatom = 0
                for b in at.bonds:
                    if hasattr(b, 'incycle') and b.incycle:
                        at.cycleatom = 1
                        continue
            for at in bestList:
                if at.cycleatom:
                    mol.ROOT = at
                    mol.autoRoot = at
                    mol.ROOT.rnum0 =0
                    foundCycle = 1
                    break
            #if bestList had a cycle atom, it's been set to root..if NOT:
            if not foundCycle:
                mol.autoRoot = bestList[0]
                mol.ROOT = bestList[0]
                mol.ROOT.rnum0 =0
        #if ties for possible root, use first entry in bestRoot list...
        elif len(bestList):
            mol.autoRoot = bestList[0]
            mol.ROOT = bestList[0]
            mol.ROOT.rnum0 =0
        else:
            #if the list is empty, use the first atom
            # this is a correction added for 'HF' ligand
            mol.autoRoot = mol.allAtoms[0]
            mol.ROOT = mol.autoRoot
            mol.ROOT.rnum0 =0
        return mol.ROOT


    def setroot(self, index):
        """
        setroot to atom at index
        """
        #if autoRoot has already been set, remove its 'rnum0'
        mol = self.molecule
        if hasattr(mol, 'autoRoot') and hasattr(mol.autoRoot, 'rnum0'):
            delattr(mol.autoRoot, 'rnum0')
        if hasattr(mol, 'ROOT') and hasattr(mol.ROOT, 'rnum0'):

            delattr(mol.ROOT, 'rnum0')
        #set the root to indicated atom
        mol.ROOT = mol.allAtoms[index]
        #set the root's rnum0 field
        mol.ROOT.rnum0 = 0
        return


    def limit_torsions(self, numTors, type='fewest'):
        """
        numTors, type='fewest'
        limit torsions to a specified number, numTors
        #changed per GMM request: invert logic of type here
        previously 'type' meant TOGGLE torsions moving 'type' atoms:
        now 'type' means KEEP torsions moving type atoms:
        options are
                'fewest' atoms which is the default
              or 'most'
        """
        mol = self.molecule
        if not hasattr(mol, 'ROOT'):

            return 'ERROR'
        if not hasattr(mol, 'torTree') or \
                not hasattr(mol.allAtoms[0], 'tt_ind'):
            mol.torTree = TorTree(mol.parser, mol.ROOT)
        self._setTorsions(numTors, type)


    def _setTorsions(self, numTors, type):
        #print "in _setTorsions ", numTors, " and type ", type
        assert type in ['fewest', 'most'], 'unknown torsion type'
        mol = self.molecule
        torsionMap = mol.torTree.torsionMap
        tNum = len(torsionMap)
        if numTors>tNum:
            numTors = tNum
        #set up rangeList to implement which kind of torsion
        # to turn off or on
        #fewest uses range 0,1,2,.... (from beginning toward end)
        #most uses range -1, -2,-3,.... (from end toward beginning)
        if type=='fewest':
            rangeList = range(numTors)
        else:
            rangeList = []
            for k in range(1, numTors+1):
                rangeList.append(-k)
        #turn them all off
        ct = self.set_torsions(mol, mol.allAtoms.bonds[0], 0)
        #print "ct = ", ct, " and mol.torscount=", mol.torscount
        #torsionMap = mol.torTree.torsionMap
        #for i in range(tNum):
        #    node = torsionMap[i]
        #    ind1, ind2 = node.bond
        #    self.toggle_torsion(ind1, ind2, use_tt_ind=1)

            #b = mol.allAtoms.get(lambda x, node=node: x.tt_ind in node.bond).bonds[0][0]
            #b.activeTors = 0
            #mol.torscount -= 1
        #turn on the right number at correct end
        for i in rangeList:
            node = torsionMap[i]
            ind1, ind2 = node.bond
            self.toggle_torsion(ind1, ind2, use_tt_ind=1)
            #b = mol.allAtoms.get(lambda x, node=node: x.tt_ind in node.bond).bonds[0][0]
            #b.activeTors = 1
            #mol.torscount += 1

