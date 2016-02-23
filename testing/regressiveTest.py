import sys
import os.path
import numpy as np
from rmgpy.chemkin import loadSpeciesDictionary
from rmgpy.species import Species
from rmgpy.tools.plot import GenericData, GenericPlot
import cantera as ct
from rmgpy.tools.plot import *

#for time timepoint1, it returns the closest value in list t2 as well as the index of this value
def getNearestTime(timepoint1, t2):
    timeDiff=1e9 #initialize to large number
    for index, timepoint2 in enumerate(t2):
        if abs(timepoint2-timepoint1) < timeDiff:
            timeDiff = abs(timepoint2-timepoint1)
        else:
            return (index-1, t2[index-1])
    else: return (index, t2[index])

"""
This function returns True if the two given curves are similar enough within tol. Otherwise returns False.

t1: time/domain of standard curve we assume to be correct
y1: values of standard curve, usually either temperature in (K) or log of a mol fraction
t2: time/domain of test curve
y2: values of test curve, usually either temperature in (K) or log of a mol fraction

The test curve is first synchronized to the standard curve using geatNearestTime function. We the calculate the value of
(y1-y2')^2/y1^2, giving us a normalized difference for every point. If the average value of these differences is less
than tol, we say the curves are similar.

We choose this criteria because it is compatible with step functions we expect to see in ignition systems.
"""
def curvesSimilar(t1, y1, t2, y2, tol):
    #Make synchornized version of t2,y2 called t3,y3.
    t3=[]
    y3=[]
    for timepoint1 in t1:
        (index, timepoint2)=getNearestTime(timepoint1, t2)
        t3.append(timepoint2)
        y3.append(y2[index])

    #get R^2 value equivalent:
    normalizedError=[(y1[x]-y3[x])**2/y1[x]**2 for x in range(len(y1))]

    if sum(normalizedError)/len(normalizedError)>tol:
        return False
    else: return True

#Finds the ignition delay based on position of max dT/dt NOT TESTED YET!, technically you could probably replace T
#with P or even OH concentration
def findIgnitionDelay(t, T):
    forTDiff=[T[x+1]-T[x] for x in range(len(t)-1)]
    forTimeDiff=[t[x+1]-t[x] for x in range(len(t)-1)]
    dTdt=[T[x]/t[x] for x in range(len(forTDiff))]

    maxIndex=dTdt.find(max(dTdt))

    return t[maxIndex]

def getNameFromRMGDict(smilesList, speciesDictPath, type="dict"):
    """
    Returns the RMG species names from a list of the SMILES.

    Args:
        smilesList: list of SMIlES for species of interest
        speciesDictPath: path to the RMG-dictionary
        type: "list" of "dict", type of output requested

    Returns: A list of RMG names in the same order or a dictionary of with SMILES keys and RMG name values

    """
    speciesDict=loadSpeciesDictionary(speciesDictPath)
    bathGases={"Ar": "Ar", "NN": "N2", "He": "He", "Ne": "Ne"}

    if type=='list':
        finalList=[]
        for smiles in smilesList:
            if smiles in bathGases:
                finalList.append(bathGases[smiles])
            else:
                species1=Species().fromSMILES(smiles)
                for name, species2 in speciesDict.iteritems():
                    if species2.isIsomorphic(species1):
                        finalList.append(name)
                        break
                else: raise KeyError("The major species {0} does not appear in the dictionary at {1}".format(smiles, speciesDictPath))
        return finalList

    elif type=='dict':
        finalDict={}
        for smiles in smilesList:
            if smiles in bathGases:
                finalDict[smiles]=bathGases[smiles]
            else:
                species1=Species().fromSMILES(smiles)
                for name, species2 in speciesDict.iteritems():
                    if species2.isIsomorphic(species1):
                        finalDict[smiles]=name
                        break
                else: raise KeyError("The major species {0} does not appear in the dictionary at {1}".format(smiles, speciesDictPath))
        return finalDict
    else: raise Exception("The 'type' argument is not compatible. Please enter either 'list' or 'dict'")

#######################################################################################
#Currently only implementing simple inputs I'm familiar with. New version should probably include more:
#eg, isothermal reaction type, option for mass fraction or conc of species, use species instead of SMILES?
class Condition:
    """
    This class organizes the inputs needed for a cantera simulation

    ======================= ====================================================
    Attribute               Description
    ======================= ====================================================
    `reactorType`           A string of the cantera reactor type. List of supported types below:
        IdealGasReactor: A constant volume, zero-dimensional reactor for ideal gas mixtures
        IdealGasConstPressureReactor: A homogeneous, constant pressure, zero-dimensional reactor for ideal gas mixtures

    `reactionTime`          A float giving the reaction time in seconds
    `molFrac`               A dictionary giving the initial mol Fractions. Keys are SMILES strings and values are floats

    To specifiy the system for an ideal gas, you must exactly 2 of the following 3 parameters:
    `T0`                    A float giving the initial temperature in K
    'P0'                    A float giving the initial pressure in Pa
    'V0'                    A float giving the initial specific volume in m^3/kg
    ======================= ====================================================


    """
    def __init__(self, reactorType, reactionTime, molFrac, T0=-1, P0=-1, V0=-1):
        self.reactorType=reactorType
        self.reactionTime=float(reactionTime)
        #Normalize initialMolFrac if not already done:
        if sum(molFrac.values())!=1.00:
            total=sum(molFrac.values())
            for species, value in molFrac.iteritems():
                molFrac[species]= value / total
        self.molFrac=molFrac
        self.T0=float(T0)
        self.P0=float(P0)
        self.V0=float(V0)

        #check to see that one of T0, P0, and V0 is unspecificed (left at default of -1)
        checkList=[self.T0, self.P0, self.V0]
        total=0
        for value in checkList:
            if value==-1: total+=value
        assert total==-1, "For a condition, one of the state variables: T0, P0, or V0 must be left out, otherwise the system is overspecified"


    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string="Condition("
        string += 'reactorType="{0}", '.format(self.reactorType)
        string += 'reactionTime={:0.10f}, '.format(self.reactionTime)
        string += 'molFrac={0}, '.format(self.molFrac.__repr__())
        if self.T0 != -1: string += 'T0={:0.10f}, '.format(self.T0)
        if self.P0 != -1: string += 'P0={:0.10f}, '.format(self.P0)
        if self.V0 != -1: string += 'V0={:0.10f}, '.format(self.V0)
        string = string[:-2] + ')'
        return string

    def __str__(self):
        """
        Return a string representation of the condition.
        """
        string=""
        string += 'Reactor Type: {0}, '.format(self.reactorType)
        string += 'Reaction Time: {:0.10f}, '.format(self.reactionTime)
        if self.T0 != -1: string += 'T0: {:0.10f}, '.format(self.T0)
        if self.P0 != -1: string += 'P0: {:0.10f}, '.format(self.P0)
        if self.V0 != -1: string += 'V0: {:0.10f}, '.format(self.V0)
        string += 'Initial Mole Fractions: {0}, '.format(self.molFrac.__repr__())
        return string

def generateAllConditions(reactorType, reactionTime, molFracList, Tlist=None, Plist=None, Vlist=None):
    """
    Creates a list of conditions making every combination from the lists provided. The user can input lists of mole
    fractions(dictionaries), temperatures, pressures, and volumes.
    """
    conditions=[]

    if Tlist is None:
        for molFrac in molFracList:
            for P in Plist:
                for V in Vlist:
                    conditions.append(Condition(reactorType, reactionTime, molFrac, P0=P, V0=V))

    elif Plist is None:
        for molFrac in molFracList:
            for T in Tlist:
                for V in Vlist:
                    conditions.append(Condition(reactorType, reactionTime, molFrac, T0=T, V0=V))

    elif Vlist is None:
        for molFrac in molFracList:
            for T in Tlist:
                for P in Plist:
                    conditions.append(Condition(reactorType, reactionTime, molFrac, T0=T, P0=P))

    else:
        for molFrac in molFracList:
            for T in Tlist:
                for P in Plist:
                    for V in Vlist:
                        conditions.append(Condition(reactorType, reactionTime, molFrac, T0=T, P0=P, V0=V))
    return conditions

#######################################################################################
class ObservablesTestCase:
    """
    We use this class to run regressive tests

    ======================= ====================================================
    Attribute               Description
    ======================= ====================================================
    Inputted attributs:
    'title'                 A string describing what is tested and the major species: e.g. "Isobutanol P-dep"
    `newDir`                A string path to chem.inp and species.txt of new model
    `oldDir`                A string path to chem.inp and species.txt of new model
    `conditions`            A list of the :class: 'Condition' objects describing reaction conditions
    `majorSpeciesSmiles`    A list of SMIlES for the major species to compare in the regressive tests
    'exptData'              GenericData objects

    Generated Attributes:
    'results'
    'oldDict'
    'newDict'
    ======================= ====================================================


    """
    def __init__(self, title, newDir, oldDir, conditions, majorSpeciesSmiles, exptData=None):
        self.title=title
        self.newDir=newDir
        self.oldDir=oldDir
        self.conditions=conditions
        self.majorSpeciesSmiles=majorSpeciesSmiles
        self.exptData=exptData

        self.results=[]
        self.oldDict=getNameFromRMGDict(majorSpeciesSmiles, os.path.join(oldDir, 'species_dictionary.txt'), type="dict")
        self.newDict=getNameFromRMGDict(majorSpeciesSmiles, os.path.join(newDir, 'species_dictionary.txt'), type="dict")
        self.oldSpeciesList=getNameFromRMGDict(majorSpeciesSmiles, os.path.join(oldDir, 'species_dictionary.txt'), type="list")
        self.newSpeciesList=getNameFromRMGDict(majorSpeciesSmiles, os.path.join(newDir, 'species_dictionary.txt'), type="list")

    def __str__(self):
        """
        Return a string representation of the species, using the label'.
        """
        return 'Observables Test Case: {0}'.format(self.title)

    #need to add a bunch of if statements to change based on condition options
    def compare(self, plot=False):
        """
        Compare an old and new model
        """
        conditionData = {}
        for modelName in ['old','new']:
                
            # Set the model, speciesList, and dictionary initial mole fractions
            if modelName == 'old':
                model = ct.Solution(os.path.join(self.oldDir,'chem.cti'))
                speciesList = self.oldSpeciesList
                speciesDict = self.oldDict
            else:
                model = ct.Solution(os.path.join(self.newDir,'chem.cti'))
                speciesList = self.newSpeciesList
                speciesDict = self.newDict
            
            # Ignore Inerts
            inertList = ['Ar','He','NN','Ne']
                    
            conditionData[modelName] = self.runSimulations(model, speciesList, speciesDict)
        if plot:
            for i in range(len(conditionData['old'])):
                time, dataList = conditionData['old'][i]
                speciesData = [data for data in dataList if data.species not in inertList]
                oldSpeciesPlot = SimulationPlot(xVar=time, yVar=speciesData, ylabel='Mole Fraction')

                time, dataList = conditionData['new'][i]
                speciesData = [data for data in dataList if data.species not in inertList]
                newSpeciesPlot = SimulationPlot(xVar=time, yVar=speciesData, ylabel='Mole Fraction')

                # Name after the index of the condition
                # though it may be better to name it after the actual conditions in T, P, etc
                oldSpeciesPlot.comparePlot(newSpeciesPlot,filename='simulation_condition_{0}.png'.format(i+1))

    def runSimulations(self, model, speciesList, speciesDict):
        """
        Run a selection of conditions in Cantera and return
        generic data objects containing the time, pressure, temperature,
        and mole fractions from the simulations.
        """
        
        conditionData = []
        for condition in self.conditions:
            # Set the mole fractions
            molFrac = {}
            for smiles, value in condition.molFrac.iteritems():
                molFrac[speciesDict[smiles]]=value
                        
            # Set Cantera simulation conditions
            if condition.T0==-1:
                model.PVX = condition.P0, condition.V0, molFrac
            elif condition.P0==-1:
                model.TVX = condition.T0, condition.V0, molFrac
            else:
                model.TPX= condition.T0, condition.P0, molFrac

            #Set Cantera Reactor Type
            if condition.reactorType=="IdealGasReactor": canteraReactor=ct.IdealGasReactor(model)
            elif condition.reactorType=="IdealGasConstPressureReactor": canteraReactor=ct.IdealConstPressureGasReactor(model)
            else: raise Exception("reactorType is not compatiable. Please set to IdealGasReactor or IdealConstPressureGasReactor")
            canteraSimulation=ct.ReactorNet([canteraReactor])

            time = 0.0
            
            # Initialize the variables to be saved
            times = np.zeros(100)
            temperature = np.zeros(100, dtype=np.float64)
            pressure = np.zeros(100, dtype=np.float64)
            speciesData = np.zeros((100,len(speciesList)),dtype=np.float64)

            #print condition
            #print('%10s %10s %10s %14s' % ('t [s]','T [K]','P [Pa]','u [J/kg]'))
            
            
            # Run the simulation over 100 time points
            for n in range(100):
                time += condition.reactionTime/100
                canteraSimulation.advance(time)
                times[n] = time * 1e3  # time in ms
                temperature[n] = canteraReactor.T
                pressure[n] = canteraReactor.thermo.P
                speciesData[n,:] = canteraReactor.thermo[speciesList].X
                
            # Resave data into generic data objects
            time = GenericData(label = 'Time', 
                               data = times,
                               units = 'ms')
            temperature = GenericData(label='Temperature',
                                      data = temperature,
                                      units = 'K')
            pressure = GenericData(label='Pressure',
                                      data = pressure,
                                      units = 'Pa')
            dataList = []
            dataList.append(temperature)
            dataList.append(pressure)
            
            for smiles in speciesDict.keys():
                speciesLabel = speciesDict[smiles]
                index = speciesList.index(speciesLabel)
                genericData = GenericData(label=speciesLabel,
                                          species = smiles,
                                          data = speciesData[:,index])
                dataList.append(genericData)
            
            conditionData.append((time,dataList))
            
        return conditionData
