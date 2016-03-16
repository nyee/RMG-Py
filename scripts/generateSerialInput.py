#!/usr/bin/env python
# encoding: utf-8

"""
This script generates a input.py using species concentrations from a previous RMG job.
In its first iteration, it copies the previous input file, but adds the old mechanism as a seed mechanism.

Then it takes the output from the LAST simulation done by RMG and uses this as a base for the new mole fractions to
put into the next RMG job.

It will then also add the additive denoted by additiveSMILES to the initial mole fraction of the new RMG job

    $ python generateSerialInput.py name /path/to/oldWorkDir /path/to/newWorkDir additiveSMILES stream2Mass

"""
from rmgpy.tools.plot import parseCSVData
from rmgpy.chemkin import getSpeciesIdentifier
from rmgpy.chemkin import loadSpeciesDictionary
from rmgpy.species import Species
import os.path
import argparse
import re
import copy as cp

"""
Still need to set up submission script that runs first RMG job, this script, importChemkinLibrary script, and then
second RMG job in sequence.

"""
################################################################################
if __name__ == '__main__':
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('name', type=str, help='the name you want associated with the job/library made')
    parser.add_argument('oldWorkDir', type=str, help='Path to the old working directory')
    parser.add_argument('newWorkDir', type=str, help='Path to the new working directory')
    parser.add_argument('additiveSMILES', type=str, help='SMILES for the additive being put into the second reactor')
    parser.add_argument('stream2Mass', type=float, help='Mass in second stream normalized over mass of first stream')
    #eg, if they are the mass flow rate, use 1 for last argument, if second stream is twice as massive, use 2.

    args = parser.parse_args()

    #Base of reactor block
    # newReactorBlockTemp=["simpleReactor(\n",
    #                   "    temperature= (1900,'K'),\n",
    #                   "    pressure=(.58,'bar'),\n",
    #                   "    initialMoleFractions={\n",
    #                   "    },\n",
    #                   "    terminationConversion={\n",
    #                   "        'methane': 0.0055,\n",
    #                   "    },\n",
    #                   "    terminationTime=(1e6,'s'),\n",
    #                   ")\n",
    # ]

    newReactorBlockTemp=["simpleReactor(\n",
                      "    temperature= (1900,'K'),\n",
                      "    pressure=(.58,'bar'),\n",
                      "    initialMoleFractions={\n",
                      "    },\n",
                      "    terminationTime=(1e6,'s'),\n",
                      ")\n",
    ]


###########Find correct simulation file#################################################################################
    #First find last reactor solver:
    lastSimulation=0
    solverPath=os.path.join(args.oldWorkDir, 'solver')
    for filename in os.listdir(solverPath):
        if filename.endswith(".csv"):
            number=re.sub('simulation\_', '', filename)
            number =int(re.sub('\_.*', '', number))
            if number > lastSimulation:
                lastSimulation=cp.copy(number)


    #Now find the correct file, it should be the solver with the highest number at the end
    lastIteration=0
    lastFile=''
    for filename in os.listdir(solverPath):
        if re.search('simulation\_'+str(lastSimulation)+'\_', filename):
            number = re.sub('.*\_', '', filename)
            number =int(re.sub('\.csv', '', number))
            if number > lastIteration:
                lastIteration=cp.copy(number)
                lastFile=cp.copy(filename)

    #Extract data from file
    time, dataList=parseCSVData(os.path.join(solverPath, lastFile))




########################Write out correct reaction stream###############################################################
    #Import species dictionary
    speciesDict=loadSpeciesDictionary(os.path.join(args.oldWorkDir, "chemkin/species_dictionary.txt"))
    bathGases={"Ar": Species().fromSMILES("[Ar]"),
               "He": Species().fromSMILES("[He]"),
               "Ne": Species().fromSMILES("[Ne]"),
               "N2": Species().fromSMILES("N#N")}
    speciesDict.update(bathGases)

    #mol fraction of stream1
    stream1={} #key is a Species Object, value is the final mol fraction from the simulation
    tol = 1E-10
    for item in dataList:
        #Check that it is a species and not volume or a sensitivity
        if item.species is not None or item.label in bathGases:
            if item.data[-1]>tol:
                stream1[speciesDict[item.label]]=item.data[-1]

    newStream={} #For the final stream that goes into input of second RMG JOB
    additiveSpecies=Species().fromSMILES(args.additiveSMILES)
    additivePresent=False
    additiveName=""

    #convert everything to number of mols in stream1 and stream 2
    stream1Mols={}
    for component in stream1:
        stream1Mols[component]=stream1[component]/component.molecule[0].getMolecularWeight()
    stream2Mols={}
    stream2Mols[additiveSpecies]=args.stream2Mass/additiveSpecies.molecule[0].getMolecularWeight()

    #get total number of mols for new stream:
    total=0
    for component, value in stream1Mols.iteritems():
        total+=value
    for component, value in stream2Mols.iteritems():
        total+=value

    #update newStream with additive
    for species in stream1:
        if species.isIsomorphic(additiveSpecies):
            newStream[species]=(stream1Mols[species]+stream2Mols[additiveSpecies])/total
            additivePresent=True
        else: newStream[species]=stream1Mols[species]/total
    if not additivePresent:
        newStream[additiveSpecies]=stream2Mols[additiveSpecies]/total
        additiveName=getSpeciesIdentifier(additiveSpecies)+"add"
    else: additiveNAme=getSpeciesIdentifier(additiveSpecies)

    #Add species into reactor block
    newReactorBlock=[]
    for line in newReactorBlockTemp:
        newReactorBlock.append(line)
        if re.search("initialMoleFractions",line):
            for species, value in newStream.iteritems():
                if species==additiveSpecies and not additivePresent:
                    newReactorBlock.append("        '"+ additiveName + "': "+ str(value)+ ",\n")
                else:
                    newReactorBlock.append("        '"+ getSpeciesIdentifier(species) + "': "+ str(value)+ ",\n")

#Construct adjList dictionary with correct indentation
    adjListDict={}
    for component in newStream:
        newAdjListList=re.split("\n", component.molecule[0].toAdjacencyList())
        newAdjList=""
        for line in newAdjListList:
            newAdjList+="        "+line+"\n"
        adjListDict[component]=newAdjList


####################Copy in, edit, write out new input file############################################################
    #Copy in the old inputFile:
    inputList=[]
    with open(os.path.join(args.oldWorkDir, 'input.py'), 'rb') as inputFile:
        for line in inputFile:
            inputList.append(line)

    newInputList=[]
    newReactorsInputted=False
    dontCopy=False
    #Start editting the input file
    for line in inputList:
        #change the reaction libraries
        if re.search('reactionLibraries', line):
            newLine=re.sub('\]', ",('" +args.name+ "',True)]", line)
            newInputList.append(newLine)
        #Change the species
        elif re.search('species\(', line):
            dontCopy=True
            for component in newStream:
                newInputList.append("species(\n")
                if component==additiveSpecies and not additivePresent:
                    newInputList.append("    label='"+additiveName+"',\n")
                else:
                    newInputList.append("    label='"+getSpeciesIdentifier(component)+"',\n")
                newInputList.append("    reactive=True,\n")
                newInputList.append('    structure=adjacencyList(\n        """\n')
                newInputList.append(adjListDict[component])
                newInputList.append('        """),\n')
                newInputList.append(")\n\n")
        #Otherwise just copy the same list
        elif re.search("simpleReactor\(", line):
            dontCopy=True
            if not newReactorsInputted:
                newReactorsInputted=True
                for line in newReactorBlock:
                    newInputList.append(line)
        else:
            if not dontCopy: newInputList.append(line)
        #Start copying again after reaching end of any input block
        if re.match('\)', line): dontCopy=False

    #Make a new directory and save new input.py
    if not os.path.exists(args.newWorkDir):
        os.makedirs(args.newWorkDir)
    with open(os.path.join(args.newWorkDir, 'input.py'), 'wb') as outFile:
        for line in newInputList:
            outFile.write(line)
