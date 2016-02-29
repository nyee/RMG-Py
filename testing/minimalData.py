import sys
import math
import numpy as np
import regressiveTest as rt
import cantera as ct
import collections
from rmgpy.tools.plot import findNearest, GenericData

title='minimal'

#5% ethane, 95% Ar, variable temperature, P=1.4-4 atm
molFracList=[{'CC': 0.05, 'Ar': 0.95}]
Plist=[278643.8]
Tlist=range(1100,1300,100)
terminationTime = 5e-5
conditions=rt.generateAllConditions("IdealGasReactor", terminationTime, molFracList, Tlist, Plist)

# majorSpecies=['C', '[CH3]', 'C[CH2]', '[H]', 'C=C', '[H][H]', 'C=[CH]', 'CC[CH2]', 'C#C', 'C=CC=C'] #real one, one below is for testing, make sure to change order
majorSpecies=['CC', 'C=C', 'C#C','C', '[CH3]', 'C[CH2]', '[H]', '[H][H]', 'C=[CH]', 'CC[CH2]']

#Experimental Data
#I want to make this compatiable with csv files so you don't have to hardcode in the experimental data in python, csv are easier to manipulate anyhow
########################################################################################################################
shockTubeEffectiveHeatingTimes=collections.OrderedDict()
shockTubeEffectiveHeatingTimes[1100]=2.730
shockTubeEffectiveHeatingTimes[1200]=2.480
shockTubeEffectiveHeatingTimes[1300]=2.220
shockTubeEffectiveHeatingTimes[1400]=1.970
shockTubeEffectiveHeatingTimes[1500]=1.719
shockTubeEffectiveHeatingTimes[1600]=1.460
shockTubeEffectiveHeatingTimes[1700]=1.200
shockTubeEffectiveHeatingTimes[1800]=0.940
shockTubeEffectiveHeatingTimes[1900]=0.690



#Taken from figure 1 of Y. Hiadka, K. Sato, H. Hoshikawa, T. Nishimori, H. Tanaka, K. Inami, N. Ito. Combust. Flame, 120 3 (2000), pp. 245-264
def getExample1Data(conditionData, conditions, shockTubeEffectiveHeatingTimes, smilesInExperiments):
    """
    Parses the simulation's data to make generic data objects to compare with shock tube data taken from figure 1 of
    #Y. Hiadka, K. Sato, H. Hoshikawa, T. Nishimori, H. Tanaka, K. Inami, N. Ito. Combust. Flame, 120 3 (2000), pp. 245-264
    """
    #First get the time points close to
    allResults=[]
    for T, timepoint in shockTubeEffectiveHeatingTimes.iteritems():
        for smiles in smilesInExperiments:
            temperatureResults=GenericData(label=label='Temperature',
                                      data = []
                                      units = 'K')
            for conIndex, data in conditionData:
                if conditions[conIndex].T0==T:
                    #index1 is the index of time from the simulation closest to timepoint
                    #data is organized as (timeArray, rest of data)
                    index1=findNearest(data[0], timepoint)
                    for speciesData in data[1]:
                        if smiles==speciesData.species:

                    if smiles not in exptCompDict1: exptCompDict1[smiles]=[]
                    #need to multiply by 20, in the experiment, they don't seem to count the Ar in the mole fraction
                    exptCompDict1[smiles].append(resultsDictionary[T][1][index1][index2+1]*20)
                (index1, timepoint2)=rt.getNearestTime(timepoint1, resultsDictionary[T][2])
                for index2, smiles in enumerate(majorSpecies[0:4]):
                    if smiles not in exptCompDict2: exptCompDict2[smiles]=[]
                    exptCompDict2[smiles].append(resultsDictionary[T][3][index1][index2+1]*20)


############Move to CSV files############
#Taken from figure 1 of Y. Hiadka, K. Sato, H. Hoshikawa, T. Nishimori, H. Tanaka, K. Inami, N. Ito. Combust. Flame, 120 3 (2000), pp. 245-264
shockTubeEffectiveHeatingTimes=collections.OrderedDict()
shockTubeEffectiveHeatingTimes[1100]=2.730
shockTubeEffectiveHeatingTimes[1200]=2.480
# shockTubeEffectiveHeatingTimes[1300]=2.220
# shockTubeEffectiveHeatingTimes[1400]=1.970
# shockTubeEffectiveHeatingTimes[1500]=1.719
# shockTubeEffectiveHeatingTimes[1600]=1.460
# shockTubeEffectiveHeatingTimes[1700]=1.200
# shockTubeEffectiveHeatingTimes[1800]=0.940
# shockTubeEffectiveHeatingTimes[1900]=0.690

C2H6molfrac=collections.OrderedDict()
C2H6molfrac[1180]=0.97
C2H6molfrac[1150]=0.79
C2H6molfrac[1310]=0.47
C2H6molfrac[1480]=0.12
C2H6molfrac[1530]=0.06
C2H6molfrac[1700]=0.03

COmolfrac=collections.OrderedDict()
COmolfrac[1310]=0.01
COmolfrac[1480]=0.35
COmolfrac[1530]=0.41
COmolfrac[1700]=0.46

C2H4molfrac=collections.OrderedDict()
C2H4molfrac[1100]=0.08
C2H4molfrac[1150]=0.19
C2H4molfrac[1310]=0.47
C2H4molfrac[1480]=0.50
C2H4molfrac[1560]=0.42
C2H4molfrac[1700]=0.26

C2H2molfrac=collections.OrderedDict()
C2H2molfrac[1340]=-0.04
C2H2molfrac[1480]=-0.16
C2H2molfrac[1560]=-0.20
C2H2molfrac[1700]=-0.26

CH4molfrac=collections.OrderedDict()
CH4molfrac[1180]=0.0
CH4molfrac[1310]=0.04
CH4molfrac[1480]=0.16
CH4molfrac[1530]=0.20
CH4molfrac[1700]=0.25
########################################################################################################################
    # Plot the results if matplotlib is installed.
    # See http://ma tplotlib.org/ to get it.
if '--plot' in sys.argv[1:]:
    import matplotlib.pyplot as plt
    #first test
    # plt.semilogy(times1, data1[:,1], 'r', label="sim1")
    # plt.semilogy(times2, data2[:,1], label='sim2')
    # plt.legend()
    # plt.show()

    exptCompDict1={}
    exptCompDict2={}

    #Plot vs experiment
    for T, timepoint1 in shockTubeEffectiveHeatingTimes.iteritems():
        (index1,timepoint2)=rt.getNearestTime(timepoint1, resultsDictionary[T][0])
        for index2, smiles in enumerate(majorSpecies[0:4]):
            if smiles not in exptCompDict1: exptCompDict1[smiles]=[]
            #need to multiply by 20, in the experiment, they don't seem to count the Ar in the mole fraction
            exptCompDict1[smiles].append(resultsDictionary[T][1][index1][index2+1]*20)
        (index1, timepoint2)=rt.getNearestTime(timepoint1, resultsDictionary[T][2])
        for index2, smiles in enumerate(majorSpecies[0:4]):
            if smiles not in exptCompDict2: exptCompDict2[smiles]=[]
            exptCompDict2[smiles].append(resultsDictionary[T][3][index1][index2+1]*20)

    x=1

    plt.clf()
    plt.subplot(2, 2, 1)
    plt.plot(C2H6molfrac.keys(), C2H6molfrac.values(), 'o', label="Expt")
    plt.plot(shockTubeEffectiveHeatingTimes.keys(),exptCompDict1['CC'], 'r', label='Trusted')
    plt.plot(shockTubeEffectiveHeatingTimes.keys(),exptCompDict2['CC'], 'g', label='New')
    plt.xlabel('Temperature (K)')
    plt.ylabel('C2H6 Mole fraction')
    plt.subplot(2, 2, 2)
    plt.plot(C2H4molfrac.keys(), C2H4molfrac.values(), 'o', label="Expt")
    plt.plot(shockTubeEffectiveHeatingTimes.keys(),exptCompDict1['C=C'], 'r', label='Trusted')
    plt.plot(shockTubeEffectiveHeatingTimes.keys(),exptCompDict2['C=C'], 'g', label='New')
    plt.xlabel('Time (ms)')
    plt.ylabel('C2H4 Mole fraction')
    plt.subplot(2, 2, 3)
    plt.plot(C2H2molfrac.keys(), C2H2molfrac.values(), 'o', label="Expt")
    plt.plot(shockTubeEffectiveHeatingTimes.keys(),exptCompDict1['C#C'], 'r', label='Trusted')
    plt.plot(shockTubeEffectiveHeatingTimes.keys(),exptCompDict2['C#C'], 'g', label='New')
    plt.xlabel('Time (ms)')
    plt.ylabel('C2H2 Mole Fraction')
    plt.subplot(2, 2, 4)
    plt.plot(CH4molfrac.keys(), CH4molfrac.values(), 'o', label="Expt")
    plt.plot(shockTubeEffectiveHeatingTimes.keys(),exptCompDict1['C'], 'r', label='Trusted')
    plt.plot(shockTubeEffectiveHeatingTimes.keys(),exptCompDict2['C'], 'g', label='New')
    plt.xlabel('Time (ms)')
    plt.ylabel('CH4 Mole Fraction')
    plt.tight_layout()
    plt.show()
else:
    print("To view a plot of these results, run this script with the option --plot")