import sys
import math
import numpy as np
import regressiveTest as rt
import cantera as ct
import collections

#5% ethane, 95% Ar, variable temperature, P=1.4-4 atm
conditions=[]


# majorSpecies=['C', '[CH3]', 'C[CH2]', '[H]', 'C=C', '[H][H]', 'C=[CH]', 'CC[CH2]', 'C#C', 'C=CC=C'] #real one, one below is for testing, make sure to change order
majorSpecies=['CC', 'C=C', 'C#C','C', '[CH3]', 'C[CH2]', '[H]', '[H][H]', 'C=[CH]', 'CC[CH2]']
speciesList1=rt.getNameFromRMGDict(majorSpecies, '/Users/Nate/Dropbox (MIT)/Research/RMG/regressiveTests/minimal/species_dictionary.txt')

# speciesList2=rt.getSpeciesFromChemkin(majorSpecies, '/Users/Nate/Dropbox (MIT)/Research/RMG/regressiveTests/minimal/species_dictionary.txt')
speciesList2=rt.getNameFromRMGDict(majorSpecies, '/Users/Nate/Dropbox (MIT)/Research/RMG/regressiveTests/butane1/species_dictionary.txt')

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


#key is temperature of run and value is four lists, t1, y1, t2, y2
resultsDictionary={}

for T in shockTubeEffectiveHeatingTimes:
    #Run simulation
    minSim1 = ct.Solution('minimal.cti')
    # minSim2 = ct.Solution('minimal_mod.cti')
    minSim2 = ct.Solution('butane1.cti')
    # minSim2 = ct.Solution('minimal.cti')
    #Variable temperature, 2.75 atm, 5% ethane
    minSim1.TPX = T, 278643.8, 'Ar:0.95, {0}:0.05'.format(speciesList1[0])
    minSim2.TPX = T, 278643.8, 'Ar:0.95, {0}:0.05'.format(speciesList2[0])
    #Constant volume reactor
    r1 = ct.IdealGasReactor(minSim1)
    r2 = ct.IdealGasReactor(minSim2)

    sim1 = ct.ReactorNet([r1])
    sim2 = ct.ReactorNet([r2])
    time = 0.0
    times1 = np.zeros(100)
    times2 = np.zeros(100)
    data1 = np.zeros((100,len(majorSpecies)+1))
    data2 = np.zeros((100,len(majorSpecies)+1))

    print('%10s %10s %10s %14s' % ('t [s]','T [K]','P [Pa]','u [J/kg]'))
    for n in range(100):
        time += 5.e-5
        sim1.advance(time)
        times1[n] = time * 1e3  # time in ms
        data1[n,0] = r1.T
        data1[n,1:] = r1.thermo[speciesList1].X
        print('%10.3e %10.3f %10.3f %14.6e' % (sim1.time, r1.T,
                                               r1.thermo.P, r1.thermo.u))
    time = 0.0
    #Not sure if its a good idea to advance simulation in same for loop, so going to do a separate one
    print('%10s %10s %10s %14s' % ('t [s]','T [K]','P [Pa]','u [J/kg]'))
    for n in range(100):
        time += 5.e-5
        sim2.advance(time)
        times2[n] = time * 1e3  # time in ms
        data2[n,0] = r2.T
        data2[n,1:] = r2.thermo[speciesList2].X
        print('%10.3e %10.3f %10.3f %14.6e' % (sim2.time, r2.T,
                                               r2.thermo.P, r2.thermo.u))

    resultsDictionary[T]=[times1, data1, times2, data2]


    for index, species in enumerate(majorSpecies):
        if not rt.curvesSimilar(times1, [math.log(x) for x in data1[:,index+1]], times2, [math.log(x) for x in data2[:,index+1]], 0.05):
            print "{0} concentration is not very similar at temperature {1}".format(species, T)
            if not type(resultsDictionary[T][-1]) is bool: resultsDictionary[T].append(True)
    else:
        if resultsDictionary[T][-1]!=True: resultsDictionary[T].append(False)

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