import regressiveTest as rt

#Test Condition class
#x=rt.Condition('IdealGasReactor', 5, {'Ar':0.95, 'X':0.05}, 900.0, 100000.0, 0.09)
#print x.__repr__()
#print x.__str__()

#rt.Condition(reactorType="IdealGasReactor", reactionTime=5, molFrac={'X': 0.05, 'Ar': 0.95}, T0=900.0000000000, P0=100000.0000000000, V0=0.0900000000)

#Test generateAllConditions
# molFracList=[{'X': 0.05, 'Ar': 0.95}, {'X': 0.025, 'Ar': 0.975}]
# Tlist=[700,800,1000]
# Plist=[10000, 20000, 30000]
# Vlist=[0.1, 0.2, 0.5]
# conditions=rt.generateAllConditions("IdealGasReactor", 5, molFracList, Tlist, Plist, Vlist)
# conditions=rt.generateAllConditions("IdealGasReactor", 5, molFracList, Tlist, Plist)
# for condition in conditions:
#     print condition
molFracList=[{'CC': 0.05, 'Ar': 0.95}]
Plist=[278643.8]
Tlist=range(1100,1300,100)
terminationTime = 5e-5
conditions=rt.generateAllConditions("IdealGasReactor", terminationTime, molFracList, Tlist, Plist)

#test ObservablesTestCase
majorSpeciesSmiles=['CC', 'C=C', 'C#C','C', '[CH3]', 'C[CH2]', '[H]', '[H][H]', 'C=[CH]', 'CC[CH2]', 'Ar']

A=rt.ObservablesTestCase("Ethane Pyrolysis (Minimal)", "regression/egA",
                    "regression/egB", conditions, majorSpeciesSmiles, exptData=None)

print A.__str__()
A.compare()
