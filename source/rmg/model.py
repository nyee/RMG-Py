#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#	RMG - Reaction Mechanism Generator
#
#	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#	RMG Team (rmg_dev@mit.edu)
#
#	Permission is hereby granted, free of charge, to any person obtaining a
#	copy of this software and associated documentation files (the 'Software'),
#	to deal in the Software without restriction, including without limitation
#	the rights to use, copy, modify, merge, publish, distribute, sublicense,
#	and/or sell copies of the Software, and to permit persons to whom the
#	Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#	DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
Contains classes for working with the reaction model generated by RMG.
"""

import logging
import numpy
import os

import constants
import settings
import reaction
import species
import unirxn.network

################################################################################

class ReactionModel:
	"""
	Represent a generic reaction model. A reaction model consists of `species`,
	a list of species, and `reactions`, a list of reactions.
	"""

	def __init__(self, species=None, reactions=None):
		self.species = species or []
		self.reactions = reactions or []

################################################################################

class CoreEdgeReactionModel:
	"""
	Represent a reaction model constructed using a rate-based screening
	algorithm. The species and reactions in the model itself are called the
	*core*; the species and reactions identified as candidates for inclusion in
	the model are called the *edge*. The attributes are:

	=========================  ==============================================================
	Attribute                  Description
	=========================  ==============================================================
	`core`                     The species and reactions of the current model core
	`edge`                     The species and reactions of the current model edge
	`absoluteTolerance`        The absolute tolerance used in the ODE/DAE solver
	`relativeTolerance`        The relative tolerance used in the ODE/DAE solver
	`fluxToleranceKeepInEdge`  The relative species flux below which species are discarded from the edge
	`fluxToleranceMoveToCore`  The relative species flux above which species are moved from the edge to the core
	`fluxToleranceInterrupt`   The relative species flux above which the simulation will halt
	`maximumEdgeSpecies`       The maximum number of edge species allowed at any time
	`termination`              A list of termination targets (i.e :class:`TerminationTime` and :class:`TerminationConversion` objects)
	`unirxnNetworks`           A list of unimolecular reaction networks (:class:`unirxn.network.Network` objects)
	`networkCount`             A counter for the number of unirxn networks created
	=========================  ==============================================================


	"""

	def __init__(self, core=None, edge=None):
		if core is None:
			self.core = ReactionModel()
		else:
			self.core = core
		if edge is None:
			self.edge = ReactionModel()
		else:
			self.edge = edge
		# The default tolerances mimic the original RMG behavior; no edge
		# pruning takes place, and the simulation is interrupted as soon as
		# a species flux higher than the validity
		self.fluxToleranceKeepInEdge = 0.0
		self.fluxToleranceMoveToCore = 1.0
		self.fluxToleranceInterrupt = 1.0
		self.absoluteTolerance = 1.0e-8
		self.relativeTolerance = 1.0e-4
		self.maximumEdgeSpecies = 1000000
		self.termination = []
		self.unirxnNetworks = []
		self.networkCount = 0

	def initialize(self, coreSpecies):
		"""
		Initialize a reaction model with a list `coreSpecies` of species to
		start out with.
		"""

		logging.info('')

		# Add all species present in input file to model core
		for spec in coreSpecies:
			self.enlarge(spec)

	def enlarge(self, newObject):
		"""
		Enlarge a reaction model by processing `newObject`. If `newObject` is a
		:class:`rmg.species.Species` object, then the species is moved from
		the edge to the core and reactions generated for that species, reacting
		with itself and with all other species in the model core. If `newObject`
		is a :class:`rmg.unirxn.network.Network` object, then reactions are
		generated for the species in the network with the largest leak flux.
		"""

		rxnList = []

		if isinstance(newObject, species.Species):

			newSpecies = newObject
			# Find reactions involving the new species as unimolecular reactant
			# or product (e.g. A <---> products)
			rxnList.extend(reaction.kineticsDatabase.getReactions([newSpecies]))
			# Find reactions involving the new species as bimolecular reactants
			# or products with itself (e.g. A + A <---> products)
			rxnList.extend(reaction.kineticsDatabase.getReactions([newSpecies, newSpecies]))
			# Find reactions involving the new species as bimolecular reactants
			# or products with other core species (e.g. A + B <---> products)
			for coreSpecies in self.core.species:
				if coreSpecies.reactive:
					rxnList.extend(reaction.kineticsDatabase.getReactions([newSpecies, coreSpecies]))

			# Add new species
			self.addSpeciesToCore(newSpecies)

			logging.info('')
			logging.info('After model enlargement:')

		elif isinstance(newObject, unirxn.network.Network) and settings.unimolecularReactionNetworks:

			network = newObject
			# Determine the species with the maximum leak flux
			maxSpecies, maxSpeciesFlux = network.getMaximumLeakSpecies()
			network.explored.append(maxSpecies)
			# Find reactions involving the found species as unimolecular
			# reactant or product (e.g. A <---> products)
			rxnList = reaction.kineticsDatabase.getReactions([maxSpecies])
			# Don't find reactions involving the new species as bimolecular
			# reactants or products with itself (e.g. A + A <---> products)
			# Don't find reactions involving the new species as bimolecular
			# reactants or products with other core species (e.g. A + B <---> products)

			logging.info('')
			logging.info('After network enlargement:')

		else:
			raise TypeError('Unable to use object %s to enlarge reaction model; expecting an object of class rmg.species.Species or rmg.unirxn.network.Network.' % newObject)

		# Add new reactions generated in above
		for rxn in rxnList:
			allSpeciesInCore = True
			for spec in rxn.reactants:
				if spec not in self.core.species: allSpeciesInCore = False
				if spec not in self.edge.species and spec not in self.core.species:
					self.addSpeciesToEdge(spec)
			for spec in rxn.products:
				if spec not in self.core.species: allSpeciesInCore = False
				if spec not in self.edge.species and spec not in self.core.species:
					self.addSpeciesToEdge(spec)
			# We only add reactions that are not unimolecular if pressure
			# dependence is on; unimolecular reactions will be added after
			# processing the associated networks
			if not settings.unimolecularReactionNetworks or not (
				rxn.isIsomerization() or rxn.isDissociation() or rxn.isAssociation()):
				if allSpeciesInCore:
					self.addReactionToCore(rxn)
				else:
					self.addReactionToEdge(rxn)
			else:
				# Update unimolecular reaction networks
				self.addReactionToUnimolecularNetworks(rxn)

		# Output current model size information after enlargement
		logging.info('\tThe model core has %s species and %s reactions' % (len(self.core.species), len(self.core.reactions)))
		logging.info('\tThe model edge has %s species and %s reactions' % (len(self.edge.species), len(self.edge.reactions)))
		logging.info('')

	def addSpeciesToCore(self, spec):
		"""
		Add a species `spec` to the reaction model core (and remove from edge if
		necessary). This function also moves any reactions in the edge that gain
		core status as a result of this change in status to the core.
		"""

		if spec in self.core.species: return

		# Add the species to the core
		self.core.species.append(spec)

		# Add it to the cantera list
		spec.toCantera()

		if spec in self.edge.species:

			# If species was in edge, remove it
			self.edge.species.remove(spec)

			# Search edge for reactions that now contain only core species;
			# these belong in the model core and will be moved there
			rxnList = []
			for rxn in self.edge.reactions:
				allCore = True
				for reactant in rxn.reactants:
					if reactant not in self.core.species: allCore = False
				for product in rxn.products:
					if product not in self.core.species: allCore = False
				if allCore: rxnList.append(rxn)

			# Move any identified reactions to the core
			for rxn in rxnList:
				self.addReactionToCore(rxn)

	def addSpeciesToEdge(self, spec):
		"""
		Add a species `spec` to the reaction model edge.
		"""
		self.edge.species.append(spec)

	def removeSpeciesFromEdge(self, spec):
		"""
		Remove species `spec` from the reaction model edge.
		"""
		# remove the species
		self.edge.species.remove(spec)
		# identify any reactions it's involved in
		rxnList = []
		for rxn in self.edge.reactions:
			if spec in rxn.reactants or spec in rxn.products:
				rxnList.append(rxn)
		# remove those reactions
		for rxn in rxnList:
			self.edge.reactions.remove(rxn)

		# Remove the species from any unirxn networks it is in
		if settings.unimolecularReactionNetworks:
			networksToDelete = []
			for network in self.unirxnNetworks:
				if spec in network.getSpeciesList():
					# Delete all path reactions involving the species
					rxnList = []
					for rxn in network.pathReactions:
						if spec in rxn.reactants or spec in rxn.products:
							rxnList.append(rxn)
					for rxn in rxnList:
						network.pathReactions.remove(rxn)
					# Delete all net reactions involving the species
					rxnList = []
					for rxn in network.netReactions:
						if spec in rxn.reactants or spec in rxn.products:
							rxnList.append(rxn)
					for rxn in rxnList:
						network.netReactions.remove(rxn)
					# Delete all isomers involving the species
					isomerList = []
					for isomer in network.isomers:
						if spec in isomer.species:
							isomerList.append(isomer)
					for isomer in isomerList:
						network.isomers.remove(isomer)
					# If no remaining reactions, delete the network (actually
					# add to list of networks to be deleted in a subsequent
					# step)
					if len(network.pathReactions) == 0 and len(network.netReactions) == 0:
						networksToDelete.append(network)

			# Complete deletion of empty networks
			for network in networksToDelete:
				logging.debug('Deleting empty unirxn network %s' % network.id)
				self.unirxnNetworks.remove(network)



		# remove from the global list of species, to free memory
		species.speciesList.remove(spec)

	def addReactionToCore(self, rxn):
		"""
		Add a reaction `rxn` to the reaction model core (and remove from edge if
		necessary). This function assumes `rxn` has already been checked to
		ensure it is supposed to be a core reaction (i.e. all of its reactants
		AND all of its products are in the list of core species).
		"""
		self.core.reactions.append(rxn)
		if rxn in self.edge.reactions:
			self.edge.reactions.remove(rxn)

		# add it to the Cantera list
		rxn.toCantera()

	def addReactionToEdge(self, rxn):
		"""
		Add a reaction `rxn` to the reaction model edge. This function assumes
		`rxn` has already been checked to ensure it is supposed to be an edge
		reaction (i.e. all of its reactants OR all of its products are in the
		list of core species, and the others are in either the core or the
		edge).
		"""
		self.edge.reactions.append(rxn)

	def getLists(self):
		"""
		Return lists of all of the species and reactions in the core and the
		edge.
		"""
		speciesList = []
		speciesList.extend(self.core.species)
		speciesList.extend(self.edge.species)
		reactionList = []
		reactionList.extend(self.core.reactions)
		reactionList.extend(self.edge.reactions)
		return speciesList, reactionList

	def getStoichiometryMatrix(self):
		"""
		Return the stoichiometry matrix for all core and edge species. The
		rows represent the species in the core and edge in order, while the
		columns represent the reactions in the core and edge in order.
		"""
		speciesList, reactionList = self.getLists()
		stoichiometry = numpy.zeros((len(speciesList), len(reactionList)), float)
		for j, rxn in enumerate(reactionList):
			for i, spec in enumerate(speciesList):
				stoichiometry[i,j] = rxn.getStoichiometricCoefficient(spec)
		return stoichiometry

	def getReactionRates(self, T, P, Ci):
		"""
		Return an array of reaction rates for each reaction in the model core
		and edge. The core reactions occupy the first rows of the array, while
		the edge reactions occupy the last rows.
		"""
		speciesList, reactionList = self.getLists()
		rxnRate = numpy.zeros(len(reactionList), float)
		totalConc = sum( Ci.values() )
		for j, rxn in enumerate(reactionList):
			rxnRate[j] = rxn.getRate(T, P, Ci, totalConc)
		return rxnRate

	def addReactionToUnimolecularNetworks(self, newReaction):
		"""
		Given a newly-created :class:`Reaction` object `newReaction`, update the
		corresponding unimolecular reaction network. If no network exists, a new
		one is created. If the new reaction is an isomerization that connects two
		existing networks, the two networks are merged. This function is called
		whenever a new high-pressure limit edge reaction is created. Returns the
		network containing the new reaction.
		"""

		def getNetworkContainingSpecies(spec, networks):
			for network in networks:
				if network.containsSpecies(spec): return network
			return None

		# If the reaction is not unimolecular, then there's nothing to do
		# Return None because we haven't added newReaction to a network
		if not newReaction.isIsomerization() and not newReaction.isDissociation() and not newReaction.isAssociation():
			return None

		# Don't add reactions that are dissociation of a diatomic or association
		# to form a diatomic
		if newReaction.isDissociation() and len(newReaction.reactants[0].structure[0].atoms()) == 2:
			return None
		elif newReaction.isAssociation() and len(newReaction.products[0].structure[0].atoms()) == 2:
			return None

		# Find networks containing either the reactant or the product as a
		# unimolecular isomer
		reactantNetwork = None; productNetwork = None
		if newReaction.isIsomerization() or newReaction.isDissociation():
			reactantNetwork = getNetworkContainingSpecies(newReaction.reactants[0], self.unirxnNetworks)
		if newReaction.isIsomerization() or newReaction.isAssociation():
			productNetwork = getNetworkContainingSpecies(newReaction.products[0], self.unirxnNetworks)

		# Action depends on results of above search
		network = None
		if reactantNetwork is not None and productNetwork is not None and reactantNetwork is productNetwork:
			# Found the same network twice, so we don't need to merge
			# This can happend when a net reaction is later generated by RMG as a
			# path reaction
			network = reactantNetwork
		elif reactantNetwork is not None and productNetwork is not None:
			# Found two different networks, so we need to merge
			network = reactantNetwork
			network.merge(productNetwork)
			self.unirxnNetworks.remove(productNetwork)
		elif reactantNetwork is not None and productNetwork is None:
			# Found one network, so we add to it
			network = reactantNetwork
		elif reactantNetwork is None and productNetwork is not None:
			# Found one network, so we add to it
			network = productNetwork
		else:
			# Didn't find any networks, so we need to make a new one
			import unirxn.network
			self.networkCount += 1
			network = unirxn.network.Network(id=self.networkCount)
			self.unirxnNetworks.append(network)

		# Add the new reaction to whatever network was found/created above
		network.addPathReaction(newReaction)

		# Return the network that the reaction was added to
		return network

	def updateUnimolecularReactionNetworks(self):
		"""
		Iterate through all of the currently-existing unimolecular reaction
		networks, updating those that have been marked as invalid. In each update,
		the phenomonological rate coefficients :math:`k(T,P)` are computed for
		each net reaction in the network, and the resulting reactions added or
		updated.
		"""

		from unirxn.network import Isomer
		from reaction import PDepReaction
		from kinetics import ChebyshevKinetics, PDepArrheniusKinetics

		for network in self.unirxnNetworks:
			if not network.valid:

				logging.debug('Updating unimolecular reaction network %s' % network.id)

				# Other inputs
				method, Tlist, Plist, grainSize, numGrains, model = settings.unimolecularReactionNetworks

				network.bathGas = [spec for spec in self.core.species if not spec.reactive][0]
				network.bathGas.expDownParam = 4.86 * 4184

				# Generate isomers
				for reaction in network.pathReactions:

					# Create isomer for the reactant
					isomer = None
					for isom in network.isomers:
						if all([spec in isom.species for spec in reaction.reactants]):
							isomer = isom
					if isomer is None:
						isomer = Isomer(reaction.reactants)
						network.isomers.append(isomer)
					reaction.reactant = isomer

					# Create isomer for the product
					isomer = None
					for isom in network.isomers:
						if all([spec in isom.species for spec in reaction.products]):
							isomer = isom
					if isomer is None:
						isomer = Isomer(reaction.products)
						network.isomers.append(isomer)
					reaction.product = isomer

				# Update list of explored isomers to include all species in core
				for isom in network.isomers:
					if isom.isUnimolecular():
						spec = isom.species[0]
						if spec not in network.explored:
							if spec in self.core.species:
								network.explored.append(spec)

				# Remove any isomers that aren't found in any path reactions
				# Ideally this block of code wouldn't be needed, but it's here
				# just in case
				isomerList = []
				for isomer in network.isomers:
					found = False
					for reaction in network.pathReactions:
						if reaction.reactant is isomer or reaction.product is isomer:
							found = True
							break
					if not found:
						isomerList.append(isomer)
				if len(isomerList) > 0:
					logging.debug('Removed %i isomer(s) from network %i.' % (len(isomerList), network.id))
					for isomer in isomerList: network.isomers.remove(isomer)

				# Sort isomers so that all unimolecular isomers come first
				isomers = [isom for isom in network.isomers if isom.isUnimolecular()]
				isomers.extend([isom for isom in network.isomers if isom.isMultimolecular()])
				network.isomers = isomers

				# Get list of species in network
				speciesList = network.getSpeciesList()

				# Calculate ground-state energy of all species in network
				# For now we assume that this is equal to the enthalpy of formation
				# of the species
				for spec in speciesList:
					spec.E0 = spec.getEnthalpy(T=298)

				# Determine isomer ground-state energies
				for isomer in network.isomers:
					isomer.E0 = sum([spec.E0 for spec in isomer.species])
				# Determine transition state ground-state energies of the reactions
				for reaction in network.pathReactions:
					E0 = sum([spec.E0 for spec in reaction.reactants])
					reaction.E0 = E0 + reaction.kinetics[0].getActivationEnergy(reaction.getEnthalpyOfReaction(T=298))

				# Shift network such that lowest-energy isomer has a ground state of 0.0
				network.shiftToZeroEnergy()

				# Determine energy grains
				Elist = network.determineEnergyGrains(grainSize, numGrains, max(Tlist))

				# Calculate density of states for all isomers in network
				network.calculateDensitiesOfStates(Elist)

				# Determine phenomenological rate coefficients
				K = network.calculateRateCoefficients(Tlist, Plist, Elist, method)

				# Generate PDepReaction objects
				for i, product in enumerate(network.isomers):
					for j, reactant in enumerate(network.isomers[0:i]):
						# Find the path reaction
						netReaction = None
						for r in network.netReactions:
							if r.hasTemplate(reactant.species, product.species):
								netReaction = r
						# If path reaction does not already exist, make a new one
						if netReaction is None:
							netReaction = PDepReaction(reactant.species, product.species, network, None)
							network.netReactions.append(netReaction)
							self.addReactionToEdge(netReaction)
						# Set its kinetics using interpolation model
						if model[0].lower() == 'chebyshev':
							modelType, degreeT, degreeP = model
							chebyshev = ChebyshevKinetics()
							chebyshev.fitToData(Tlist, Plist, K[:,:,i,j], degreeT, degreeP)
							netReaction.kinetics = chebyshev
						elif model.lower() == 'pdeparrhenius':
							pDepArrhenius = PDepArrheniusKinetics()
							pDepArrhenius.fitToData(Tlist, Plist, K[:,:,i,j])
							netReaction.kinetics = pDepArrhenius
						else:
							pass

						# Update cantera if this is a core reaction
						if netReaction in self.core.reactions:
							netReaction.toCantera()

				for spec in speciesList:
					del spec.E0
				for reaction in network.pathReactions:
					del reaction.reactant
					del reaction.product
					del reaction.E0

				network.valid = True

	def loadSeedMechanism(self, path):
		"""
		Loads a seed mechanism from the folder indicated by `path` into the
		core-edge reaction model.
		"""

		import os.path
		import quantities as pq
		import data
		import species
		import kinetics
		import reaction
		
		# Load the species list from the file species.txt
		# This file has the format of a standard RMG dictionary
		d = data.Dictionary()
		d.load(os.path.join(path, 'species.txt'))
		d.toStructure(addH=True)
		speciesDict = {}
		for label, struct in d.iteritems():
			speciesDict[label] = species.makeNewSpecies(struct, label, reactive=True)
		print species.speciesList
		
		# Load the reactions from the file reaction.txt
		reactionList = []
		f = open(os.path.join(path, 'reactions.txt'), 'r')
		for line in f:
			line = data.removeCommentFromLine(line)
			line.strip()
			if len(line) > 0:
				items = line.split()
				if len(items) > 0:
					rxn = items[0:-6]

					# Extract reactants and products
					if '<=>' in rxn: arrow = rxn.index('<=>')
					elif '=>' in rxn: arrow = rxn.index('=>')
					else: raise IOError('No arrow found in reaction equation from line %s' % line)
					reactants = rxn[0:arrow:2]
					products = rxn[arrow+1::2]

					# Remove third body 'M' if present
					thirdBody = False
					if 'M' in reactants and 'M' in products:
						thirdBody = True
						reactants.remove('M')
						products.remove('M')
					
					# Convert strings to species objects
					reactants = [speciesDict[r] for r in reactants]
					products = [speciesDict[r] for r in products]

					# Process Arrhenius parameters
					order = len(reactants)
					if (thirdBody): order += 1
					Aunits = 'cm^%i/(mol^%i*s)' % (3*(order-1), order-1)
					A = float(pq.Quantity(float(items[-6]), Aunits).simplified)
					n = float(items[-5])			# dimensionless
					Ea = float(pq.Quantity(float(items[-4]), 'cal/mol').simplified)
					kin = [kinetics.ArrheniusKinetics(A=A, n=n, Ea=Ea)]

					# Create reaction object and add to list
					rxn = reaction.Reaction(reactants, products, 'seed', kin, thirdBody=thirdBody)
					reaction.processNewReaction(rxn)
					reactionList.append(rxn)

		f.close()
		
		# Add species to core
		for label, spec in speciesDict.iteritems():
			self.addSpeciesToCore(spec)
		# Add reactions to core
		for rxn in reactionList:
			self.addReactionToCore(rxn)

################################################################################

class TemperatureModel:
	"""
	Represent a temperature profile. Currently the only implemented model is
	isothermal (constant temperature).
	"""

	def __init__(self):
		self.type = ''
		self.temperatures = []

	def isIsothermal(self):
		return self.type == 'isothermal'

	def setIsothermal(self, temperature):
		self.type = 'isothermal'
		self.temperatures = [ [0.0, temperature] ]

	def getTemperature(self, time=None):
		if self.isIsothermal():
			return self.temperatures[0][1]
		else:
			return None

	def __str__(self):
		string = 'Temperature model: ' + self.type + ' '
		if self.isIsothermal():
			string += str(self.getTemperature(0))
		return string

################################################################################

class PressureModel:
	"""
	Represent a pressure profile. Currently the only implemented model is
	isobaric (constant pressure).
	"""

	def __init__(self):
		self.type = ''
		self.pressures = []

	def isIsobaric(self):
		return self.type == 'isobaric'

	def setIsobaric(self, pressure):
		self.type = 'isobaric'
		self.pressures = [ [0.0, pressure] ]

	def getPressure(self, time=None):
		if self.isIsobaric():
			return self.pressures[0][1]
		else:
			return None

	def __str__(self):
		string = 'Pressure model: ' + self.type + ' '
		if self.isIsobaric():
			string += str(self.getPressure(0))
		return string

################################################################################

class IdealGas:
	"""
	An equation of state based on the ideal gas approximation

	.. math::

		f(P, V, T, \\mathbf{N}) = NRT - PV

	where :math:`N \\equiv \\sum_i N_i` is the total number of moles.

	The ideal gas approximation is generally valid for gases at low pressures
	and moderate tempertaures; it does not predict the gas-liquid phase
	transition and is not applicable to liquids.
	"""

	def getTemperature(self, P, V, N):
		"""
		Return the temperature associated with pressure `P`, volume `V`, and
		numbers of moles `N`.
		"""
		return P * V / (sum(N) * constants.R)

	def getPressure(self, T, V, N):
		"""
		Return the pressure associated with temperature `T`, volume `V`, and
		numbers of moles `N`.
		"""
		return sum(N) * constants.R * T / V

	def getVolume(self, T, P, N):
		"""
		Return the volume associated with temperature `T`, pressure `P`, and
		numbers of moles `N` (which may be a list, in which case it's summed).
		"""
		try:
		    N = sum(N)
		except TypeError: #can't iterate; N probably a float not a list
		    pass
		return N * constants.R * T / P

	def getdPdV(self, P, V, T, N):
		"""
		Return the derivative :math:`\\frac{dP}{dV}\\big|_{T,\mathbf{N}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`.
		"""
		return (-P/V)

	def getdPdT(self, P, V, T, N):
		"""
		Return the derivative :math:`\\frac{dP}{dT}\\big|_{V,\mathbf{N}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`.
		"""
		return (P/T)

	def getdVdT(self, P, V, T, N):
		"""
		Return the derivative :math:`\\frac{dV}{dT}\\big|_{P,\mathbf{N}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`.
		"""
		return (V/T)

	def getdVdP(self, P, V, T, N):
		"""
		Return the derivative :math:`\\frac{dV}{dP}\\big|_{T,\mathbf{N}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`.
		"""
		return 1.0 / self.getdPdV(P, V, T, N)

	def getdTdP(self, P, V, T, N):
		"""
		Return the derivative :math:`\\frac{dT}{dP}\\big|_{V,\mathbf{N}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`.
		"""
		return 1.0 / self.getdPdT(P, V, T, N)

	def getdTdV(self, P, V, T, N):
		"""
		Return the derivative :math:`\\frac{dT}{dV}\\big|_{P,\mathbf{N}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`.
		"""
		return 1.0 / self.getdVdT(P, V, T, N)

	def getdPdNi(self, P, V, T, N, i):
		"""
		Return the derivative :math:`\\frac{dP}{dN_i}\\big|_{T, V,\mathbf{N}_{j \\ne i}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`. The final parameter `i` is used to determine which
		species to use; if `N` is a list, then `i` is an index, while if `N` is
		a dictionary, `i` is a key.
		"""
		if type(N) is dict: return (P/numpy.sum(N.values()))
		else: return (P/numpy.sum(N))

	def getdVdNi(self, P, V, T, N, i):
		"""
		Return the derivative :math:`\\frac{dV}{dN_i}\\big|_{T, P,\mathbf{N}_{j \\ne i}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`. The final parameter `i` is the index of the
		species of interest, corresponding to an index into the list `N`.
		"""
		if type(N) is dict: return (V/numpy.sum(N.values()))
		else: return (V/numpy.sum(N))

	def getdTdNi(self, P, V, T, N, i):
		"""
		Return the derivative :math:`\\frac{dT}{dN_i}\\big|_{P, V,\mathbf{N}_{j \\ne i}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`. The final parameter `i` is the index of the
		species of interest, corresponding to an index into the list `N`.
		"""
		if type(N) is dict: return (-T/numpy.sum(N.values()))
		else: return (-T/numpy.sum(N))

################################################################################

class IncompressibleLiquid:
	"""
	An equation of state for incompressible liquids

	.. math:: f(P, V, T, \\mathbf{N}) = ?

	where :math:`N \\equiv \\sum_i N_i` is the total number of moles.

	Initialise with keyword arguments::
		il = IncompressibleLiquid(T=298, P=1E5, V=, N=, Vmol=)
	"""
	def __init__(self, T=None, P=None, V=None, N=None, Vmol=None):
		self.T = T
		self.P = P
		self.V = V
		self.N = N
		self.Vmol = Vmol # Molar volume

	def getTemperature(self, P, V, N):
		"""
		Return the temperature associated with pressure `P`, volume `V`, and
		numbers of moles `N`.
		"""
		if self.T:
			return self.T
		else:
			raise Exception("I'm a liquid, and can't deduce T from P,V,N.")

	def getPressure(self, T, V, N):
		"""
		Return the pressure associated with temperature `T`, volume `V`, and
		numbers of moles `N`.
		"""
		if self.P:
			return self.P
		else:
			raise Exception("I'm a liquid, and can't deduce P from T,V,N.")

	def getVolume(self, T, P, N):
		"""
		Return the volume associated with temperature `T`, pressure `P`, and
		numbers of moles `N` (which may be a list, in which case it's summed).
		"""
		if self.V:
			logging.debug("Using explicit volume.")
			return self.V
		try:
		    N = sum(N)
		except TypeError: #can't iterate; N probably a float not a list
		    pass
		if self.Vmol:
			return N * self.Vmol
		else:
			raise Exception("I'm a liquid with unknown molarVolume, and can't deduce V from T,P,N.")

	def getdPdV(self, P, V, T, N):
		"""
		Return the derivative :math:`\\frac{dP}{dV}\\big|_{T,\mathbf{N}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`.
		"""
		return  0

	def getdPdT(self, P, V, T, N):
		"""
		Return the derivative :math:`\\frac{dP}{dT}\\big|_{V,\mathbf{N}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`.
		"""
		return 0

	def getdVdT(self, P, V, T, N):
		"""
		Return the derivative :math:`\\frac{dV}{dT}\\big|_{P,\mathbf{N}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`.
		"""
		return 0

	def getdVdP(self, P, V, T, N):
		"""
		Return the derivative :math:`\\frac{dV}{dP}\\big|_{T,\mathbf{N}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`.
		"""
		return 0

	def getdTdP(self, P, V, T, N):
		"""
		Return the derivative :math:`\\frac{dT}{dP}\\bigg|_{V,\mathbf{N}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`.
		"""
		return 0

	def getdTdV(self, P, V, T, N):
		"""
		Return the derivative :math:`\\frac{dT}{dV}\\big|_{P,\mathbf{N}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`.
		"""
		return 0

	def getdPdNi(self, P, V, T, N, i):
		"""
		Return the derivative :math:`\\frac{dP}{dN_i}\\big|_{T, V,\mathbf{N}_{j \\ne i}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`. The final parameter `i` is used to determine which
		species to use; if `N` is a list, then `i` is an index, while if `N` is
		a dictionary, `i` is a key.
		"""
		return numpy.inf
		### Warning: may be inconsistent with dVdNi and dVdP
		## if dVdNi>0 and dVdP=0 then is'nt dPdNi = infinity ?

	def getdTdNi(self, P, V, T, N, i):
		"""
		Return the derivative :math:`\\frac{dT}{dN_i}\\big|_{P, V,\mathbf{N}_{j \\ne i}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`. The final parameter `i` is the index of the
		species of interest, corresponding to an index into the list `N`.
		"""
		return 0

	def getdVdNi(self, P, V, T, N, i):
		"""
		Return the derivative :math:`\\frac{dV}{dN_i}\\big|_{T, P,\mathbf{N}_{j \\ne i}}`
		evaluated at a given pressure `P`, volume `V`, temperature `T`, and
		numbers of moles `N`. The final parameter `i` is the index of the
		species of interest, corresponding to an index into the list `N`.

		For lack of better information,
		we assume that the partial molar volume of species `i`
		is equal to the average molar volume of the mixture as a whole.
		"""

		if self.Vmol: # molar volume is set; assume it's the same for all components, and use it
			return self.Vmol

		# Else assume that the partial molar volume of species i
		# is equal to the average molar volume of all species
		if type(N) is dict: return (V/numpy.sum(N.values()))
		else: return (V/numpy.sum(N))



################################################################################

class InvalidReactionSystemException(Exception):
	"""
	An exception used when an invalid reaction system is encountered.
	"""

	def __init__(self, label):
		self.label = label

	def __str__(self):
		return 'Invalid reaction system: ' + self.label

################################################################################

class ReactionSystem:
	"""
	Represent a generic reaction system, e.g. a chemical reactor. A reaction
	system is defined by a temperature model `temperatureModel`, a pressure
	model `pressureModel`, a volume model `volumeModel`, and a dictionary of
	initial and constant concentrations `initialConcentration`. Only two of
	`temperatureModel`, `pressureModel`, and `volumeModel` are independent; the
	remaining one must be set to :data:`None`.

	Each RMG job can handle multiple reaction systems; the resulting model
	will generally be the union of the models that would have been generated
	via individual RMG jobs, and will therefore be valid across all reaction
	systems provided.
	"""

	def __init__(self, temperatureModel=None, pressureModel=None,
				 volumeModel=None, initialConcentration=None):
		self.setModels(temperatureModel, pressureModel, volumeModel)
		self.initialConcentration = initialConcentration or {}

	def setModels(self, temperatureModel, pressureModel, volumeModel):
		"""
		Set any two of the reactor's `temperatureModel`, `pressureModel` or
		`volumeModel`. Attempting to set all three to non-None will cause an
		:class:`InvalidReactorModelException` to be raised.
		"""
		# Fail if all three models are specified
		if temperatureModel is not None and pressureModel is not None and volumeModel is not None:
			raise InvalidReactionSystemException('Attempted to specify temperature, pressure, and volume models; can only specify two of these at a time.')
		# Otherwise set models
		self.temperatureModel = temperatureModel
		self.pressureModel = pressureModel
		self.volumeModel = volumeModel

################################################################################

class TerminationTime:
	"""
	Represent a time at which the simulation should be terminated. This class
	has one attribute: the termination `time` in seconds.
	"""

	def __init__(self, time=0.0):
		self.time = time

################################################################################

class TerminationConversion:
	"""
	Represent a conversion at which the simulation should be terminated. This
	class has two attributes: the `species` to monitor and the fractional
	`conversion` at which to terminate.
	"""

	def __init__(self, spec=None, conv=0.0):
		self.species = spec
		self.conversion = conv

################################################################################

if __name__ == '__main__':

	import chem
	import data
	import species
	import reaction
	import thermo

	import os.path
	import main
	main.initializeLog(logging.DEBUG)

	datapath = '../data/RMG_database/'

	logging.debug('General database: ' + os.path.abspath(datapath))
	species.thermoDatabase = species.ThermoDatabaseSet()
	species.thermoDatabase.load(datapath)
	thermo.forbiddenStructures = data.Dictionary()
	thermo.forbiddenStructures.load(datapath + 'forbiddenStructure.txt')
	thermo.forbiddenStructures.toStructure()
	#reaction.kineticsDatabase = reaction.ReactionFamilySet()
	#reaction.kineticsDatabase.load(datapath)

	structure = chem.Structure(); structure.fromSMILES('C')
	CH4 = species.makeNewSpecies(structure)

	structure = chem.Structure(); structure.fromSMILES('[H]')
	H = species.makeNewSpecies(structure)

	structure = chem.Structure(); structure.fromSMILES('[CH3]')
	CH3 = species.makeNewSpecies(structure)

	forward = reaction.Reaction([CH3, H], [CH4])
	reverse = reaction.Reaction([CH4], [CH3, H])
	forward.reverse = reverse
	reverse.reverse = forward

	kinetics = reaction.ArrheniusEPKinetics()
	kinetics.fromDatabase([300, 2000, 1.93E14, 0, 0, 0.27, 0, 0, 0, 0, 3], '', 2)
	forward.kinetics = [kinetics]

	speciesList = [CH3, H, CH4]
	reactionList = [forward]

	reactionSystem = BatchReactor()
	reactionSystem.temperatureModel = TemperatureModel()
	reactionSystem.temperatureModel.setIsothermal(pq.Quantity(1000, 'K'))
	reactionSystem.pressureModel = PressureModel()
	reactionSystem.pressureModel.setIsobaric(pq.Quantity(1, 'bar'))
	reactionSystem.equationOfState = IdealGas()
	reactionSystem.initialConcentration[CH4] = pq.Quantity(1, 'mol/m**3')

	reactionSystem.solve(0.0, 1.0e0, speciesList, reactionList)

	