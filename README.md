# HLBinEurope

Use HLBModelConfiguration.m to input parameters and run the model.

This will call the matlab function files:
HLB_model.m - called from HLBModelConfiguration.m, uses the input parameters to run the model multiple times. It calls teh remaining functions:
localArea.m - calculates the cells within the 'local area' (within the maximum radius of the local dispersal kernel), for each cell and calculates the local dispersal kernel.
nextEvent.m - generates the which event to occur and where, given the rates of all events, using the Gillespie algorithm.
updateRates.m - updates all event rates after the last event occured.

There is an option to choose to generate figures in HLBModelConfiguration.m. If selected, the figures will be saved in a 'Figures' directory.
These are created by the following functions:
plotMaps.m - plots maps of the landscape used for the simulations.
plotMultiTotals.m - plots averages from multiple simulations.
plotSingleInf.m - plots infection results from single simulations.
plotSingleVec.m - plots vector dispersal results from single simulations.
plotTimeSeries.m - plots a storyboard of the infection of a single simulation at 5 regular intervals.

Additionally, text files will be saved in a 'ModelOutputs' directory:
Inputs - The model parameters used for the simulation.
Outputs - The number of cells/units in each compartment of the model at regular intervals.
CellInfections - The status of each cell at the end of the simulation.
InfectionTimes - The time that each infection event occurred.
