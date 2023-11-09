# Framework for modelling HLB in Europe
  
Use HLBModelConfiguration.m to input parameters and run the model.  
  
This will call the matlab function files:  
HLB_model.m - called from HLBModelConfiguration.m, uses the input parameters to run the model multiple times. It calls the remaining functions:  
localArea.m - calculates the cells within the 'local area' (within the maximum radius of the local dispersal kernel), for each cell and calculates the local dispersal kernel.  
nextEvent.m - generates the which event to occur and where, given the rates of all events, using the Gillespie algorithm.  
updateRates.m - updates all event rates after the last event occured.  

In the data folder is two files:
citrusMatrices.mat - contains matrices of size 1136x1001 which indicate how much of each type of citrus is in each cell (total citrus, residential, commercial, abandoned, organic), the trends in citrus production and whether the vector is currently present.
climateMatrix.mat - contains two matrices of size 1136x1001 which give the climate suitability for AfCP and ACP (note that we only use AfCP in the model so far and the calculation for ACP is a rough estimate).
  
There is an option to choose to generate figures in HLBModelConfiguration.m. If selected, the figures will be saved in a 'Figures' directory.
These are created by the following functions:  
plotMaps.m - plots maps of the landscape used for the simulations.  
plotMultiTotals.m - plots averages from multiple simulations.  
plotStoryboard.m - plots a storyboard of the infection of a single simulation at 5 regular intervals.  
plotVectorStoryboard.m - plots a storyboard of the vector infestation of a single simulation at 5 regular intervals.  
plotMultiStoryboard.m - plots a storyboard of the infection from averages of multiple simulations at 5 regular intervals.  
  
Additionally, text files will be saved in a 'ModelOutputs' directory (see README in the directory for more info):  
Inputs - The model parameters used for the simulation.  
Outputs - The number of cells/units in each compartment of the model at regular intervals.  
CellInfections - The status of each cell at the end of the simulation.  
InfectionTimes - The time that each infection event occurred.  
