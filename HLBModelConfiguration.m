

% KEY INPUTS
name = 'test';  % Simulation name (used for saving output txt files and figures)
nRuns = 100;      % Number of repeated simulations
nYrs = 20;      % Length of simulation
spyr = 12;       % Number of data save points per year
tStepMax = 1e7; % Maximum number of time steps

%FIGURES  - SET =1 IF FIGURES ARE REQUIRED
mapFigure=0;                % Maps of the landscape
multiFigure=1;              % Figures showing averages from multiple simulations
individualFigures=1;        % Figures showing results from individual simulations
nIndividualSimFigs = 1;     % number of individual simulation figures required

% INITIAL CONDITIONS
fullInfestation = 1; % Set the vector to be established everywhere
initInfNum = 1; %Initially infected units of citrus
nInfected0 = 1; % Number of cells infected
nInfested0 = 1; % Number of cells infested

%initial infection - can select "random2 or choose a number (has to
%correspond to a cell number in the region)
infected0 = "random";
%initial infestion - can select "random" or choose a number (has to
%correspond to a cell number in the region). 
% other options are "current" and "galacia", which correspond to the
% current vector presence or a random selection of 10 from Galacia which
% was used for dispersal fitting. These both require simulation of the
% whole of Iberia
infested0 = "random";

% Detection parameters
nInspCom = [1 100]; % percentage of commercial cells to sample
nInspRsd = [0 0]; % percentage of residential cells to sample
nSamp = [5 100]; % number of samples per cell
inspInt = [1 1/2]; % Before and after the first detection
detProb = 0%0.5;
complianceProb = 0.9;
roguingProb = 0.9;
chemicalIncrease = 0; %Proportion decrease of remaining vector pop 
stopAtDetection = 0;

% GEOGRAPHIC LIMITS
region="valencia"
% Choose from: "valencia", "murcia", "seville", "malaga",
% "huelva", "bilbao", "algarve", "lisbon", "porto"
% Input 0 to use the whole of Iberia
% Alternatively input a matrix structured: [xmin xmax; ymin ymax]


% DISPERSAL PARAMETERS
delta = 1000;               % DISPERSAL RATE
zeta = 1/1000;              % PROPORTION OF DISPERSALS THAT ARE LONG DISTANCE
estRate=1/365;              % VECTOR ESTABLISHMENT RATE
alphaLoc = 1.96;            % LOCAL DISPERSAL DISTANCE PARAMETER
maxDist = 6;                % MAXIMUM DISTANCE OF LOCAL DISPERSAL
alphaLD = 130;              % LONG DISTANCE DISPERSAL SCALE
extV = 0;                   % EXTERNAL INTRODUCTION OF PSYLLIDs per year
m = 0.9;                    % POPULATION REDUCTION DUE TO COMMERCIAL PEST MANAGEMENT
p = 0.99;                   % POPULATION REDUCTION PARASITE PRESENCE

% INFECTION PARAMETERS
beta = 10;                  % INFECTION RATE
rho = 0.4;                  % PROPORTION OF INFECTION WITHIN CELL
lRate = 1/365;              % RATE OF TRANSITION FROM L->C
sRate = 5/365;              % RATE OF TRANSITION FROM C->I
extI = 0;                   % EXTERNAL INTRODUCTIONS OF INFECTED CITRUS per year
phi = 1;                    %EFFECT OF CITRUS TRENDS ON INTRODUCTION


% MAX NUMBER OF CITRUS UNITS
UCcom = 100;
UCrsd = UCcom*10;


% Run the model and generate the figures. Returns each event and the time
% that it occured (event and t) for each simulation and the vector state (A
% and N correspond to Y and Z) and the disease state (ECIR) of each cell at regular intervals.
[event,t,ArsdInt,AcomInt,NrsdInt,NcomInt,ErsdInt,EcomInt,CrsdInt,CcomInt,IrsdInt,IcomInt,RrsdInt]...
    = HLB_model(name,nRuns,nYrs,spyr,tStepMax,region,fullInfestation,infested0,infected0,nInfected0,...
    nInfested0,initInfNum,delta,zeta,estRate,alphaLoc,maxDist,alphaLD,extV,m,p,beta,rho,lRate,sRate,...
    extI,phi,UCrsd,UCcom,nInspCom,nInspRsd,nSamp,inspInt,detProb,complianceProb,roguingProb,...
    chemicalIncrease,stopAtDetection,mapFigure,multiFigure,individualFigures,nIndividualSimFigs);

%{
% event list:
First column:
1 - residential
2 - commercial

Second column:
1 - Vector local dispersal 
2 - Vector long distance dispersal with no infection
3 - Vector long distance dispersal with infection
4 - Pathogen long distance dispersal with no vector establishment
5 - Vector long distance dispersal with no infection into other citrus type
(rsd to com or com to rsd)
6 - Vector long distance dispersal with infection into other citrus type
7 - Pathogen long distance dispersal with no vector establishment into other citrus type
8 - New vector primary infestation
9 - Vector transition from arrived to established (Y->Z)
10 - Pathogen local dispersal
11 - Pathogen primary infection
12 - Pathogen E->C transition
13 - Pathogen C->I transition
14 - Symptomatic citrus (I) detected

Negative numbers indicate events that do not change things:
-1 - Long distance dispersal without establishment or vector 
-2 - Long distance dispersal to somewhere without citrus 
-3 - Failed primary infestation or infection because cell is fully infested/infected
%}