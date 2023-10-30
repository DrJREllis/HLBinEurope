

% KEY INPUTS
name = 'test';  % Simulation name (used for saving output txt files and figures)
nRuns = 5;      % Number of repeated simulations
nYrs = 10;      % Length of simulation
spyr = 1;       % Number of data save points per year
tStepMax = 1e7; % Maximum number of time steps

%FIGURES  - SET =1 IF FIGURES ARE REQUIRED
mapFigure=0;                % Maps of the landscape
multiFigure=0;              % Figures showing averages from multiple simulations
individualFigures=0;        % Figures showing results from individual simulations
nIndividualSimFigs = 5;     % number of individual simulation figures required

% INITIAL CONDITIONS
fullInfestation = 1;
initInfNum = 20; %Initially infected units of citrus
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
nInspCom = [5 100]; % percentage of commercial cells to sample
nInspRsd = [0 0]; % percentage of residential cells to sample
nSamp = [0 100]; % number of samples per cell
inspInt = [1 1/2]; % Before and after the first detection
detProb = 0.5;
complianceProb = 0.95;
stopAtDetection = 0;

% GEOGRAPHIC LIMITS
region="valencia"
% Choose from: "valencia", "murcia", "seville", "malaga",
% "huelva", "bilbao", "algarve", "lisbon", "porto"
% Input 0 to use the whole of Iberia
% Alternatively input a matrix structured: [xmin xmax; ymin ymax]


% DISPERSAL PARAMETERS
delta = 1600;               % DISPERSAL RATE
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
rho = 0.7;                  % PROPORTION OF INFECTION WITHIN CELL
lRate = 1/365;              % RATE OF TRANSITION FROM L->C
sRate = 5/365;              % RATE OF TRANSITION FROM C->I
extI = 0;                   % EXTERNAL INTRODUCTIONs OF INFECTED CITRUS per year
phi = 1;                    %EFFECT OF CITRUS TRENDS ON INTRODUCTION


% MAX NUMBER OF CITRUS UNITS
UCcom = 100;
UCrsd = UCcom*10;



HLB_model(name,nRuns,nYrs,spyr,tStepMax,region,fullInfestation,infested0,infected0,delta,zeta,estRate,alphaLoc,maxDist,alphaLD,extV,m,p,...
    beta,rho,lRate,sRate,extI,phi,UCrsd,UCcom,nInspCom,nInspRsd,nSamp,inspInt,detProb,complianceProb,stopAtDetection,mapFigure,multiFigure,individualFigures,nIndividualSimFigs )

