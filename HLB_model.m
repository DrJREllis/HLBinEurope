function [event,t,ArsdInt,AcomInt,NrsdInt,NcomInt,ErsdInt,EcomInt,CrsdInt,CcomInt,...
    IrsdInt,IcomInt,RrsdInt] = HLB_model(name,nRuns,nYrs,spyr,tStepMax,region,fullInfestation,...
    infested0,infected0,nInfected0,nInfested0,initInfNum,delta,zeta,estRate,...
    alphaLoc,maxDist,alphaLD,extV,m0,p,beta,rho,lRate,sRate,extI,phi,UCrsd,...
    UCcom,percInspCom,percInspRsd,nSamp,inspInt,detProb,complianceProb,roguingProb,...
    chemicalIncrease,stopAtDetection,mapFigure,multiFigure,individualFigures,...
    nIndividualSimFigs )


load("Data/citrusMatrices.mat");
load("Data/climateMatrix.mat");
climMat = climMatAfCP;

if isstring(region) || min(size(region)==[2,2])==1
    name=append(name,'_',region);
    if region=="valencia"
        limits = [750 800; 350 400];
    else if region== "murcia"
            limits= [700 750; 230 280];
    else if region=="seville"
            limits= [260 310; 250 300]; % could extend to [260 310; 250 370]
    else if region=="malaga"
            limits= [345 395; 150 200];
    else if region=="huelva" 
            limits= [145 195; 270 320];
    else if region=="bilbao"
            limits= [620 670; 840 890];
    else if region=="algarve"
            limits = [40 90; 280 330];
    else if region=="lisbon" 
            limits= [1 51; 470 520];
    else if region=="porto"
            limits= [120 170; 710 760];
    end
    end
    end
    end
    end
    end
    end
    end
    end

    x=limits(1,:); y=limits(2,:);
    citMat = citMat(x(1):x(2),y(1):y(2));
    comMat = comMat(x(1):x(2),y(1):y(2));
    rsdMat = rsdMat(x(1):x(2),y(1):y(2));
    abnMat = abnMat(x(1):x(2),y(1):y(2));
    orgMat = orgMat(x(1):x(2),y(1):y(2));
    climMat = climMat(x(1):x(2),y(1):y(2));
    vecMat = vecMat(x(1):x(2),y(1):y(2));

end

gridX=size(citMat,1); gridY=size(citMat,2);
ocean = nan(gridX,gridY); ocean(~isnan(citMat)) = 1;

paraPresence = zeros(gridX,gridY);

tMax = 365*nYrs;


citMat(isnan(citMat)) = 0;

nCells = sum(citMat>0,"all");
citNums = nan(size(citMat)); citNums((citMat>0)) = (1:nCells);
climArray=climMat(citMat>0); 
realComCit=comMat(citMat>0);
realRsdCit=rsdMat(citMat>0);
abnArray=abnMat(citMat>0);
orgArray=orgMat(citMat>0);
comArray = ceil(realComCit*UCcom);
rsdArray = ceil(realRsdCit*UCrsd);
% trendArray = (trendMat(citMat>0));
trendArray = (trendMat(citMat>0)*100+2);
trendArray = trendArray/max(trendArray);

localDispRate = delta*(1-zeta) ;               % LOCAL DISPERSAL
longDispRate = delta*zeta;                % LONG DISTANCE DISPERSAL

% Calculate dispersal kernel
[Fi, Kd, Ki, coord] = localArea(alphaLoc,rho,maxDist,gridX,gridY,citNums);



NrsdInt=nan(nCells,nYrs*spyr+1,nRuns); IrsdInt=nan(nCells,nYrs*spyr+1,nRuns); 
NcomInt=nan(nCells,nYrs*spyr+1,nRuns); IcomInt=nan(nCells,nYrs*spyr+1,nRuns); 
ArsdInt=nan(nCells,nYrs*spyr+1,nRuns); AcomInt=nan(nCells,nYrs*spyr+1,nRuns); 
ErsdInt=nan(nCells,nYrs*spyr+1,nRuns); EcomInt=nan(nCells,nYrs*spyr+1,nRuns); 
CrsdInt=nan(nCells,nYrs*spyr+1,nRuns); CcomInt=nan(nCells,nYrs*spyr+1,nRuns); 
RrsdInt=nan(nCells,nYrs*spyr+1,nRuns); RcomInt=nan(nCells,nYrs*spyr+1,nRuns); 
firstDetection=zeros(nRuns,1);
detectionCitrusIncidence=zeros(nRuns,1);
detectionCellIncidence=zeros(nRuns,1);


residentialOutputs=zeros(nRuns,6); commercialOutputs=zeros(nRuns,6); eventCount=zeros(16,2,nRuns);
rsdVecTime=nan(nCells,nRuns); rsdInfTime=nan(nCells,nRuns);
comVecTime=nan(nCells,nRuns); comInfTime=nan(nCells,nRuns);
stops=zeros(nRuns,4);

externalIRate = (extI/365)* realComCit.*( trendArray.^phi) /sum(realComCit.*( trendArray.^phi));      % EXTERNAL INTRODUCTION OF INFECTED CITRUS
externalVRate = (extV/365)* realComCit.*( trendArray.^phi) /sum(realComCit.*( trendArray.^phi));  

paraArray = paraPresence(citMat>0)*(1-p); paraArray(paraArray==0)=1;
comCC =  ceil(climArray.*realComCit.*paraArray );
rsdCC =  ceil(climArray.*realRsdCit.*paraArray );

abnArray(isnan(abnArray)) = 0; orgArray(isnan(orgArray)) = 0;
M = 1-(m0*(1-(abnArray+orgArray)./realComCit) + (1-paraArray).*(abnArray+orgArray)./realComCit);

extInfCSum = cumsum(reshape(externalIRate,nCells,1));
extVCSum = cumsum(reshape(externalVRate,nCells,1));


nInspCom = ceil(percInspCom*sum(comArray>0)/100);
nInspRsd = ceil(percInspRsd*sum(rsdArray>0)/100);

if not(isfolder("ModelOutputs"))
    mkdir("ModelOutputs")
end
InputLabel = ["nRuns"; "nYrs"; "spyr"; "tMax"; "tStepMax"; "No. of cells"; "delta"; "zeta"; "alphaLoc"; "maxDist"; "alphaLD"; "estRate"; "beta"; "rho"; "lRate"; "sRate"; "extV"; "extI"; "p"; "m"; "phi"];
Inputs = [nRuns; nYrs; spyr; tMax; tStepMax; nCells; delta; zeta; alphaLoc; maxDist; alphaLD; estRate; beta; rho; lRate; sRate; extV; extI; p; m0; phi];
InputsTab = table(nRuns, nYrs, spyr, tMax, tStepMax, nCells, delta, zeta, alphaLoc, maxDist, alphaLD, estRate, beta, rho, lRate, sRate, extV, extI, p, m0, phi);
writetable(InputsTab,append('ModelOutputs/',name,'_Inputs.txt'),'Delimiter','tab')  

Cell = (1:nCells)'; x=coord(:,1); y=coord(:,2);

infected0Indicator = 'NA';
if isstring(infected0)
infected0Indicator = infected0;
end
if isstring(infested0)
infested0Indicator = infested0;
end

for k=1:nRuns
tic
m=m0;
M = 1-(m0*(1-(abnArray+orgArray)./realComCit) + (1-paraArray).*(abnArray+orgArray)./realComCit);

adjNrsd = realRsdCit.*climArray.*paraArray./rsdCC; adjNrsd=max(0,adjNrsd); 
adjNcom = realComCit.*climArray.*M./comCC; adjNcom=max(0,adjNcom); 

if infected0Indicator == "random"
 infected0= randsample(Cell(comArray>0),nInfected0,true,comArray(comArray>0)); 
end

if fullInfestation~=1
if infested0Indicator == "random"
    infested0 = randsample(Cell(comArray>0),nInfested0);

else if infested0Indicator == "current"
    infested0=find(vecMat(citMat>0)==1); 

else if infested0Indicator == "galacia"
    infested0 = infested0(coord(infested0,2)>835 & coord(infested0,1)<250); 
    infested0(randperm(size(infested0,1),size(infested0,1)-10)) = []; 
end
end
end
end


Ecom=zeros(nCells,1); Ersd=zeros(nCells,1); 
Ccom=zeros(nCells,1); Crsd=zeros(nCells,1); 
Icom=zeros(nCells,1); Irsd=zeros(nCells,1); 
Rcom=zeros(nCells,1); Rrsd=zeros(nCells,1); 
Nrsd=zeros(nCells,1); Ncom=zeros(nCells,1); 
Arsd=zeros(nCells,1); Acom=zeros(nCells,1); 
forceVectorRsd = zeros(nCells,1); forceVectorCom = zeros(nCells,1);
forceInfectRsd = zeros(nCells,1); forceInfectCom = zeros(nCells,1);
longDistDispRsd = zeros(nCells,1); longDistDispCom = zeros(nCells,1);
estabRsd = zeros(nCells,1); estabCom = zeros(nCells,1);
latentRateRsd = zeros(nCells,1); latentRateCom = zeros(nCells,1);
symptomRateRsd = zeros(nCells,1); symptomRateCom = zeros(nCells,1);
lambdaVectorSumRsd = 0; lambdaVectorSumCom = extVCSum(end);
lambdaInfectionSumRsd = 0; lambdaInfectionSumCom = extInfCSum(end);


if fullInfestation==1
    Nrsd(rsdArray>0)=1;
    Ncom(comArray>0)=1;
else
for i=1:size(infested0,1)
     cell = infested0(i);
    Acom(cell) = 1;
    comVecTime(cell,k) = 1;
    [forceInfectCom,forceVectorCom,longDistDispCom,estabCom,...
        latentRateCom,symptomRateCom,forceVectorRsd,forceInfectRsd,...
         lambdaInfectionSumCom,lambdaVectorSumCom,lambdaInfectionSumRsd,lambdaVectorSumRsd,~,~,adjNcom] = ...
        updateRates(Ncom,Ecom,Ccom,Icom,Rcom,0,fullInfestation,forceInfectCom,forceVectorCom,1,M,...
        forceVectorRsd,forceInfectRsd,lambdaInfectionSumCom,lambdaVectorSumCom,lambdaInfectionSumRsd,lambdaVectorSumRsd,...
        longDistDispCom,estabCom,longDistDispRsd,extInfCSum,extVCSum,0,0,estabRsd,...
        latentRateCom,symptomRateCom,latentRateRsd,symptomRateRsd,beta,localDispRate,longDispRate,estRate,lRate,sRate,...
        comArray,adjNcom,climArray,comCC,...
        rsdCC,Nrsd,Ersd,Crsd,Irsd,Rrsd,rsdArray,comVecTime(:,1),comInfTime(:,1),Fi,Kd,Ki,cell,1,0);

end
end

if infected0~=0
    for i=1:size(infected0,1)
         cell = infected0(i);

         Ecom(cell) = min(initInfNum,comArray(cell));
         comInfTime(cell,k) = 0;         
        [forceInfectCom,forceVectorCom,longDistDispCom,estabCom,...
            latentRateCom,symptomRateCom,forceVectorRsd,forceInfectRsd,...
             lambdaInfectionSumCom,lambdaVectorSumCom,lambdaInfectionSumRsd,lambdaVectorSumRsd,~,~,adjNcom] = ...
            updateRates(Ncom,Ecom,Ccom,Icom,Rcom,0,fullInfestation,forceInfectCom,forceVectorCom,1,M,...
            forceVectorRsd,forceInfectRsd,lambdaInfectionSumCom,lambdaVectorSumCom,lambdaInfectionSumRsd,lambdaVectorSumRsd,...
            longDistDispCom,estabCom,longDistDispRsd,extInfCSum,extVCSum,0,0,estabRsd,...
            latentRateCom,symptomRateCom,latentRateRsd,symptomRateRsd,beta,localDispRate,longDispRate,estRate,lRate,sRate,...
            comArray,adjNcom,climArray,comCC,...
            rsdCC,Nrsd,Ersd,Crsd,Irsd,Rrsd,rsdArray,comVecTime(:,1),comInfTime(:,1),Fi,Kd,Ki,cell,11,0);

    end
end


t{k} = zeros(1,tStepMax); event{k}=nan(tStepMax,2);


ArsdInt(:,1,k) = Arsd; AcomInt(:,1,k) = Acom;
NrsdInt(:,1,k) = Nrsd; NcomInt(:,1,k) = Ncom;
ErsdInt(:,1,k) = Ersd; EcomInt(:,1,k) = Ecom;
CrsdInt(:,1,k) = Crsd; CcomInt(:,1,k) = Ccom;
IrsdInt(:,1,k) = Irsd; IcomInt(:,1,k) = Icom;
RrsdInt(:,1,k) = Rrsd; RcomInt(:,1,k) = Rcom;

numInfested{k}=nan(tMax,1); citInfected{k}=nan(tMax,1);  cellInfected{k}=nan(tMax,1);  
if fullInfestation~=1
numInfested{k}(1) = size(infested0,1);
else
    numInfested{k}(1) = sum(rsdArray>0)+sum(comArray>0);
end
citInfected{k}(1) = size(infected0,1);
cellInfected{k}(1) = size(infected0,1);

tStep = 1;

cellnum=zeros(10,1); xpos=zeros(10,1); ypos=zeros(10,1);

firstInspection = rand*365*inspInt(1);
complianceArray = binornd(1,complianceProb,length(comArray),1);
tooLong = 0;

%% Begin simulation
 while  t{k}(tStep)<tMax && tStep<tStepMax



    lambdaRsd = lambdaVectorSumRsd + lambdaInfectionSumRsd ;
    lambdaCom = lambdaVectorSumCom + lambdaInfectionSumCom ;
    lambda = lambdaRsd + lambdaCom;

    randNum=rand(2,1);
    tau = log(1/randNum(1)) / lambda;


    if tau > (nYrs+1)*365 || lambda<=0
        tooLong = 1;
        tau = (nYrs+1)*365;
    end       

    t{k}(tStep+1) = t{k}(tStep) + tau;

    % DETECTION:
    if firstDetection(k)==0
    if floor((t{k}(tStep+1)-firstInspection)/(365*inspInt(1))) - floor((t{k}(tStep)-firstInspection)/(365*inspInt(1))) >0
        %For each inspection interval that has occured before the next event:
        for i=1:floor((t{k}(tStep+1)-firstInspection)/(365*inspInt(1))) - floor((t{k}(tStep)-firstInspection)/(365*inspInt(1)))
            event{k}(tStep,2)=15;
            rsdInsp = randsample(Cell(rsdArray>0),min(sum(rsdArray>0),nInspRsd(1)));
            comInsp = randsample(Cell(comArray>0),min(sum(comArray>0),nInspCom(1)));

            % HOW TO MODEL SAMPLING, CAN UNITS BE REPEATED?
            rsdFound = hygernd(rsdArray(rsdInsp),Irsd(rsdInsp),min(nSamp(1),rsdArray(rsdInsp)));
            comFound = hygernd(comArray(comInsp),Icom(comInsp),min(nSamp(1),comArray(comInsp)));
            rsdDetected = binornd(rsdFound,detProb);
            comDetected = binornd(comFound,detProb);

            if sum(rsdDetected)+sum(comDetected)>0

               
                firstDetection(k)=ceil((t{k}(tStep)-firstInspection)/(365*inspInt(1)))*inspInt(1)*365+firstInspection + (i-1)*inspInt(1);
                detectionCitrusIncidence(k) = (sum(Icom) + sum(Ccom) + sum(Ecom))/sum(comArray);
                detectionCellIncidence(k) = (sum(Icom + Ccom + Ecom>0))/sum(comArray>0);
                t{k}(tStep+1) = firstDetection(k);
                event{k}(tStep,2) = 14;                    
                
                rsdRogued = binornd(rsdDetected,roguingProb);

                if sum(rsdRogued)>0

                    event{k}(tStep,1)=1;
                    cellArray=rsdInsp(rsdRogued>0);
                    cellArray = cellArray(complianceArray(cellArray)==1);

                    % Remove infected trees 
                    Irsd(cellArray) = Irsd(cellArray) - rsdRogued(ismember(rsdInsp,cellArray));
                    Rrsd(cellArray) = Rrsd(cellArray) + rsdRogued(ismember(rsdInsp,cellArray));
                    
                    % and update rates 
                    cellArray=rsdInsp(rsdRogued>0);
                    cellArray = cellArray(complianceArray(cellArray)==1);
                    for c=1:length(cellArray)
                        cell=cellArray(c);
                        newR = rsdRogued(rsdInsp==cell);
                    [forceInfectRsd,forceVectorRsd,longDistDispRsd,estabRsd,...
                    latentRateRsd,symptomRateRsd,forceVectorCom,forceInfectCom,...
                    lambdaInfectionSumRsd,lambdaVectorSumRsd,lambdaInfectionSumCom,lambdaVectorSumCom, rsdVecTime(:,k),rsdInfTime(:,k),adjNrsd] = ...
                    updateRates(Nrsd,Ersd,Crsd,Irsd,Rrsd,newR,fullInfestation,forceInfectRsd,forceVectorRsd,0,M,...
                    forceVectorCom,forceInfectCom,lambdaInfectionSumRsd,lambdaVectorSumRsd,lambdaInfectionSumCom,lambdaVectorSumCom,...
                    longDistDispRsd,estabRsd,longDistDispCom,0,0,extInfCSum,extVCSum,estabCom,...
                    latentRateRsd,symptomRateRsd,latentRateCom,symptomRateCom,beta,localDispRate,longDispRate,estRate,lRate,sRate,rsdArray,adjNrsd,climArray,rsdCC,...
                    comCC,Ncom,Ecom,Ccom,Icom,Rcom,comArray,rsdVecTime(:,k),rsdInfTime(:,k),Fi,Kd,Ki,cell,event{k}(tStep,2),t{k}(tStep+1));
                    end
                end
                    
                comRogued = binornd(comDetected,roguingProb);

                if sum(comRogued)>0

                    event{k}(tStep,1)=2;
                    cellArray=comInsp(comRogued>0);
                    cellArray = cellArray(complianceArray(cellArray)==1);
                    
                    % Remove infected trees 
                    Icom(cellArray) = Icom(cellArray) - comRogued(ismember(comInsp,cellArray));
                    Rcom(cellArray) = Rcom(cellArray) + comRogued(ismember(comInsp,cellArray));
                    
                    % and update rates 
                    for c=1:length(cellArray)
                        cell=cellArray(c);
                        newR = comRogued(comInsp==cell);
                    [forceInfectCom,forceVectorCom,longDistDispCom,estabCom,...
                    latentRateCom,symptomRateCom,forceVectorRsd,forceInfectRsd,...
                    lambdaInfectionSumCom,lambdaVectorSumCom,lambdaInfectionSumRsd,lambdaVectorSumRsd,comVecTime(:,k),comInfTime(:,k),adjNcom] = ...
                    updateRates(Ncom,Ecom,Ccom,Icom,Rcom,newR,fullInfestation,forceInfectCom,forceVectorCom,1,M,...
                    forceVectorRsd,forceInfectRsd,lambdaInfectionSumCom,lambdaVectorSumCom,lambdaInfectionSumRsd,lambdaVectorSumRsd,...
                    longDistDispCom,estabCom,longDistDispRsd,extInfCSum,extVCSum,0,0,estabRsd,...
                    latentRateCom,symptomRateCom,latentRateRsd,symptomRateRsd,beta,localDispRate,longDispRate,estRate,lRate,sRate,...
                    comArray,adjNcom,climArray,comCC,rsdCC,Nrsd,Ersd,Crsd,Irsd,Rrsd,rsdArray,...
                    comVecTime(:,k),comInfTime(:,k),Fi,Kd,Ki,cell,event{k}(tStep,2),t{k}(tStep+1));
                    end
                end

                if chemicalIncrease >0
                    M_old = M;
                    adjNcom_old = adjNcom;

                    m = 1-(1-m)*(1-chemicalIncrease);
                    M = 1-(m*(1-(abnArray+orgArray)./realComCit) + (1-paraArray).*(abnArray+orgArray)./realComCit);
                    adjNcom = realComCit.*climArray.*M./comCC; adjNcom=max(0,adjNcom); 

                    %UPDATE RATES
                    longDistDispCom = longDistDispCom .* max(0,M ./ M_old);

                    fIC_old = forceInfectCom;
                    Scom = comArray- Ecom-Ccom-Icom-Rcom;
                    Srsd = rsdArray- Ersd-Crsd-Irsd-Rrsd;
                    for c=1:length(comArray)
                        if fullInfestation~=1
                            forceVectorCom(c) = forceVectorCom(c)- forceVectorCom(c)*chemicalIncrease*max(0,sum(Kd{c}.*Ncom(Fi{c}).*adjNcom_old(Fi{c}).* (1-(abnArray(Fi{c})+orgArray(Fi{c}))./comArray(Fi{c})) )/sum(Kd{c}.*(adjNcom_old(Fi{c}).*Ncom(Fi{c})+adjNrsd(Fi{c}).*Nrsd(Fi{c})))) ;
                            forceVectorRsd(c) = forceVectorRsd(c)- forceVectorRsd(c)*chemicalIncrease*max(0,sum(Kd{c}.*Ncom(Fi{c}).*adjNcom_old(Fi{c}).* (1-(abnArray(Fi{c})+orgArray(Fi{c}))./comArray(Fi{c})) )/sum(Kd{c}.*(adjNcom_old(Fi{c}).*Ncom(Fi{c})+adjNrsd(Fi{c}).*Nrsd(Fi{c})))) ;
                        end

                        forceInfectCom(c) = beta * (nansum(Ki{c}.*( ((Ccom(Fi{c})+Icom(Fi{c}))./comArray(Fi{c})).*adjNcom(Fi{c})))...
                            +nansum(Ki{c}.*((Crsd(Fi{c})+Irsd(Fi{c}))./rsdArray(Fi{c})).*adjNrsd(Fi{c})))* max(0,Scom(c)./(comArray(c))) ;

                        forceInfectRsd(c) = beta * (nansum(Ki{c}.*( ((Ccom(Fi{c})+Icom(Fi{c}))./comArray(Fi{c})).*adjNcom(Fi{c})))...
                            +nansum(Ki{c}.*((Crsd(Fi{c})+Irsd(Fi{c}))./rsdArray(Fi{c})).*adjNrsd(Fi{c})))* max(0,Srsd(c)./(rsdArray(c))) ;
                    end
                    lambdaVectorSumCom = sum(forceVectorCom) + sum(longDistDispCom) + extVCSum(end) + sum(estabCom);
                    lambdaVectorSumRsd = sum(forceVectorRsd) + sum(longDistDispRsd) + sum(estabRsd);
                    lambdaInfectionSumCom = sum(forceInfectCom) + extInfCSum(end) + sum(latentRateCom) + sum(symptomRateCom);
                    lambdaInfectionSumRsd = sum(forceInfectRsd) + sum(latentRateRsd) + sum(symptomRateRsd);

                end

                if floor(t{k}(tStep+1)*spyr/(365)) - floor(t{k}(tStep)*spyr/(365)) >0
                    for ii=1:floor(t{k}(tStep+1)*spyr/(365)) - floor(t{k}(tStep)*spyr/(365))
                        intNum = floor(t{k}(tStep)*spyr/(365))+ii ;
                        if intNum<=nYrs*spyr
                            ArsdInt(:,intNum+1,k) = Arsd; AcomInt(:,intNum+1,k) = Acom;
                            NrsdInt(:,intNum+1,k) = Nrsd; NcomInt(:,intNum+1,k) = Ncom;
                            ErsdInt(:,intNum+1,k) = Ersd; EcomInt(:,intNum+1,k) = Ecom;
                            CrsdInt(:,intNum+1,k) = Crsd; CcomInt(:,intNum+1,k) = Ccom;
                            IrsdInt(:,intNum+1,k) = Irsd; IcomInt(:,intNum+1,k) = Icom;
                            RrsdInt(:,intNum+1,k) = Rrsd; RcomInt(:,intNum+1,k) = Rcom;
                        end
                    end
                end

                tStep=tStep+1;
                lambdaRsd = lambdaVectorSumRsd + lambdaInfectionSumRsd ;
                lambdaCom = lambdaVectorSumCom + lambdaInfectionSumCom ;
                lambda = lambdaRsd + lambdaCom;
                randNum=rand(2,1);
                tau = log(1./randNum(1)) ./ lambda;
                if tau > (nYrs+1)*365 || lambda<=0
                    tooLong = 1;
                end       
               
                t{k}(tStep+1) = t{k}(tStep) + tau;
        
            end
        end
    end
    else 
        if firstDetection(k)>0
            if floor((t{k}(tStep+1)-firstInspection)/(365*inspInt(2))) - floor((t{k}(tStep)-firstInspection)/(365*inspInt(2))) >0
                %For each inspection interval that has occured before the next event:
                for i=1:floor((t{k}(tStep+1)-firstInspection)/(365*inspInt(2))) - floor((t{k}(tStep)-firstInspection)/(365*inspInt(2)))
                    event{k}(tStep,2)=15;
                    rsdInsp = randsample(Cell(rsdArray>0),min(sum(rsdArray>0),nInspRsd(2)));
                    comInsp = randsample(Cell(comArray>0),min(sum(comArray>0),nInspCom(2)));
        
                    % HOW TO MODEL SAMPLING, CAN UNITS BE REPEATED?
                    rsdFound = hygernd(rsdArray(rsdInsp)-Rrsd(rsdInsp),Irsd(rsdInsp),min(nSamp(2),rsdArray(rsdInsp)-Rrsd(rsdInsp)));
                    comFound = hygernd(comArray(comInsp)-Rcom(comInsp),Icom(comInsp),min(nSamp(2),comArray(comInsp)-Rcom(comInsp)));
                    rsdDetected = binornd(rsdFound,detProb);
                    comDetected = binornd(comFound,detProb);

                    rsdRogued = binornd(rsdDetected,roguingProb);
                    comRogued = binornd(comDetected,roguingProb);
                    if sum(rsdRogued)+sum(comRogued)>0
                        
                        t{k}(tStep+1) = ceil((t{k}(tStep)-firstInspection)/(365*inspInt(2)))*inspInt(2)*365+firstInspection + (i-1)*inspInt(2);
                        event{k}(tStep,2)=14;
                        
                        if sum(rsdRogued)>0

                            event{k}(tStep,1)=1;
                            cellArray=rsdInsp(rsdRogued>0);
                            cellArray = cellArray(complianceArray(cellArray)==1);

                            % Remove infected trees 
                            Irsd(cellArray) = Irsd(cellArray) - rsdRogued(ismember(rsdInsp,cellArray));
                            Rrsd(cellArray) = Rrsd(cellArray) + rsdRogued(ismember(rsdInsp,cellArray));
                            
                            % and update rates 
                            for c=1:length(cellArray)
                                cell=cellArray(c);
                                newR = rsdRogued(rsdInsp==cell);
                            [forceInfectRsd,forceVectorRsd,longDistDispRsd,estabRsd,...
                            latentRateRsd,symptomRateRsd,forceVectorCom,forceInfectCom,...
                            lambdaInfectionSumRsd,lambdaVectorSumRsd,lambdaInfectionSumCom,lambdaVectorSumCom, rsdVecTime(:,k),rsdInfTime(:,k),adjNrsd] = ...
                            updateRates(Nrsd,Ersd,Crsd,Irsd,Rrsd,newR,fullInfestation,forceInfectRsd,forceVectorRsd,0,M,...
                            forceVectorCom,forceInfectCom,lambdaInfectionSumRsd,lambdaVectorSumRsd,lambdaInfectionSumCom,lambdaVectorSumCom,...
                            longDistDispRsd,estabRsd,longDistDispCom,0,0,extInfCSum,extVCSum,estabCom,...
                            latentRateRsd,symptomRateRsd,latentRateCom,symptomRateCom,beta,localDispRate,longDispRate,estRate,lRate,sRate,rsdArray,adjNrsd,climArray,rsdCC,...
                            comCC,Ncom,Ecom,Ccom,Icom,Rcom,comArray,rsdVecTime(:,k),rsdInfTime(:,k),Fi,Kd,Ki,cell,event{k}(tStep,2),t{k}(tStep+1));
                            end
                        end
        
                        if sum(comRogued)>0

                            event{k}(tStep,1)=2;
                            cellArray=comInsp(comRogued>0);
                            cellArray = cellArray(complianceArray(cellArray)==1);

                            % Remove infected trees 
                            Icom(cellArray) = Icom(cellArray) - comRogued(ismember(comInsp,cellArray));
                            Rcom(cellArray) = Rcom(cellArray) + comRogued(ismember(comInsp,cellArray));
                            
                            % and update rates 
                            for c=1:length(cellArray)
                                cell=cellArray(c);
                                newR = comRogued(comInsp==cell);
                            [forceInfectCom,forceVectorCom,longDistDispCom,estabCom,...
                            latentRateCom,symptomRateCom,forceVectorRsd,forceInfectRsd,...
                            lambdaInfectionSumCom,lambdaVectorSumCom,lambdaInfectionSumRsd,lambdaVectorSumRsd,comVecTime(:,k),comInfTime(:,k),adjNcom] = ...
                            updateRates(Ncom,Ecom,Ccom,Icom,Rcom,newR,fullInfestation,forceInfectCom,forceVectorCom,1,M,...
                            forceVectorRsd,forceInfectRsd,lambdaInfectionSumCom,lambdaVectorSumCom,lambdaInfectionSumRsd,lambdaVectorSumRsd,...
                            longDistDispCom,estabCom,longDistDispRsd,extInfCSum,extVCSum,0,0,estabRsd,...
                            latentRateCom,symptomRateCom,latentRateRsd,symptomRateRsd,beta,localDispRate,longDispRate,estRate,lRate,sRate,...
                            comArray,adjNcom,climArray,comCC,rsdCC,Nrsd,Ersd,Crsd,Irsd,Rrsd,rsdArray,...
                            comVecTime(:,k),comInfTime(:,k),Fi,Kd,Ki,cell,event{k}(tStep,2),t{k}(tStep+1));
                            end

                        end
                        if floor(t{k}(tStep+1)*spyr/(365)) - floor(t{k}(tStep)*spyr/(365)) >0
                            for ii=1:floor(t{k}(tStep+1)*spyr/(365)) - floor(t{k}(tStep)*spyr/(365))
                                intNum = floor(t{k}(tStep)*spyr/(365))+ii ;
                                if intNum<=nYrs*spyr
                                    ArsdInt(:,intNum+1,k) = Arsd; AcomInt(:,intNum+1,k) = Acom;
                                    NrsdInt(:,intNum+1,k) = Nrsd; NcomInt(:,intNum+1,k) = Ncom;
                                    ErsdInt(:,intNum+1,k) = Ersd; EcomInt(:,intNum+1,k) = Ecom;
                                    CrsdInt(:,intNum+1,k) = Crsd; CcomInt(:,intNum+1,k) = Ccom;
                                    IrsdInt(:,intNum+1,k) = Irsd; IcomInt(:,intNum+1,k) = Icom;
                                    RrsdInt(:,intNum+1,k) = Rrsd; RcomInt(:,intNum+1,k) = Rcom;
                                end
                            end
                        end
        
                        tStep=tStep+1;
                        lambdaRsd = lambdaVectorSumRsd + lambdaInfectionSumRsd ;
                        lambdaCom = lambdaVectorSumCom + lambdaInfectionSumCom ;
                        lambda = lambdaRsd + lambdaCom;
                        randNum=rand(2,1);
                        tau = log(1./randNum(1)) ./ lambda;
                        t{k}(tStep+1) = t{k}(tStep) + tau;
                        if tau > (nYrs+1)*365 || lambda<=0
                            tooLong = 1;
                        end       
                
                    
                    end

                end
            end
        end
    end

   if (stopAtDetection==1 && firstDetection(k)>0) || tooLong==1
        ArsdInt(:,floor(t{k}(tStep)*spyr/(365))+1:nYrs*spyr+1,k) = Arsd*ones(1,nYrs*spyr+1-floor(t{k}(tStep)*spyr/(365))); AcomInt(:,floor(t{k}(tStep)*spyr/(365))+1:nYrs*spyr+1,k) = Acom*ones(1,nYrs*spyr+1-floor(t{k}(tStep)*spyr/(365)));
        NrsdInt(:,floor(t{k}(tStep)*spyr/(365))+1:nYrs*spyr+1,k) = Nrsd*ones(1,nYrs*spyr+1-floor(t{k}(tStep)*spyr/(365))); NcomInt(:,floor(t{k}(tStep)*spyr/(365))+1:nYrs*spyr+1,k) = Ncom*ones(1,nYrs*spyr+1-floor(t{k}(tStep)*spyr/(365)));
        ErsdInt(:,floor(t{k}(tStep)*spyr/(365))+1:nYrs*spyr+1,k) = Ersd*ones(1,nYrs*spyr+1-floor(t{k}(tStep)*spyr/(365))); EcomInt(:,floor(t{k}(tStep)*spyr/(365))+1:nYrs*spyr+1,k) = Ecom*ones(1,nYrs*spyr+1-floor(t{k}(tStep)*spyr/(365)));
        CrsdInt(:,floor(t{k}(tStep)*spyr/(365))+1:nYrs*spyr+1,k) = Crsd*ones(1,nYrs*spyr+1-floor(t{k}(tStep)*spyr/(365))); CcomInt(:,floor(t{k}(tStep)*spyr/(365))+1:nYrs*spyr+1,k) = Ccom*ones(1,nYrs*spyr+1-floor(t{k}(tStep)*spyr/(365)));
        IrsdInt(:,floor(t{k}(tStep)*spyr/(365))+1:nYrs*spyr+1,k) = Irsd*ones(1,nYrs*spyr+1-floor(t{k}(tStep)*spyr/(365))); IcomInt(:,floor(t{k}(tStep)*spyr/(365))+1:nYrs*spyr+1,k) = Icom*ones(1,nYrs*spyr+1-floor(t{k}(tStep)*spyr/(365)));
        RrsdInt(:,floor(t{k}(tStep)*spyr/(365))+1:nYrs*spyr+1,k) = Rrsd*ones(1,nYrs*spyr+1-floor(t{k}(tStep)*spyr/(365))); RcomInt(:,floor(t{k}(tStep)*spyr/(365))+1:nYrs*spyr+1,k) = Rcom*ones(1,nYrs*spyr+1-floor(t{k}(tStep)*spyr/(365)));
       break
   end

% GENERATE NEXT STEP:
    if randNum(2)*lambda < lambdaRsd
        event{k}(tStep,1) = 1;

        [Nrsd,Arsd,Ersd,Crsd,Irsd,Rrsd,Ncom,Acom,Ecom,cell,event{k}(tStep,2)] = ...
            nextEvent(Nrsd,Arsd,Ersd,Crsd,Irsd,Rrsd,Ncom,Acom,Ecom,Ccom,Icom,Rcom,forceInfectRsd,forceVectorRsd,0,0,...
            lambdaInfectionSumRsd,lambdaVectorSumRsd,longDistDispRsd,estabRsd,latentRateRsd,symptomRateRsd,...
            alphaLD,rsdArray,comArray,rsdCC,comCC,climArray,M,0,coord,citNums,gridX,gridY,fullInfestation);

    else 
        if randNum(2)*lambda < lambdaRsd + lambdaCom
            event{k}(tStep,1) = 2;
            [Ncom,Acom,Ecom,Ccom,Icom,Rcom,Nrsd,Arsd,Ersd,cell,event{k}(tStep,2)] =...
                nextEvent(Ncom,Acom,Ecom,Ccom,Icom,Rcom,Nrsd,Arsd,Ersd,Crsd,Irsd,Rrsd,forceInfectCom,forceVectorCom,extInfCSum,extVCSum,...
                lambdaInfectionSumCom,lambdaVectorSumCom,longDistDispCom,estabCom,latentRateCom,symptomRateCom,...
                alphaLD,comArray,rsdArray,comCC,rsdCC,climArray,M,1,coord,citNums,gridX,gridY,fullInfestation);


        end
    end

    % UPDATE RATES AFTER EVENT:
    if event{k}(tStep,2)>0
        if (event{k}(tStep,1)==1 && sum(event{k}(tStep,2)-[5,6,7]==0)==0) || (event{k}(tStep,1)==2 && sum(event{k}(tStep,2)-[5,6,7]==0)>0)


       
            [forceInfectRsd,forceVectorRsd,longDistDispRsd,estabRsd,...
                latentRateRsd,symptomRateRsd,forceVectorCom,forceInfectCom,...
                lambdaInfectionSumRsd,lambdaVectorSumRsd,lambdaInfectionSumCom,lambdaVectorSumCom, rsdVecTime(:,k),rsdInfTime(:,k),adjNrsd] = ...
                updateRates(Nrsd,Ersd,Crsd,Irsd,Rrsd,0,fullInfestation,forceInfectRsd,forceVectorRsd,0,M,...
                forceVectorCom,forceInfectCom,lambdaInfectionSumRsd,lambdaVectorSumRsd,lambdaInfectionSumCom,lambdaVectorSumCom,...
                longDistDispRsd,estabRsd,longDistDispCom,0,0,extInfCSum,extVCSum,estabCom,...
                latentRateRsd,symptomRateRsd,latentRateCom,symptomRateCom,beta,localDispRate,longDispRate,estRate,lRate,sRate,rsdArray,adjNrsd,climArray,rsdCC,...
                comCC,Ncom,Ecom,Ccom,Icom,Rcom,comArray,rsdVecTime(:,k),rsdInfTime(:,k),Fi,Kd,Ki,cell,event{k}(tStep,2),t{k}(tStep+1));

        else
            if (event{k}(tStep,1)==2 && sum(event{k}(tStep,2)-[5,6,7]==0)==0) || (event{k}(tStep,1)==1 && sum(event{k}(tStep,2)-[5,6,7]==0)>0)
            

                [forceInfectCom,forceVectorCom,longDistDispCom,estabCom,...
                    latentRateCom,symptomRateCom,forceVectorRsd,forceInfectRsd,...
                    lambdaInfectionSumCom,lambdaVectorSumCom,lambdaInfectionSumRsd,lambdaVectorSumRsd,comVecTime(:,k),comInfTime(:,k),adjNcom] = ...
                    updateRates(Ncom,Ecom,Ccom,Icom,Rcom,0,fullInfestation,forceInfectCom,forceVectorCom,1,M,...
                    forceVectorRsd,forceInfectRsd,lambdaInfectionSumCom,lambdaVectorSumCom,lambdaInfectionSumRsd,lambdaVectorSumRsd,...
                    longDistDispCom,estabCom,longDistDispRsd,extInfCSum,extVCSum,0,0,estabRsd,...
                    latentRateCom,symptomRateCom,latentRateRsd,symptomRateRsd,beta,localDispRate,longDispRate,estRate,lRate,sRate,...
                    comArray,adjNcom,climArray,comCC,rsdCC,Nrsd,Ersd,Crsd,Irsd,Rrsd,rsdArray,...
                    comVecTime(:,k),comInfTime(:,k),Fi,Kd,Ki,cell,event{k}(tStep,2),t{k}(tStep+1));

            end
        end
    end


    % Track number of infested and infected cells at each daily interval
    if floor(t{k}(tStep+1)) - floor(t{k}(tStep))>0
        numInfested{k}(floor(t{k}(tStep+1))) = sum(Nrsd) + sum(Ncom);% + sum(Arsd) + sum(Acom);
        citInfected{k}(floor(t{k}(tStep+1))) = sum(Irsd) + sum(Icom) + sum(Crsd) + sum(Ccom) + sum(Ersd) + sum(Ecom) + sum(Rrsd) + sum(Rcom);
        cellInfected{k}(floor(t{k}(tStep+1))) = sum(Irsd+Icom+Crsd+Ccom+Ersd+Ecom+Rrsd+Rcom>0) ;
    end

    % Track individual status of cells at each interval (spyr - saves per year)
    if floor(t{k}(tStep+1)*spyr/(365)) - floor(t{k}(tStep)*spyr/(365)) >0
        for i=1:floor(t{k}(tStep+1)*spyr/(365)) - floor(t{k}(tStep)*spyr/(365))
            intNum = floor(t{k}(tStep)*spyr/(365))+i ;
            if intNum<=nYrs*spyr
                ArsdInt(:,intNum+1,k) = Arsd; AcomInt(:,intNum+1,k) = Acom;
                NrsdInt(:,intNum+1,k) = Nrsd; NcomInt(:,intNum+1,k) = Ncom;
                ErsdInt(:,intNum+1,k) = Ersd; EcomInt(:,intNum+1,k) = Ecom;
                CrsdInt(:,intNum+1,k) = Crsd; CcomInt(:,intNum+1,k) = Ccom;
                IrsdInt(:,intNum+1,k) = Irsd; IcomInt(:,intNum+1,k) = Icom;
                RrsdInt(:,intNum+1,k) = Rrsd; RcomInt(:,intNum+1,k) = Rcom;
            end
        end
    end


    if ~isnan(cell)
    cellnum(tStep)=cell;
    xpos(tStep) = x(cell); ypos(tStep) = y(cell);
    end
    tStep=tStep+1;
end

t{k}=t{k}(2:tStep); event{k}=event{k}(1:tStep-1,:);
infEvent = event{k}(ismember(event{k}(:,2),[3,4,6,7,10,11,12,13]) ,:);
infEvent(ismember(infEvent(:,2),[3,4,6,7,10,11]) ,2) = 1;
infEvent(infEvent(:,2)==12 ,2) = 2;
infEvent(infEvent(:,2)==13 ,2) = 3;
eventtime = t{k}(ismember(event{k}(:,2),[3,4,6,7,10,11,12,13]))'/365;
cellnum = cellnum(ismember(event{k}(:,2),[3,4,6,7,10,11,12,13]));
ypos = ypos(ismember(event{k}(:,2),[3,4,6,7,10,11,12,13]));
xpos = xpos(ismember(event{k}(:,2),[3,4,6,7,10,11,12,13]));
if infected0>0
    eventtime=[zeros(length(infected0),1); eventtime];
    infEvent = [[2 1].*ones(length(infected0),1); infEvent];
    cellnum = [infected0; cellnum];
    xpos = [x(infected0); xpos];
    ypos = [y(infected0); ypos];
end
InfectionTimes=table(eventtime,cellnum,xpos,ypos,infEvent);
writetable(InfectionTimes,append('ModelOutputs/',name,'_',num2str(k),'_InfectionTimes.txt'),'Delimiter','tab')  

numInfested{k}(tMax) = sum(Nrsd) + sum(Ncom);% + sum(Arsd) + sum(Acom);
citInfected{k}(tMax) = sum(Irsd) + sum(Icom) + sum(Crsd) + sum(Ccom) + sum(Ersd) + sum(Ecom) + sum(Rrsd) + sum(Rcom);

cellInfected{k}(tMax) = sum(Irsd+Icom+Crsd+Ccom+Ersd+Ecom+Rrsd+Rcom>0) ;

time=toc;
if mod(k,1)==0
    disp(['Rep complete: ',num2str(k),', rep time = ',num2str(time/60),' mins'])
end

eventCount(:,:,k)=[hist(event{k}(event{k}(:,1)==1,2),-2:13)' hist(event{k}(event{k}(:,1)==2,2),-2:13)'];
residentialOutputs(k,:) = [sum(Nrsd), sum(Nrsd>0), sum(Irsd), sum(Irsd>0), time, tStep];
commercialOutputs(k,:) = [sum(Ncom), sum(Ncom>0), sum(Icom), sum(Icom>0), time, tStep];

ArsdInt(:,nYrs*spyr+1,k) = Arsd; AcomInt(:,nYrs*spyr+1,k) = Acom;
NrsdInt(:,nYrs*spyr+1,k) = Nrsd; NcomInt(:,nYrs*spyr+1,k) = Ncom;
ErsdInt(:,nYrs*spyr+1,k) = Ersd; EcomInt(:,nYrs*spyr+1,k) = Ecom;
CrsdInt(:,nYrs*spyr+1,k) = Crsd; CcomInt(:,nYrs*spyr+1,k) = Ccom;
IrsdInt(:,nYrs*spyr+1,k) = Irsd; IcomInt(:,nYrs*spyr+1,k) = Icom;
RrsdInt(:,nYrs*spyr+1,k) = Rrsd; RcomInt(:,nYrs*spyr+1,k) = Rcom;


Time = (0:1/spyr:nYrs)';
rsdOutputs = [sum(ArsdInt(:,:,k))' sum(NrsdInt(:,:,k))' sum(ErsdInt(:,:,k))' sum(CrsdInt(:,:,k))' sum(IrsdInt(:,:,k))'];
comOutputs = [sum(AcomInt(:,:,k))' sum(NcomInt(:,:,k))' sum(EcomInt(:,:,k))' sum(CcomInt(:,:,k))' sum(IcomInt(:,:,k))'];
OutputsTab = splitvars(table(Time,rsdOutputs,comOutputs));
OutputsTab.Properties.VariableNames = ["Time","Arsd","Nrsd","Lrsd","Crsd","Irsd","Acom","Ncom","Lcom","Ccom","Icom"];
writetable(OutputsTab,append('ModelOutputs/',name,'_Outputs_',num2str(k),'.txt'),'Delimiter','tab')  

vecTime(:,k) = min(rsdVecTime(:,k),comVecTime(:,k));
infTime(:,k) = min(rsdInfTime(:,k),comInfTime(:,k));

rsd_infecttime = rsdInfTime(:,k)/365;
com_infecttime = comInfTime(:,k)/365;
rsd_infesttime = rsdVecTime(:,k)/365;
com_infesttime = comVecTime(:,k)/365;
com_cit = comArray; rsd_cit = rsdArray;
CellInfections=table(Cell,x,y,rsd_cit,com_cit,rsd_infecttime,com_infecttime,Ersd,Crsd,Irsd,Ecom,Ccom,Icom);
writetable(CellInfections,append('ModelOutputs/',name,'_',num2str(k),'_CellInfections.txt'),'Delimiter','tab')  


end

%%

for k=1:nRuns
    nInf(:,k)=numInfested{k}(1:tMax) ;
    citInf(:,k)=citInfected{k}(1:tMax);
    cellInf(:,k)=cellInfected{k}(1:tMax);
    while sum(isnan(citInf(365*10,k)),'all')>0
        citInf(isnan(citInf(:,k)),k)=citInf(find(isnan(citInf(:,k)))-1,k);
    end
end


nInf(nInf==0)=nan;
if fullInfestation~=1
nInf(1,:) = ones(1,nRuns)*0;%size(infested0,1);
else
    nInf(1,:) = ones(1,nRuns)* (sum(rsdArray>0)+sum(comArray>0));
end

%% Export inputs and outputs


cellRsdOutputs = [[sum(NrsdInt(:,end,:)>0,3) sum(NrsdInt(:,end,:)+ArsdInt(:,end,:)>0,3)...
    sum(IrsdInt(:,end,:)>0,3) sum(IrsdInt(:,end,:)+CrsdInt(:,end,:)+ErsdInt(:,end,:)>0,3)]/nRuns...
    [mean(NrsdInt(:,end,:),3,'omitnan') mean(IrsdInt(:,end,:),3,'omitnan') mean(IrsdInt(:,end,:)+CrsdInt(:,end,:)+ErsdInt(:,end,:),3,'omitnan')...
    mean(rsdVecTime,2,'omitnan') mean(rsdInfTime,2,'omitnan')]];
cellComOutputs = [[sum(NcomInt(:,end,:)>0,3) sum(NcomInt(:,end,:)+AcomInt(:,end,:)>0,3)...
    sum(IcomInt(:,end,:)>0,3) sum(IcomInt(:,end,:)+CcomInt(:,end,:)+EcomInt(:,end,:)>0,3)]/nRuns...
    [mean(NcomInt(:,end,:),3,'omitnan') mean(IcomInt(:,end,:),3,'omitnan') mean(IcomInt(:,end,:)+CcomInt(:,end,:)+EcomInt(:,end,:),3,'omitnan')...
    mean(comVecTime,2,'omitnan') mean(comInfTime,2,'omitnan')]];
cellTotOutputs= [[sum(NcomInt(:,end,:)+NrsdInt(:,end,:)>0,3) sum(NcomInt(:,end,:)+AcomInt(:,end,:)+NrsdInt(:,end,:)+ArsdInt(:,end,:)>0,3)...
    sum(IcomInt(:,end,:)+IrsdInt(:,end,:)>0,3) sum(IcomInt(:,end,:)+CcomInt(:,end,:)+EcomInt(:,end,:)+IrsdInt(:,end,:)+CrsdInt(:,end,:)+ErsdInt(:,end,:)>0,3)]/nRuns...
    [mean(NcomInt(:,end,:)+NrsdInt(:,end,:),3,'omitnan') mean(IcomInt(:,end,:)+IrsdInt(:,end,:),3,'omitnan') mean(IcomInt(:,end,:)+CcomInt(:,end,:)+EcomInt(:,end,:)+IrsdInt(:,end,:)+CrsdInt(:,end,:)+ErsdInt(:,end,:),3,'omitnan')...
    mean(min(comVecTime,rsdVecTime),2,'omitnan') mean(min(comInfTime,rsdInfTime),2,'omitnan')]];
cellOutputsTab = splitvars(table(Cell,x,y,cellTotOutputs,cellRsdOutputs,cellComOutputs));
OutputsTab.Properties.VariableNames = ["Cell","x","y","Vrsd cells","V+A cells","I cells","I+C+L cells","Vcom","V+A","I","I+C+L"];
writetable(cellOutputsTab,append('ModelOutputs/',name,'_CellOutputs.txt'),'Delimiter','tab')  
 




%% Plot figures

if sum([mapFigure,multiFigure,individualFigures])>0
    if not(isfolder("Figures"))
        mkdir("Figures")
    end
end

if mapFigure>0
    plotMaps(citMat,ocean,rsdArray,comArray,name);
end

if multiFigure>0
    plotMultiTotals(citMat,ocean,rsdArray,comArray,NrsdInt(:,end,:),NcomInt(:,end,:),...
        IrsdInt(:,end,:),IcomInt(:,end,:),ErsdInt(:,end,:),EcomInt(:,end,:),CrsdInt(:,end,:),CcomInt(:,end,:),...
        nInf,citInf,cellInf,rsdInfTime,comInfTime,rsdVecTime,comVecTime,tMax,name)
    plotMultiStoryboard(citMat,ocean,spyr,rsdArray,comArray,NrsdInt,NcomInt,...
        IrsdInt,IcomInt,ErsdInt,EcomInt,CrsdInt,CcomInt,...
        nInf,citInf,cellInf,rsdInfTime,comInfTime,rsdVecTime,comVecTime,tMax,name)
end

if individualFigures>0
    if not(isfolder("Figures/SingleSims"))
        mkdir("Figures/SingleSims")
    end
    
    for k=1:min(nRuns,nIndividualSimFigs)
        plotStoryboard(k,citMat,spyr,nYrs,rsdArray,comArray,ocean,ErsdInt(:,:,k),EcomInt(:,:,k),CrsdInt(:,:,k),CcomInt(:,:,k),...
             IrsdInt(:,:,k),IcomInt(:,:,k),RrsdInt(:,:,k),RcomInt(:,:,k),name)
        plotVectorStoryboard(k,citMat,spyr,nYrs,rsdArray,comArray,ocean,NrsdInt(:,:,k),NcomInt(:,:,k),adjNrsd,adjNcom,name)
    end
end
