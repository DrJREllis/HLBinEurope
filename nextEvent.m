function [N,A,E,C,I,R,Noth,Aoth,Eoth,cell,event] = nextEvent(N,A,E,C,I,R,Noth,Aoth,Eoth,Coth,Ioth,Roth,forceInfect,forceVector,...
    extInfCSum,extVCSum,lambdaInfectionSum,lambdaVectorSum,ldd,estab,latentRate,symptomRate,...
    alphaLD,citArray,othArray,carrCap,othCC,climArray,M,comInd,coord,citNums,gridX,gridY,fullInfestation)
% Function to generate the next event in the gillespie algorithm

lambda = lambdaVectorSum + lambdaInfectionSum ;

randNum = rand(1,3);
randLam = randNum(1)*lambda;

event=0;
cell=0;


forceVCSum = cumsum(forceVector); 

 %SUSCEPTIBLE DISPERSAL EVENT
if randLam < forceVCSum(end)
    event = 1;

    cell = find(randNum(3)*forceVCSum(end)<forceVCSum,1); 

    A(cell) = 1;
 

else % LONG DISTANCE DISPERSAL
    if randLam < forceVCSum(end) + sum(ldd)

        onLS = false;
        while(onLS == false)
            oCell = find(randNum(2)*sum(ldd)<cumsum(ldd),1); 
            dist = trnd(3)*alphaLD;
            angle = rand*2*pi;
    
            movX = dist*cos(angle); %step length in x direction
            movY = dist*sin(angle); %step length in y direction
            Ix=coord(oCell,1); Iy=coord(oCell,2);
    
            newX = Ix+round(movX); newY = Iy+round(movY);
    
            if min(newX,newY)<=0 || newX>gridX || newY>gridY 
                cell=nan;
            else
                onLS = true;
                cell=citNums(newX,newY);
            end
        end
        if ~isnan(cell) && citArray(cell)+othArray(cell)>0

            % SELECT THE TYPE OF HOST RANDOMLY
            hostType = rand<citArray(cell)/(othArray(cell)+citArray(cell));

            if hostType==1
                survivalProb = climArray(cell) * max(M(cell),(1-comInd));
                if rand<survivalProb
                    S=citArray-E-C-I-R;
        
                if fullInfestation==1 
                    infProb=(S(cell))/citArray(cell);
        
                else
                    infProb = ((I(oCell)+C(oCell))/citArray(oCell)) * (S(cell))/citArray(cell);
                end

                if N(cell)<carrCap(cell)
                    event=2;
                    A(cell) = 1;
                    if infProb>rand
                        event=3;
        
                        E(cell) = E(cell) + 1;
                           
                    end
                else
                     if infProb>rand
                        event=4;
        
                        E(cell) = E(cell) + 1;
                           
                    else
                        event=-1;
                    end
                end
                else
                    event=-1;
                end
            else 
                if hostType==0
                    survivalProb = climArray(cell) * max(M(cell),comInd);
                    if rand<survivalProb
                        Soth=othArray-Eoth-Coth-Ioth-Roth;
                    if fullInfestation==1 
                        infProb = (Soth(cell))/othArray(cell);
                    else
                        infProb = ((I(oCell)+C(oCell))/citArray(oCell)) * (Soth(cell))/othArray(cell);
                    end
                    if Noth(cell)<othCC(cell)
                        event=5;
                        Aoth(cell) = Aoth(cell) + 1;
                        if infProb>rand
                            event=6;
            
                            Eoth(cell) = Eoth(cell) + 1;
                               
                        end
                    else
                         if infProb>rand
                            event=7;
            
                            Eoth(cell) = Eoth(cell) + 1;
                               
                        else
                            event=-1;
                        end
                    end
                    else 
                        event=-1;
                    end
                end
            end
        else
            event=-2;
        end
   

    else % NEW PSYLLID INTRODUCTION
        if randLam < forceVCSum(end) + sum(ldd) + extVCSum(end)
    
            cell = find(randNum(2)*extVCSum(end)<extVCSum,1);
    
            if N(cell)< carrCap(cell)
                event=8;
                A(cell) = 1;    
            else
                if N(cell)>= carrCap(cell)
                    event=-3;
                end
            end

        else % PSYLLIDS START DISPERSING
            if randLam < lambdaVectorSum
    
            cell = find(randNum(2)*sum(estab)<cumsum(estab),1);
    
            if A(cell)>0 && N(cell)==0
                event=9;
                N(cell) = 1;    
                A(cell) = 0;    
            else
               event=-4;
            end
    
                
    
        else % DISPERSAL INFECTION EVENT
            if randLam < lambdaVectorSum + sum(forceInfect)
                event=10;
    
                cell = find(randNum(2)* sum(forceInfect)< cumsum(forceInfect),1);
                E(cell) = E(cell) + 1;
    
            
        else % EXTERNAL TREE INFECTION EVENT 
            if randLam < lambdaVectorSum + sum(forceInfect) + extInfCSum(end)
                cell = find(randNum(2)*extInfCSum(end)<extInfCSum,1);
                S=citArray-E-C-I-R;
                if S(cell)~=0
                    event=11;
    
                    E(cell) = E(cell) + 1;
                else
                    event=-3;
                end
            else 
                if randLam < lambdaVectorSum + sum(forceInfect) + extInfCSum(end) + sum(latentRate)
                    event=12;
                    cell = find(randNum(2)*sum(latentRate)<cumsum(latentRate),1);
                    E(cell)=E(cell)-1;
                    C(cell)=C(cell)+1;

                else 
                if randLam < lambda
                    event=13;
                    cell = find(randNum(2)*sum(symptomRate)<cumsum(symptomRate),1);
                    C(cell)=C(cell)-1;
                    I(cell)=I(cell)+1;        
                end
                end
            end
            end
            end
        end
    end
end