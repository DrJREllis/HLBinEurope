function [Fi, Kd, Ki, coord] = localArea(localDispPar,sameCellInf,maxDist,gridX,gridY,citNums)
% Function to generate the local dispersal kernel for the vector and pathogen

Fi1 = cell(gridX,gridY);

nums=reshape(1:gridX*gridY,gridX,gridY);
X = repmat((1:gridX)',1,gridY);
Y = repmat(1:gridY,gridX,1);

if gridX<=maxDist && gridY<=maxDist
for j=[1:gridY]
    for i=[1:gridX]
        dist = sqrt((i-X).^2+(j-Y).^2);
        Fi1{i,j} = nums(dist<=maxDist);
        Kd1{i,j} = exp(-dist(Fi1{i,j})/localDispPar) ;%./ (pi*localDispPar^2);
        Ki1{i,j} = Kd1{i,j};
        Ki1{i,j}(dist(Fi1{i,j})==0) = 0;
    end
end
normd = 19.6516;
normi = 18.6516;
else

for j=[1:maxDist+1, max(1,gridY-maxDist):gridY]
    for i=[1:maxDist+1, max(1,gridX-maxDist):gridX]
        dist = sqrt((i-X).^2+(j-Y).^2);
        Fi1{i,j} = nums(dist<=maxDist);
        Kd1{i,j} = exp(-dist(Fi1{i,j})/localDispPar) ;%./ (pi*localDispPar^2);
        Ki1{i,j} = Kd1{i,j};
        Ki1{i,j}(dist(Fi1{i,j})==0) = 0;
    end
    for i=maxDist + 2: gridX-maxDist-1

        Fi1{i,j} = Fi1{i-1,j} + 1;
        Kd1{i,j} = Kd1{i-1,j};
        Ki1{i,j} = Ki1{i-1,j};
    end
end

for j=maxDist+2:gridY-maxDist-1
    for i=1:gridX

        Fi1{i,j} = Fi1{i,j-1}+gridX;
        Kd1{i,j} = Kd1{i,j-1};
        Ki1{i,j} = Ki1{i,j-1};
    end
end
normd = sum(Kd1{maxDist+1,maxDist+1});
normi = sum(Ki1{maxDist+1,maxDist+1});
end

k=1;
    for j=1:gridY
for i=1:gridX
        if ~isnan(citNums(i,j))
            test(k)=citNums(i,j)';
            Fi{k} = citNums(Fi1{i,j});
            Kd{k} =  Kd1{i,j}(~isnan(Fi{k}))/ normd;

            Ki{k} = (1-sameCellInf) * Ki1{i,j}(~isnan(Fi{k}))/ normi;
            Fi{k} = Fi{k}(~isnan(Fi{k}));
            Ki{k}(Fi{k}==citNums(i,j)) = sameCellInf;

            k=k+1;
        end
    end
end

xMat=repmat((1:gridX)',1,gridY);
yMat = repmat(1:gridY,gridX,1);
coord = [xMat(citNums>0), yMat(citNums>0)];
