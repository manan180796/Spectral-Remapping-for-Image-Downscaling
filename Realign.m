function [ Result ] = Realign( olderWaves,newWaves )
%REALIGN Summary of this function goes here
%   Detailed explanation goes here
waveNum = 0;
imageSize = size(newWaves);
imageSize = imageSize(1:2);
for i = 1:imageSize(1)
    for j = 1:imageSize(2)
        waves = olderWaves{i,j};
        number_of_waves = size(waves);
        number_of_waves = number_of_waves(1);
        o = olderWaves{i,j};
        if size(o,1)==0
            continue;
        end
        o = [o,zeros(size(o,1),1)];
        o(:,5)=1:size(o,1);
        o(:,5) = o(:,5)+waveNum;
        olderWaves{i,j} = o;
        
        o = newWaves{i,j};
        o = [o,zeros(size(o,1),1)];
        o(:,5)=1:size(o,1);
        o(:,5) = o(:,5)+waveNum;
        newWaves{i,j} = o;
        
        waveNum = waveNum + size(o,1);
    end
end
totalEquations = 0;
for i = 1:imageSize(1)
    for j = 1:imageSize(2)
        exterior = 0;
        if(i>1)
            temp = olderWaves{i-1,j};
            exterior = exterior + size(temp,1);
        end
        if(j>1)
            temp = olderWaves{i,j-1};
            exterior = exterior + size(temp,1);
        end
        c = olderWaves{i,j};
        current = (size(c,1)*(size(c,1)-1))/2;
        totalEquations = totalEquations+current;
        totalEquations = totalEquations+exterior*(size(c,1));
    end
end
values = zeros(2*totalEquations,3);

B = zeros(totalEquations,1);
Equation = 0;
for i = 1:imageSize(1)
    for j = 1:imageSize(2)
%         [i,j]
        oWaves = olderWaves{i,j};
        nWaves = newWaves{i,j};
        for k = 1:size(oWaves,1)

            currentOldWave = oWaves(k,:)';
            currentNewWave = nWaves(k,:)';
            if(i>1)
                ouWaves = olderWaves{i-1,j};
                nuWaves = newWaves{i-1,j};
                for l = 1:size(ouWaves,1)
                    tempOldWave = ouWaves(l,:)';
                    tempNewWave = nuWaves(l,:)';
                    [a,b] = equation(currentOldWave,tempOldWave,currentNewWave,tempNewWave,[i-0.5,j]);
                    if(a==0)
                        continue;
                    end
                    Equation = Equation+1;
                    values(2*Equation-1,:) = [Equation,tempOldWave(5),-a];
                    values(2*Equation,:) = [Equation,currentOldWave(5),a];
                    B(Equation) = b;
                    
                end
            end
            if(j>1)
                ouWaves = olderWaves{i,j-1};
                nuWaves = newWaves{i,j-1};
                for l = 1:size(ouWaves,1)
                    tempOldWave = ouWaves(l,:)';
                    tempNewWave = nuWaves(l,:)';
                    [a,b] = equation(currentOldWave,tempOldWave,currentNewWave,tempNewWave,[i,j-0.5]);
                    if(a==0)
                        continue;
                    end
                    Equation = Equation+1;
                    values(2*Equation-1,:) = [Equation,tempOldWave(5),-a];
                    values(2*Equation,:) = [Equation,currentOldWave(5),a];
                    B(Equation) = b;
                    
                end
            end
            for l = k+1:size(oWaves,1)
                tempOldWave = oWaves(l,:)';
                tempNewWave = nWaves(l,:)';
                [a,b] = equation(currentOldWave,tempOldWave,currentNewWave,tempNewWave,[i,j]);
                if(a==0)
                   continue;
                end
                Equation = Equation+1;
                values(2*Equation-1,:) = [Equation,tempOldWave(5),-a];
                values(2*Equation,:) = [Equation,currentOldWave(5),a];
                B(Equation) = b;
            end
        end
    end
end

i1  =find(values(:,1)>0);
i2 = values(:,2);
A = sparse(values(values(:,1)>0,1),values(values(:,1)>0,2),values(values(:,1)>0,3), Equation,waveNum);
B = B(1:Equation);
 x = lsqr(A,B,[],400);
for i = 1:imageSize(1)
    for j = 1:imageSize(2)
        waves = newWaves{i,j};
        if size(waves,1)==0
            continue;
        end
        waves(:,4) = x(waves(1,5):waves(end,5));
        newWaves{i,j} = waves;
    end
end
Result = newWaves;
end

