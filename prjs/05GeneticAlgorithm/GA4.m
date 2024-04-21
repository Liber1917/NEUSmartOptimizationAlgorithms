%%%%%%%%%%%%%%%?????? 0-1 ????%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%?????%%%%%%%%%%%%%%%%%%%%%% 
clear all; %??? 
close all; %?? 
clc; %?? 

NP = 50; % ????
L = 10; % ?????
Pc = 0.8; % ????
Pm = 0.05; % ????
G = 100; % ????
V = 300; % ????
C = [95, 75, 23, 73, 50, 22, 6, 57, 89, 98]; % ?????
W = [89, 59, 19, 43, 100, 72, 44, 16, 7, 64]; % ?????
afa = 5; % ????

f = randi([0, 1], NP, L); % ????????

%%%%%%%%%%%%%%%%%%%%%??????%%%%%%%%%%%%%%%%%%%%% 
for k = 1:G 
    %%%%%%%%%%%%%%%%%%?????%%%%%%%%%%%%%%%%% 
    for i = 1:NP 
        Fit(i) = func4(f(i,:), C, W, V, afa); 
    end 

    maxFit = max(Fit); % ???
    minFit = min(Fit); % ???
    rr = find(Fit == maxFit);
    fBest = f(rr(1, 1), :); % ??????
    Fit = (Fit - minFit) / (maxFit - minFit); % ???????

    %%%%%%%%%%%%%%??????????%%%%%%%%%%%%% 
    sum_Fit = sum(Fit);
    fitvalue = Fit ./ sum_Fit;
    fitvalue = cumsum(fitvalue);
    ms = sort(rand(NP, 1));
    fiti = 1;
    newi = 1;
    while newi <= NP 
        if (ms(newi)) < fitvalue(fiti) 
            nf(newi, :) = f(fiti, :);
            newi = newi + 1;
        else
            fiti = fiti + 1;
        end 
    end 

    %%%%%%%%%%%%%%%?????????%%%%%%%%%%%%% 
    for i = 1:2:NP 
        p = rand;
        if p < Pc 
            q = randi([0, 1], 1, L);
            for j = 1:L 
                if q(j) == 1
                    temp = nf(i + 1, j);
                    nf(i + 1, j) = nf(i, j);
                    nf(i, j) = temp;
                end 
            end 
        end 
    end 

    %%%%%%%%%%%%%?????????%%%%%%%%%%%%%%
    for m = 1:NP 
        for n = 1:L 
            r = rand(1, 1);
            if r < Pm 
                nf(m, n) = ~nf(m, n);
            end 
        end 
    end 

    f = nf;
    f(1, :) = fBest; % ???????????
    trace(k) = maxFit; % ???????
end 

fBest; % ????
figure 
plot(trace) 
xlabel('IterNum') 
ylabel('TargetValue') 
title('FEC') 

function result = func4(f, C, W, V, afa) 
fit = sum(f .* W); 
TotalSize = sum(f .* C); 
if TotalSize <= V 
    fit = fit; 
else 
    fit = fit - afa * (TotalSize - V); 
end 
result = fit; 
end
