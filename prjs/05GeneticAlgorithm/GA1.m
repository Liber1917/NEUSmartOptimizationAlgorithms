%%%%%%%%%%%%%%%%%%%%???????????%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%?????%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear all; %?????? 
close all; %?? 
clc; %?? 
NP = 500; %???? 
L = 20; %??????? 
Pc = 0.4; %??? 
Pm = 0.2; %??? 
G = 50; %?????? 
Xs = 10; %?? 
Xx = 0; %?? 
f = randi([0,1], NP, L); %???????? 
%%%%%%%%%%%%%%%%%%%%%%%%%??????%%%%%%%%%%%%%%%%%%%%%%%% 
for k = 1:G 
    %%%%%%%%%%%%????????????????%%%%%%%%%%%%%% 
    for i = 1:NP 
        U = f(i,:); 
        m = 0; 
        for j = 1:L 
            m = U(j)*2^(j-1) + m; 
        end 
        x(i) = Xx + m*(Xs-Xx)/(2^L-1); 
        Fit(i) = func1(x(i)); 
    end 
    maxFit = max(Fit); %??? 
    minFit = min(Fit); %??? 
    rr = find(Fit == maxFit); 
    fBest = f(rr(1,1),:); %?????? 
    xBest = x(rr(1,1)); 
    Fit = (Fit - minFit) / (maxFit - minFit); %??????? 
    %%%%%%%%%%%%%%%%%%??????????%%%%%%%%%%%%%%%%%%% 
    sum_Fit = sum(Fit); 
    fitvalue = Fit / sum_Fit; 
    fitvalue = cumsum(fitvalue); 
    ms = sort(rand(NP,1)); 
    fiti = 1;  
    newi = 1; 
    while newi <= NP 
        if (ms(newi)) < fitvalue(fiti) 
            nf(newi,:) = f(fiti,:); 
            newi = newi + 1; 
        else 
            fiti = fiti + 1; 
        end 
    end 
    %%%%%%%%%%%%%%%%%%%%%%?????????%%%%%%%%%%%%%%%%%% 
    for i = 1:2:NP 
        p = rand; 
        if p < Pc 
            q = randi([0,1], 1, L); 
            for j = 1:L 
                if q(j) == 1 
                    temp = nf(i+1,j); 
                    nf(i+1,j) = nf(i,j); 
                    nf(i,j) = temp; 
                end 
            end 
        end 
    end 
    %%%%%%%%%%%%%%%%%%%?????????%%%%%%%%%%%%%%%%%%%%%%% 
    i = 1; 
    while i <= round(NP * Pm) 
        h = randi([1,NP], 1, 1); %?????????????? 
        for j = 1:round(L * Pm) 
            g = randi([1,L], 1, 1); %?????????? 
            nf(h,g) = ~nf(h,g); 
        end 
        i = i + 1; 
    end 
    f = nf; 
    f(1,:) = fBest; %??????????? 
    trace(k) = maxFit; %??????? 
end 
maxTargetValue = max(trace);
fprintf('MAX_TargetValue: %f\n', maxTargetValue);

xBest; %???? 
figure 
plot(trace) 
xlabel('IterNum') 
ylabel('TargetValue') 
title('FEC')


function result = func1(x) 
    fit = x + 10*sin(5*x) + 7*cos(4*x); 
    result = fit; 
end
