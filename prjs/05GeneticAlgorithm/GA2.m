%%%%%%%%%%%%%%%%%%%%???????????%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%???%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear all; %?????? 
close all; %?? 
clc; %?? 
D = 10; %???? 
NP = 100; %????? 
Xs = 20; %?? 
Xx = -20; %?? 
G = 1000; %?????? 
f = zeros(D, NP); %??????? 
nf = zeros(D, NP); %?????? 
Pc = 0.8; %???? 
Pm = 0.1; %???? 
f = rand(D, NP) * (Xs - Xx) + Xx; %???????? 
%%%%%%%%%%%%%%%%%%%%%%????????%%%%%%%%%%%%%%%%%%%%%%% 
for np = 1:NP 
    FIT(np) = func2(f(:, np)); 
end 
[SortFIT, Index] = sort(FIT); 
Sortf = f(:, Index); 
%%%%%%%%%%%%%%%%%%%%%%%??????%%%%%%%%%%%%%%%%%%%%%%%%%% 
trace = zeros(1, G); % ??????????
for gen = 1:G 
    %%%%%%%%%%%%%%?????????????%%%%%%%%%%%%%%%% 
    Emper = Sortf(:, 1); %????? 
    nf = Sortf; 
    
    % ???????
    total_fit = sum(1./FIT);
    fit_ratio = (1./FIT) / total_fit;
    
    % ?????
    selected_indices = zeros(1, NP);
    for i = 1:NP
        r = rand();
        acc = 0;
        for j = 1:NP
            acc = acc + fit_ratio(j);
            if r <= acc
                selected_indices(i) = j;
                break;
            end
        end
    end
    
    % ????????????
    for i = 1:NP
        nf(:, i) = Sortf(:, selected_indices(i));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%????%%%%%%%%%%%%%%%%%%%%%%%%%  
    for m = 1:NP 
        for n = 1:D 
            r = rand(1, 1); 
            if r < Pm 
                nf(n, m) = rand(1, 1) * (Xs - Xx) + Xx; 
            end 
        end 
    end 
    
    %%%%%%%%%%%%%%%%%%%%%???????????%%%%%%%%%%%%%%%%%% 
    for np = 1:NP 
        NFIT(np) = func2(nf(:, np)); 
    end 
    [NSortFIT, Index] = sort(NFIT); 
    NSortf = nf(:, Index); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%?????%%%%%%%%%%%%%%%%%%%%%%%%%% 
    f1 = [Sortf, NSortf]; %??????? 
    FIT1 = [SortFIT, NSortFIT]; %???????????? 
    [SortFIT1, Index] = sort(FIT1); %???????? 
    Sortf1 = f1(:, Index); %???????? 
    SortFIT = SortFIT1(1:NP); %?? NP ????? 
    Sortf = Sortf1(:, 1:NP); %?? NP ??? 
    trace(gen) = SortFIT(1); %???????? 
end 
Bestf = Sortf(:, 1); %???? 
trace(end) %??? 

% ??????
BestFit = func2(Bestf);
fprintf('????: %f\n', BestFit);

% ?????????????
figure 
plot(trace) 
xlabel('IterNum') 
ylabel('TargetValue') 
title('FEC') 

function result = func2(x) 
    summ = sum(x.^2); 
    result = summ; 
end