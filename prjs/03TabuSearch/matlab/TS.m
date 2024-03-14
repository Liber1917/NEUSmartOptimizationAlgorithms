%%%%%%%%%%%%%%%%%%%%%????????TSP??%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%???%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;                        %??????
close all;                        %?? 
clc;                              %??
C = [1304 2312;3639 1315;4177 2244;3712 1399;3488 1535;3326 1556;...
    3238 1229;4196 1044;4312  790;4386  570;3007 1970;2562 1756;...
    2788 1491;2381 1676;1332  695;3715 1678;3918 2179;4061 2370;...
    3780 2212;3676 2578;4029 2838;4263 2931;3429 1908;3507 2376;...
    3394 2643;3439 3201;2935 3240;3140 3550;2545 2357;2778 2826;...
    2370 2975];                   %31???????
N=size(C,1);                      %TSP?????,?????
D=zeros(N);                       %??????????
%%%%%%%%%%%%%%%%%%%%%???????????????%%%%%%%%%%%%%%%%%%%%%
for i=1:N
    for j=1:N
        D(i,j)=((C(i,1)-C(j,1))^2+...
            (C(i,2)-C(j,2))^2)^0.5;
    end
end
Tabu=zeros(N);                      %??????????????????????
TabuL=round((N*(N-1)/2)^0.5);       %???????????????????
Ca=230;                             %?????
CaNum=zeros(Ca,N);                  %????????
%S0=randperm(N);                     %???????
% ???????????
S0 = greedyAlgorithm(C);

% ??????????
S0 = localOptimization(D, S0);

bestsofar=S0;                       %?????
BestL=Inf;                          %???????
figure(1);
p=1;
Gmax=1200;                          %??????   


%%%%%%%%%%%%%%%%%%%%%%%%%%%??????%%%%%%%%%%%%%%%%%%%%%%%%%%
while p<Gmax
    ALong(p)=func1(D,S0);            %??????
    %%%%%%%%%%%%%%%%%%%%%%%%%%%????%%%%%%%%%%%%%%%%%%%%%%%%%%
    i=1;
    A=zeros(Ca,2);                   %????????????????????
    %%%%%%%%%%%%%%%%%????????????%%%%%%%%%%%%%%%%%%%%%
    while i<=Ca
        M=N*rand(1,2); %%%????? 1?2? 0.9213 0.5334
        M=ceil(M); %%%% w=ceil(z)?????z????????w????????????   
        if M(1)~=M(2) 
            A(i,1)=max(M(1),M(2));
            A(i,2)=min(M(1),M(2));
            if i==1
                isa=0;
            else
                for j=1:i-1
                    if A(i,1)==A(j,1) && A(i,2)==A(j,2)
                        isa=1;
                        break;
                    else
                        isa=0;
                    end
                end
            end
            if ~isa
                i=i+1;
            else
            end
        else
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%?????%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%???BestCaNum??????%%%%%%%%%%%%%%%%%%%
    BestCaNum=Ca/2; %%????
    BestCa=Inf*ones(BestCaNum,4);
    F=zeros(1,Ca);
    for i=1:Ca
        CaNum(i,:)=S0;%%??????????
        CaNum(i,[A(i,2),A(i,1)])=S0([A(i,1),A(i,2)]);%%??A(i, 1)?A(i, 2)????
        F(i)=func1(D,CaNum(i,:));%%?????
        if i<=BestCaNum %?????
            BestCa(i,1)=i; %% ?1?????   
            BestCa(i,2)=F(i);%%  ?2???????(????)   
            BestCa(i,3)=S0(A(i,1));%% ?3?????????
            BestCa(i,4)=S0(A(i,2));%% ?4?????????
        else
            for j=1:BestCaNum % ?????
                if F(i)<BestCa(j,2)
                    BestCa(j,2)=F(i);
                    BestCa(j,1)=i;
                    BestCa(j,3)=S0(A(i,1));
                    BestCa(j,4)=S0(A(i,2));
                    break;
                end
            end
        end
    end
    % ????????????????
    [JL,Index]=sort(BestCa(:,2));
    SBest=BestCa(Index,:);
    BestCa=SBest; % ???100????????
    %%%%%%%%%%%%%%%%%%%%%%%%????%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if BestCa(1,2)<BestL % ??????????????????????????????
        % ?????????1???????????
        % ??????????1??????
        BestL=BestCa(1,2); % BestL????????
        S0=CaNum(BestCa(1,1),:);% ??????
        bestsofar=S0;
        for m=1:N
            for n=1:N
                if Tabu(m,n)~=0
                    Tabu(m,n)=Tabu(m,n)-1;    
                    %?????
                end
            end
        end
        Tabu(BestCa(1,3),BestCa(1,4))=TabuL; % ?????????????????
        %?????
    else % ???????? ???????????????
        for  i=1:BestCaNum % ?????
            if  Tabu(BestCa(i,3),BestCa(i,4))==0% BestCa??????????????????????????????0
                S0=CaNum(BestCa(i,1),:);% ????????????????
                for m=1:N
                    for n=1:N
                        if Tabu(m,n)~=0
                            Tabu(m,n)=Tabu(m,n)-1;
                            %?????
                        end
                    end
                end
                Tabu(BestCa(i,3),BestCa(i,4))=TabuL;% ????????
                %?????
                break;%  ????for????????????????????
            end
        end
    end
    ArrBestL(p)=BestL;
    p=p+1;   
    for i=1:N-1
        plot([C(bestsofar(i),1),C(bestsofar(i+1),1)],...
            [C(bestsofar(i),2),C(bestsofar(i+1),2)],'bo-');
        hold on;
    end
    plot([C(bestsofar(N),1),C(bestsofar(1),1)],...
        [C(bestsofar(N),2),C(bestsofar(1),2)],'ro-');
    title(['Optimized minimum distance:',num2str(BestL)]);
    hold off;
    pause(0.005);
end
BestShortcut=bestsofar;            %????
theMinDistance=BestL;              %??????
figure(2);
plot(ArrBestL); 
xlabel('Iteration times')
ylabel('Objective function value')
title('Fitness evolution curve')



function total_distance = func1(distance_matrix, solution)
    N = length(solution);
    total_distance = 0;
    for i = 1:N-1
        total_distance = total_distance + distance_matrix(solution(i), solution(i+1));
    end
    % Add distance from last city to first city to complete the cycle
    total_distance = total_distance + distance_matrix(solution(N), solution(1));
end


function initialSolution = greedyAlgorithm(C)
    N = size(C, 1);
    initialSolution = zeros(1, N);
    visited = zeros(1, N);
    
    % ??????????
    currentCity = randi(N);
    initialSolution(1) = currentCity;
    visited(currentCity) = 1;
    
    % ?????????????
    for i = 2:N
        minDistance = Inf;
        nextCity = -1;
        for j = 1:N
            if ~visited(j)
                distance = norm(C(currentCity, :) - C(j, :));
                if distance < minDistance
                    minDistance = distance;
                    nextCity = j;
                end
            end
        end
        initialSolution(i) = nextCity;
        visited(nextCity) = 1;
        currentCity = nextCity;
    end
end

function improvedSolution = localOptimization(D, solution)
    improvedSolution = solution;
    N = length(solution);
    improved = true;
    
    while improved
        improved = false;
        for i = 1:N-1
            for j = i+1:N
                newSolution = swap(solution, i, j);
                if calculateDistance(D, newSolution) < calculateDistance(D, improvedSolution)
                    improvedSolution = newSolution;
                    improved = true;
                end
            end
        end
    end
end

function newSolution = swap(solution, i, j)
    newSolution = solution;
    temp = newSolution(i);
    newSolution(i) = newSolution(j);
    newSolution(j) = temp;
end

function distance = calculateDistance(D, solution)
    distance = 0;
    N = length(solution);
    for i = 1:N-1
        distance = distance + D(solution(i), solution(i+1));
    end
    % Add distance from last city to first city to complete the cycle
    distance = distance + D(solution(N), solution(1));
end

