%%%%%%%%%%%%?????? TSP ??%%%%%%%%%%%%%%%%%%%%%%% 
clear all; %?????? 
close all; %?? 
clc; %?? 

% 31 ???????
C = [1304 2312; 3639 1315; 4177 2244; 3712 1399; 3488 1535; 3326 1556; ...
    3238 1229; 4196 1044; 4312 790; 4386 570; 3007 1970; 2562 1756; ...
    2788 1491; 2381 1676; 1332 695; 3715 1678; 3918 2179; 4061 2370; ...
    3780 2212; 3676 2578; 4029 2838; 4263 2931; 3429 1908; 3507 2376; ...
    3394 2643; 3439 3201; 2935 3240; 3140 3550; 2545 2357; 2778 2826; ...
    2370 2975];

% TSP ?????,?????
N = size(C, 1);

% ????????????
D = zeros(N);

% ?????????????
for i = 1:N 
    for j = 1:N 
        D(i, j) = sqrt((C(i, 1) - C(j, 1))^2 + (C(i, 2) - C(j, 2))^2); 
    end 
end 

NP = 400; % ????
G = 2000; % ??????

f = zeros(NP, N); % ??????
F = []; % ????????

for i = 1:NP 
    f(i, :) = randperm(N); % ????????
end 

R = f(1, :); % ??????
len = zeros(NP, 1); % ??????
fitness = zeros(NP, 1); % ????????
gen = 0;

% ??????
while gen < G 
    % ??????
    for i = 1:NP 
        len(i, 1) = D(f(i, N), f(i, 1)); 
        for j = 1:(N - 1) 
            len(i, 1) = len(i, 1) + D(f(i, j), f(i, j + 1)); 
        end 
    end 
    maxlen = max(len); % ????
    minlen = min(len); % ????

    % ??????
    rr = find(len == minlen); 
    R = f(rr(1, 1), :); 

    % ????????
    for i = 1:length(len) 
        fitness(i, 1) = 1 - ((len(i, 1) - minlen) / (maxlen - minlen + 0.001)); 
    end 

    % ????
    nn = 0; 
    for i = 1:NP 
        if fitness(i, 1) >= rand 
            nn = nn + 1; 
            F(nn, :) = f(i, :); 
        end 
    end 

    [aa, bb] = size(F); 
    while aa < NP 
        nnper = randperm(nn); 
        A = F(nnper(1), :); 
        B = F(nnper(2), :); 

        % ????
        for i = 1:2:NP
            A = F(i, :);
            B = F(i+1, :);

            % ??????????
            method = randi([1, 3]); % ?????????1-PMX, 2-OX, 3-CX

            % ???????????????
            switch method
                case 1 % PMX
                    [A, B] = PMX(A, B);
                case 2 % OX
                    [A, B] = OX(A, B);
                case 3 % CX
                    [A, B] = CX(A, B);
            end

            % ??????????????
            F(i, :) = A;
            F(i+1, :) = B;
        end

        % ????
        p1 = floor(1 + N * rand()); 
        p2 = floor(1 + N * rand()); 
        while p1 == p2 
            p1 = floor(1 + N * rand());  
            p2 = floor(1 + N * rand()); 
        end 

        tmp = A(p1); 
        A(p1) = A(p2); 
        A(p2) = tmp; 
        tmp = B(p1); 
        B(p1) = B(p2); 
        B(p2) = tmp; 

        F = [F; A; B]; 
        [aa, bb] = size(F); 
    end 

    if aa > NP 
        F = F(1:NP, :); % ??????? NP
    end 

    f = F; % ????
    f(1, :) = R; % ????????
    clear F; 

    gen = gen + 1; 
    Rlength(gen) = minlen; 
end 

figure 
for i = 1:N - 1 
    plot([C(R(i), 1), C(R(i + 1), 1)], [C(R(i), 2), C(R(i + 1), 2)], 'bo-'); 
    hold on; 
end 
plot([C(R(N), 1), C(R(1), 1)], [C(R(N), 2), C(R(1), 2)], 'ro-'); 
title(['OMD:', num2str(minlen)]); 

figure 
plot(Rlength) 
xlabel('IterNum') 
ylabel('TargetValue') 
title('FEC') 


function [A, B] = PMX(A, B)
    N = length(A);
    % ??????????
    p1 = randi([1, N]);
    p2 = randi([1, N]);
    while p1 == p2
        p2 = randi([1, N]);
    end
    
    % ?? p1 < p2
    if p1 > p2
        temp = p1;
        p1 = p2;
        p2 = temp;
    end
    
    % PMX??
    tempA = A;
    tempB = B;
    A(p1:p2) = tempB(p1:p2);
    B(p1:p2) = tempA(p1:p2);
    
    for i = 1:N
        if i < p1 || i > p2
            % ????
            while ismember(tempA(i), A(p1:p2))
                idx = find(tempA == tempB(i));
                tempA(i) = tempA(idx);
            end
            while ismember(tempB(i), B(p1:p2))
                idx = find(tempB == tempA(i));
                tempB(i) = tempB(idx);
            end
        end
    end
end

function [A, B] = OX(A, B)
    N = length(A);
    % ??????????
    p1 = randi([1, N]);
    p2 = randi([1, N]);
    while p1 == p2
        p2 = randi([1, N]);
    end
    
    % ?? p1 < p2
    if p1 > p2
        temp = p1;
        p1 = p2;
        p2 = temp;
    end
    
    % ??????
    middleA = A(p1:p2);
    middleB = B(p1:p2);
    
    % ???????
    A(p1:p2) = middleB;
    B(p1:p2) = middleA;
    
    % ??????
    restA = setdiff(B, middleA);
    restB = setdiff(A, middleB);
    lenRestA = length(restA);
    lenRestB = length(restB);
    A(1:p1-1) = [restA(1:min(p1-1, lenRestA)), A(p1:p2-1)];
    A(p2:N) = restA(min(p1, lenRestA)+1:end);
    B(1:p1-1) = [restB(1:min(p1-1, lenRestB)), B(p1:p2-1)];
    B(p2:N) = restB(min(p1, lenRestB)+1:end);
end


function [A, B] = CX(A, B)
    N = length(A);
    visited = zeros(1, N);
    cycle = zeros(1, N);
    idx = 1;
    
    while ~all(visited)
        % ?????????
        while visited(idx)
            idx = idx + 1;
            if idx > N
                idx = 1;
            end
        end
        
        % ??????
        start = idx;
        current = start;
        while ~visited(current)
            cycle(current) = B(current);
            visited(current) = 1;
            current = find(A == B(current));
        end
        
        % ????????
        temp = A;
        A = B;
        B = temp;
    end
    
    % ?????????????
    A(visited) = cycle(visited);
    B(visited) = cycle(visited);
end
