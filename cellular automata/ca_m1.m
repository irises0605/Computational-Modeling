%{
    Iris Liu
    June 2020
    Project: SIR Model
    To run: Type script name in command line
%}

clear;
close all;

%% Timestep
timestep = 1;
totalsteps = 100;

%% Global Constant Values
% n: grid size n x n
n = 20;

% probSusceptible: probability the individual is initially susceptible
% probInfectious: probability an individual that is infectious initially
% probImmune: probability the individual is initially immune
% probRandom: probability the individual being random infected
probSusceptible = 0.50;
probInfectious = 0.50;
probImmune = 1 - probSusceptible - probInfectious;
probRandom = 0.5;

% vac_thres: threshold for vaccine implementation
% protect_rate: protection rate of the vaccine
% close_thres: threshold for vaccine implementation
vac_thres = 0.2;
protect_rate = 0.5;
close_thres = 0.1;
date = 60;

%% Values
% Type: 8 types of phases
Type = 8;

% 0: Susceptible
% 1: Infectious Day 1
% 2: Infectious Day 2
% 3: Immune Day 1
% 4: Immune Day 2
% 5: Immune Day 3
% 6: Immune Day 4
% 7: Immune Day 5
S = 0;
IN_1 = 1;
IN_2 = 2;
IM_1 = 3;
IM_2 = 4;
IM_3 = 5;
IM_4 = 6;
IM_5 = 7;

% WALL: edges 
WALL = 8;

%% Initial Grid

% initial n x n matrix
grid = zeros(n);
for x = 1:n
    for y = 1:n
        if rand <= probInfectious/2
            grid(x,y) = IN_1;
        elseif probInfectious/2 <= rand < probInfectious
            grid(x,y) = IN_2;
        elseif probInfectious <= rand < probInfectious + probImmune/5
            grid(x,y) = IM_1;
        elseif probInfectious + probImmune/5 <= rand < probInfectious + 2*probImmune/5
            grid(x,y) = IM_2;
        elseif probInfectious + 2*probImmune/5 <= rand < probInfectious + 3*probImmune/5
            grid(x,y) = IM_3;
        elseif probInfectious + 3*probImmune/5 <= rand < probInfectious + 4*probImmune/5
            grid(x,y) = IM_4;
        elseif probInfectious + 4*probImmune/5 <= rand < probInfectious + probImmune
            grid(x,y) = IM_5;
        end
    end
end

% Edge values of the grid to white space
grid(:,[1 n]) = WALL;
grid([1 n],:) = WALL;

%% Simulation

% stats: infected/immune population stats
% S_stat: susceptible population stat
stats = zeros(Type,totalsteps);
S_stat = zeros(1,totalsteps);
vac = false;

for days = 1:timestep:totalsteps
    display(grid);
    nextgrid = grid;
    
    % If the infection reaches the threshold for vaccine implementation
    if (sum(sum(grid(2:n-1,2:n-1)==IN_1)) + ...
        sum(sum(grid(2:n-1,2:n-1)==IN_2)))/n^2 > vac_thres && days > date
        vac = true;
    
    % If the infection rate is low enough to discard vaccine
    elseif (sum(sum(grid(2:n-1,2:n-1)==IN_1)) + ...
        sum(sum(grid(2:n-1,2:n-1)==IN_2)))/n^2 < close_thres
        vac = false;
    end
    
    for x = 2:n-1
        for y = 2:n-1
            
            % If susceptile the individual will be infected when 
            if grid(x,y) == S
                
                % Infected by neighbors
                neighnors = [grid(x-1,y), grid(x,y-1), ...
                             grid(x+1,y), grid(x,y+1)];
                if sum(neighnors == IN_1) +  sum(neighnors == IN_2) > 0
                    nextgrid(x,y) = IN_1;
                end
                
                % Randomly Infected
                if rand < probRandom
                    nextgrid(x,y) = IN_1;
                end
                
                % Implementation of vaccine
                if vac == true && rand < protect_rate
                    nextgrid(x,y) = IM_1;
                end
                
                
            end
            
            % If infected or immune
            if grid(x,y) > S && grid(x,y) < WALL
                
                % Individual must have 2-day infection period and 
                % 5-day immune period
                nextgrid(x,y) = mod(grid(x,y)+1, Type);
            end
        end
    end
    
    % Record the number of people in each phase
    S_stat(:,days) = sum(sum(grid(2:n-1,2:n-1)==S));
    stats(IN_1,days) = sum(sum(grid(2:n-1,2:n-1)==IN_1));
    stats(IN_2,days) = sum(sum(grid(2:n-1,2:n-1)==IN_2));
    stats(IM_1,days) = sum(sum(grid(2:n-1,2:n-1)==IM_1));
    stats(IM_2,days) = sum(sum(grid(2:n-1,2:n-1)==IM_2));
    stats(IM_3,days) = sum(sum(grid(2:n-1,2:n-1)==IM_3));
    stats(IM_4,days) = sum(sum(grid(2:n-1,2:n-1)==IM_4));
    stats(IM_5,days) = sum(sum(grid(2:n-1,2:n-1)==IM_5));
    
    % Renew the grid
    grid=nextgrid;
end

%% Plot & Stats
figure;
IN_total = stats(IN_1,:)+stats(IN_2,:);
IM_total = stats(IM_1,:)+stats(IM_2,:)+stats(IM_3,:)+stats(IM_4,:)+stats(IM_5,:);

plot(1:totalsteps, S_stat, 'r-', 1:totalsteps, IN_total, 'g-', ...
     1:totalsteps, IM_total, 'b-');
legend('Total Susceptible', 'Total Infected', 'Total Immune');

title('Population Stats vs. Time')
xlabel('Time (days)')
ylabel('Number of People')

%% Display
function display(x)

imagesc(x);
axis off;

% Susceptible: Full green
% Infectious Day 1: Full red
% Infectious Day 2: Paler red (orange)
% Immune Day 1: Blue
% Immune Day 2-5: Paler Blue
% Wall: Red
clist = [0 1 0; 1 0 0; 1 0.5 0;
         0 0 1; 0 0.3 1; 0 0.5 1; 0 0.7 1; 0 0.9 1
         1 1 1];   
colormap(clist);
pause(0.001);

end       