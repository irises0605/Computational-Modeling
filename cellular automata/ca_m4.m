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
totalsteps = 200;

%% Global Constant Values
% n: grid size n x n
n = 20;

% probSusceptible: probability the individual is initially susceptible
% probInfectious: probability an individual that is infectious initially
% probImmune: probability the individual is initially immune
probSusceptible = 0.70;
probInfectious = 0.10;
probImmune = 1 - probSusceptible - probInfectious;

% probBeSusceptible: probability the individual will become susceptible
probBeSusceptible = 0.9;

% max_infect: maximum infection rate in the flu season
% period: how long it would last each year
% probLongInfect: probility an individual will have infection period longer
%                 than 2
max_infect = 0.8;
period = 60;
probLongInfect = @(day) max_infect*cos(-0.5*pi+day*pi/period)...
                    *(mod(day,2*period) <= period);

%% Phases

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
Type = 8;

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
grid(:,[1 n]) = 8;
grid([1 n],:) = 8;

%% Simulation
stats = zeros(Type,totalsteps);
S_stat = zeros(1,totalsteps);
for days = 1:timestep:totalsteps
    display(grid);
    nextgrid = grid;
    for x = 2:n-1
        for y = 2:n-1
            
            % If susceptible
            if grid(x,y) == S
                neighnors = [grid(x-1,y), grid(x,y-1), ...
                             grid(x+1,y), grid(x,y+1)];
                infected_neigh = sum(neighnors == IN_1) + ...
                                 sum(neighnors == IN_2);
                
                % Probability being infected based on the percentage of
                % infected neighbors
                if rand < infected_neigh/length(neighnors)
                    nextgrid(x,y) = IN_1;
                end
            end
            
            % If infected/immune
            if grid(x,y) > S
                
                % If immune for 5 days, there's a chance to be susceptile
                % or the individual would have an immune period longer than
                % 5 days
                if grid(x,y) == IM_5 && rand < (1-probBeSusceptible)
                    nextgrid(x,y) = IM_5;
                else
                    nextgrid(x,y) = mod(grid(x,y)+1, Type);
                end
                
                % Probability having a longer infection period
                if grid(x,y) == IN_2 && rand < probLongInfect(days)
                    nextgrid(x,y) = IN_2;
                else
                    nextgrid(x,y) = mod(grid(x,y)+1, Type);
                end
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

clist = [0 1 0; 1 0 0; 1 0.5 0;
         0 0 1; 0 0.3 1; 0 0.5 1; 0 0.7 1; 0 0.9 1
         1 1 1];   
colormap(clist);
pause(0.01);

end       