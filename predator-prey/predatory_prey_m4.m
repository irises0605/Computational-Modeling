%{
    Iris Liu
    May 2020
    Model: Predator & Prey Model for P, Y and H
           where P hunts Y and H hunts both P and Y
    To run: Type script name in command line
%}

clear

%% Time Step Setup

% dt: time step (month)
% t_start: start time (month)
% t_final: end time (month)
% N: iteration #
dt = 0.0001;
t_start = 0;
t_final = 12;
N = (t_final - t_start)/dt;

%% Global Constant Values

% Y_birth_fraction: birth fraction of tuna
% P_birth_fraction: birth fraction of sharks
Y_birth_fraction = 2;
P_birth_fraction = 0.01;

% Y_death_constant: death fraction of tuna
% P_death_constant: death fraction of sharks
Y_death_constant = 0.02;
P_death_constant = 1.06;

% fishing_rate
fishing_rate = 0.5;

%% Initial Conditions

% t: time
% Y_pop: population of prey (tuna)
% P_pop: population of predator (shark)
t(1) = 0;
Y_pop(1) = 100;
P_pop(1) = 15;

%% Helper Functions

Y_births = @(Y_population) ...
           Y_birth_fraction * Y_population;
P_births = @(Y_population, P_population) ...
           P_birth_fraction * Y_population * P_population;

Y_deaths = @(Y_population, P_population) ...
           Y_death_constant * P_population * Y_population;
P_deaths = @(P_population) ...
           P_death_constant * P_population;
       
Y_hunts = @(Y_population) ...
           fishing_rate * Y_population;
P_hunts = @(P_population) ...
           fishing_rate * P_population;

%% Simulation

% main simulation loop over assigned time period
% i: simulation number 
for i = 1:N-1
    
    t(i+1) = t(i) + dt;
    Y_pop(i+1) = Y_pop(i) + (Y_births(Y_pop(i))...
                 - Y_deaths(Y_pop(i), P_pop(i)) - Y_hunts(Y_pop(i)))*dt;
    P_pop(i+1) = P_pop(i) + (P_births(Y_pop(i), P_pop(i))...
                 - P_deaths(P_pop(i)) - P_hunts(P_pop(i)))*dt;
            
end

%% Plot
figure;
plot(t,Y_pop,t,P_pop,'r-')
title('Tuna & Shark Population vs. Time')
xlabel('Time (month)')
legend('Tuna','Shark')