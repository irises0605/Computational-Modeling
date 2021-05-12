%{
    Iris Liu
    May 2020
    Model: Predator & Prey Model for P, Y and H
           where P hunts Y and H hunts P
    To run: Type script name in command line
%}

clear

%% Time Step Setup

% dt: time step (month)
% t_start: start time (month)
% t_final: end time (month)
% N: iteration #
dt = 0.001;
t_start = 0;
t_final = 12*12;
N = (t_final - t_start)/dt;

%% Global Constant Values

% Y_birth_fraction: birth fraction of tuna
% P_birth_fraction: birth fraction of sharks
% H_birth_fraction: birth fraction of human
Y_birth_fraction = 2;
P_birth_fraction = 0.01;
H_birth_fraction = 0.002;

% Y_death_constant: death constant of tuna
% P_death_constant: death constant of sharks
% H_death_constant: death constant of human
Y_death_constant = 0.02;
P_death_constant = 0.0106;
H_death_constant = 0.199;

%% Initial Conditions

% t: time
% Y_pop: population of tuna
% P_pop: population of shark
% P_pop: population of human
t(1) = 0;
Y_pop(1) = 500;
P_pop(1) = 150;
H_pop(1) = 20;

%% Helper Functions

Y_births = @(Y_population) ...
           Y_birth_fraction * Y_population;
P_births = @(Y_population, P_population) ...
           P_birth_fraction * (Y_population + 1) * P_population;
H_births = @(Y_population, P_population, H_population) ...
           H_birth_fraction * (Y_population + P_population) * H_population;

Y_deaths = @(Y_population, P_population, H_population) ...
           Y_death_constant * (P_population + H_population) * Y_population;
P_deaths = @(P_population, H_population) ...
           P_death_constant * (H_population + 1)* P_population;
H_deaths = @(H_population) ...
           H_death_constant * H_population;
       
%% Simulation

% main simulation loop over assigned time period
% i: simulation number 
for i = 1:N-1
    
    t(i+1) = t(i) + dt;
    Y_pop(i+1) = Y_pop(i) + (Y_births(Y_pop(i))...
                 - Y_deaths(Y_pop(i), P_pop(i), H_pop(i)))*dt;
    P_pop(i+1) = P_pop(i) + (P_births(Y_pop(i), P_pop(i))...
                 - P_deaths(P_pop(i), H_pop(i)))*dt;
    H_pop(i+1) = H_pop(i) + (H_births(Y_pop(i), P_pop(i), H_pop(i))...
                 - H_deaths(H_pop(i)))*dt;
end

%% Plot
figure;
plot(t,Y_pop,t,P_pop,'r-',t,H_pop,'g-')
title('Tuna & Shark Population vs. Time')
xlabel('Time (month)')
legend('Tuna','Shark','Human')

% figure;
% plot(t, Y_pop)