%{
    Iris Liu
    May 2020
    Model: Euler's simulation of a ball 
           being thrown striaght up from
           the bridge
    To run: Type script name in command line
%}

clear

%% Time Step Setup

% dt: time step (s)
% t_start: start time (s)
% t_final: end time (s)
% N: iteration #
dt = .25;
t_start = 0;
t_final = 4;
N = (t_final - t_start)/dt;

%% Global Constant Values

% gravity: acceleration due to gravity (m/s^2)
gravity = -9.81;

% v_initial: initial velocity (m/s)
% p_initial: initial position (m)
v_initial = 15;
p_initial = 11;

%% Differential Equations & Helper Functions

% calc_acceleration: calculate instant acceleration
% calc_speed: calculate instant speed
calc_acceleration = @(velocity) gravity;
calc_speed = @(velocity) abs(velocity);

% dvdt: instant rate of change of velocity
% dsdt: instant rate of change of position
dvdt = @(velocity) calc_acceleration(velocity);
dsdt = @(velocity) velocity;

%% Initial Conditions

% time: initial time
% aceeleration: initial acceleration
% velocity: initial velocity
% position: initial velocity
time(1) = t_start;
acceleration(1) = calc_acceleration(v_initial);
velocity(1) = v_initial;
position(1) = p_initial;

%% Simulation

% Euler simulation loop
% i: simulation number 
for i = 1:N
    
    time(i+1) = time(i) + dt;
    acceleration(i+1) = calc_acceleration(velocity(i));
    velocity(i+1) = velocity(i) + dvdt(velocity(i))*dt;
    position(i+1) = position(i) + dsdt(velocity(i))*dt;

end

%% Plot
figure;
plot(time,position,'g-',time,velocity,'b-')
grid on
title('Position and Velocity vs. Time')
xlabel('Time (s)')
ylabel('Height (m) or Speed (m/s)')
legend('position','velocity')