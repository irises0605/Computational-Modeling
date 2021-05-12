%{
    Iris Liu
    CS346 Optional Homework
    May 2020
    Model: RK4 simulation of a ball 
           being dropped from 400m
    To run: Type script name in command line
%}

clear

%% Time Step Setup

% dt: time step (s)
% t_start: start time (s)
% t_final: end time (s)
% N: iteration #
dt = .01;
t_start = 0;
t_final = 15;
N = (t_final - t_start)/dt;

%% Global Constant Values

% gravity: acceleration due to gravity (m/s^2)
gravity = -9.81;

% mass: mass of the ball (kg)
% weight: weight of the ball (N)
% radius: radius of the ball (m)
% proj_area: projected area of the ball (m^2)
mass = .5;
radius = .05;
weight = gravity*mass;
proj_area = 3.14159*(radius^2);

% v_initial: initial velocity (m/s)
% p_initial: initial position (m)
v_initial = 0;
p_initial = 400;

%% Differential Equations & Helper Functions

% calc_air_friction: calculate air friction
% calc_total_force: calculate total force
% calc_acceleration: calculate instant acceleration
% calc_speed: calculate instant speed

calc_air_friction = @(velocity) ...
               -0.65 * proj_area * velocity * abs(velocity);
calc_total_force = @(velocity) ...
               weight + calc_air_friction(velocity);
calc_acceleration = @(velocity) ...
               calc_total_force(velocity)/mass;
calc_speed = @(velocity) abs(velocity);

% dvdt: instant rate of change of velocity
% dsdt: instant rate of change of position
dvdt = @(time, velocity) calc_acceleration(velocity);
dsdt = @(time, velocity) velocity;

%% Initial Conditions

% time: initial time
% aceeleration: initial acceleration
% velocity: initial velocity
% position: initial velocity
time(1) = t_start;
acceleration(1) = calc_acceleration(v_initial);
velocity(1) = v_initial;
speed(1) = calc_speed(v_initial);
position(1) = p_initial;

%% Simulation

% RK4 simulation loop
% i: simulation number 
for i = 1:N
    
    time(i+1) = time(i) + dt;
    
    a1 = calc_acceleration(velocity(i));
    dv1 = dt * dvdt(time(i), velocity(i));
    ds1 = dt * dsdt(time(i), velocity(i));
    
    a2 = calc_acceleration(velocity(i) + 0.5*dv1);
    dv2 = dt * dvdt(time(i) + 0.5*dt, velocity(i) + 0.5*dv1);
    ds2 = dt * dsdt(time(i) + 0.5*dt, velocity(i) + 0.5*dv1);
    
    a3 = calc_acceleration(velocity(i) + 0.5*dv2);
    dv3 = dt * dvdt(time(i) + 0.5*dt, velocity(i) + 0.5*dv2);
    ds3 = dt * dsdt(time(i) + 0.5*dt, velocity(i) + 0.5*dv2);
    
    a4 = calc_acceleration(velocity(i) + dv3);
    dv4 = dt * dvdt(time(i) + dt,velocity(i) + dv3);
    ds4 = dt * dsdt(time(i) + dt,velocity(i) + dv3);
    
    acceleration(i+1) = (a1 + 2*a2 + 2*a3 + a4)/6;
    velocity(i+1) = velocity(i) + (dv1 + 2*dv2 + 2*dv3 + dv4)/6;
    speed(i+1) = calc_speed(velocity(i+1));
    position(i+1) = position(i) + (ds1 + 2*ds2 + 2*ds3 + ds4)/6;

end

%% Plot
figure;
plot(time,position)
grid on
title('Position vs. Time')
xlabel('Time (s)')
ylabel('Height (m)')

figure;
plot(time,speed)
grid on
title('Speed vs. Time')
xlabel('Time (s)')
ylabel('Speed (m/s)')