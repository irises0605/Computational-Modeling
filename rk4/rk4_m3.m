%{
    Iris Liu
    CS346 Optional Homework
    May 2020
    Model: RK4 simulation of a ball 
           being thrown striaght up from
           the bridge
    To run: Type script name in command line
%}

clear

%% Time Step Setup

% dt: time step (mS)
% t_start: start time (mS)
% t_final: end time (mS)
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
calc_acceleration = @(time, velocity, position) ...
                    gravity + 0.01*(velocity + position)+0.3*time^2;
calc_speed = @(velocity) abs(velocity);

% dvdt: instant rate of change of velocity
% dsdt: instant rate of change of position
dvdt = @(time, velocity, position) ...
        calc_acceleration(time, velocity, position);
dsdt = @(time, velocity) velocity;

%% Initial Conditions

% time: initial time
% aceeleration: initial acceleration
% velocity: initial velocity
% position: initial velocity
time(1) = t_start;
acceleration(1) = calc_acceleration(t_start,v_initial,p_initial);
velocity(1) = v_initial;
position(1) = p_initial;

%% Simulation

% RK4 simulation loop
% i: simulation number 
for i = 1:N
    
    time(i+1) = time(i) + dt;
    
    a1 = calc_acceleration(time(i),velocity(i),position(i));
    dv1 = dt * dvdt(time(i),velocity(i),position(i));
    ds1 = dt * dsdt(time(i),velocity(i));
    
    a2 = calc_acceleration(time(i)+0.5*dt,velocity(i)+0.5*dv1,...
                           position(i)+0.5*ds1);
    dv2 = dt * dvdt(time(i)+0.5*dt,velocity(i)+0.5*dv1,position(i)+0.5*ds1);
    ds2 = dt * dsdt(time(i)+0.5*dt,velocity(i)+0.5*dv1);
    
    a3 = calc_acceleration(time(i)+0.5*dt,velocity(i)+0.5*dv2,...
                           position(i)+0.5*ds2);
    dv3 = dt * dvdt(time(i)+0.5*dt,velocity(i)+0.5*dv2,position(i)+0.5*ds2);
    ds3 = dt * dsdt(time(i)+0.5*dt,velocity(i)+0.5*dv2);
    
    a4 = calc_acceleration(time(i)+dt,velocity(i)+dv3,position(i)+ds3);
    dv4 = dt * dvdt(time(i)+dt,velocity(i)+dv3,position(i)+ds3);
    ds4 = dt * dsdt(time(i)+dt,velocity(i)+dv3);
    
    acceleration(i+1) = (a1 + 2*a2 + 2*a3 + a4)/6;
    velocity(i+1) = velocity(i) + (dv1 + 2*dv2 + 2*dv3 + dv4)/6;
    position(i+1) = position(i) + (ds1 + 2*ds2 + 2*ds3 + ds4)/6;

end

%% Plot
figure;
plot(time,position,'g-',time,velocity,'b-')
grid on
title('Position and Velocity vs. Time')
xlabel('Time (s)')
ylabel('Height (m) or Speed (m/s)')
legend('position','velocity')