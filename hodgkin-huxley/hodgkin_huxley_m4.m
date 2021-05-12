%{
    Iris Liu
    May 2020
    Model: Hodgkin-Huxley Model with
           voltage-sensitive sodium and potassium channels, 
           leakage channels, and concentration-sensitive Na-K pump
           implementing RK4 simulations
    To run: Type script name in command line
%}

clear

%% Time Step Setup

% dt: time step (mS)
% t_final: start time (mS)
% t_final: end time (mS)
% N: iteration #
dt = .001;
t_start = 0;
t_final = 3;
N = (t_final - t_start)/dt;

%% Global Constant Values

% maximum conductance (mS/cm^2)
% g_K: maximum potassuium conductance
% g_Na: maximum sodium conductance
% g_L: maximum leakage conductance
g_K = 36;
g_Na = 120;
g_L = 0.3;

% equilibirum potentials (mV)
% V_rest: membrane potential at rest
% V_K: potassium equilibrium potential
% V_Na: sodium equilibrium potential
% V_L: leakage equilibrium potential
V_rest = -65;
V_K = -77;
V_Na = 50;
V_L = -54.4;

% C_m: membrane capacitance (uF/cm^2)
C_m = 0.1;

% thresholds potentials for channel flows (mV)
% o_thres: opening for sodium
% c_thres: closing for sodium & opening for potassium
open_thres = -55;
close_thres = 50;

%% Initial Conditions

% V: initial action potential
% n: initial active potassium channels activation probability
% m: initial active sodium channels activation probability 
% h: initial active sodium channels inactivation probability
V(1) = V_rest;
n(1) = 0.317;
m(1) = 0.05;
h(1) = 0.6;
t(1) = 0;

% ion concentrations (mM/L)
% K_ion_in: potassium ion concentration inside
% K_ion_out: potassium ion concentration outside
% Na_ion_in: sodium ion concentration inside
% Na_ion_out: sodium ion concentration outside
K_ion_in(1) = 150;
K_ion_out(1) = 5.5;
Na_ion_in(1) = 15;
Na_ion_out(1) = 150;

% I: applied current
% apply 15 nA current (15 mV voltage)
% start at .5mS, last for .5mS
x = 0:dt:.5;
I = zeros(1,N);
I(length(x):2*length(x)-1) = 15;

%% Differential Equations & Helper Functions

% instant opening rate constants
calc_alpha_n = @(V) 0.01*(V+55)/(1-exp(-(V+55)/10));
calc_alpha_m = @(V) 0.1*(V+40)/(1-exp(-(V+40)/10));
calc_alpha_h = @(V) 0.07*exp(-(V+65)/20);

% instant closing rate constants
calc_beta_n = @(V) 0.125*exp(-(V+65)/80);
calc_beta_m = @(V) 4*exp(-(V+65)/18);
calc_beta_h = @(V) (1/(exp(-(V+35)/10)+1));

% instant currents
calc_I_K = @(n, V) ...
            (g_K.*(n^4)*(V - V_K));
calc_I_Na = @(m, h, V) ...
            (g_Na*(m^3)*h*(V - V_Na));
calc_I_L = @(V) ...
            (g_L*(V - V_L));
calc_I_P = @(V, Nain, Kout) ...
            calc_I_L(V_rest)*(Nain > 0)*(Kout > 0);
        
% dndt: rate of change of sodium channels activation probability
% dmdt: rate of change of potassium channels activation probability
% dhdt: rate of change of leakage channels activation probability
dndt = @(t, n, V) ...
        (calc_alpha_n(V)*(1-n)) - (calc_beta_n(V)*n);
dmdt = @(t, m, V) ...
        (calc_alpha_m(V)*(1-m)) - (calc_beta_m(V)*m);
dhdt = @(t, h, V) ...
        (calc_alpha_h(V)*(1-h)) - (calc_beta_h(V)*h);

% dK_indt: rate of change of potassium concentration inside
% dK_outdt: rate of change of potassium concentration outside
% dNa_indt: rate of change of sodium concentration inside
% dNa_outdt: rate of change of sodium concentration outside
dK_indt = @(t)...
           2*calc_I_L(V_rest);
dK_outdt = @(t)...
           -2*calc_I_L(V_rest);
dNa_indt = @(t)...
           -3*calc_I_L(V_rest);
dNa_outdt = @(t)...
            3*calc_I_L(V_rest);

% dVdt: rate of change of membrane potential
% dVdt_close: leakage and concentration-sensitive Na-K pump flows
% dVdt_Na_open: leakage, concentration-sensitive Na-K pump, and Na flow
% dVdt_K_open: leakage, concentration-sensitive Na-K pump, and K flow
% dVdt_open: leakage, concentration-sensitive Na-K pump, Na and K flow
dVdt_close = @(t, V, n, m, h, I_i, Nain, Kout) ...
              (I_i - calc_I_L(V) + calc_I_P(V,Nain,Kout))/C_m;
dVdt_Na_open = @(t, V, n, m, h, I_i, Nain, Kout) ...
                (I_i - calc_I_Na(m,h,V) - calc_I_L(V) + calc_I_P(V,Nain,Kout))/C_m;
dVdt_K_open = @(t, V, n, m, h, I_i, Nain, Kout) ...
               (I_i - calc_I_K(n,V) - calc_I_L(V) + calc_I_P(V,Nain,Kout))/C_m;
dVdt_open = @(t, V, n, m, h, I_i, Nain, Kout) ...
             (I_i - calc_I_Na(m,h,V) - calc_I_K(n,V) - calc_I_L(V) ...
             + calc_I_P(V,Nain,Kout))/C_m;

%% Simulation

% main simulation loop over assigned time period
% i: simulation number 
for i = 1:N-1
    
    % initial gates conditions
    % leakage current only
    if i == 1
        Na_opened = false;
        K_opened = false;
        dVdt = dVdt_close;
    end
    
    % membrane potential greater than Na opening threshold
    % start Na current
    % NOTE: in one simulation, this should only happen once
    if V(i) >= open_thres && Na_opened == false
        dVdt = dVdt_Na_open;
        Na_opened = true;
        fprintf('sodium channel current start time: %g.mS\n',t(i))
    
    % membrane potential greater than Na closing/K opening threshold
    % close Na current, start K current
    % NOTE: in one simulation, this should only happen once
    elseif V(i) >= close_thres && K_opened == false
        dVdt = dVdt_K_open;
        K_opened = true;
        fprintf('sodium channel current close time: %g.mS\n',t(i))
        fprintf('potassium channel current open time: %g.mS\n',t(i))
    
    % during hyperpolorization, when membrane potential starts increasing
    % close K current, back to leakage ONLY
    elseif V(i) < V(1) && dV1 >= 0 && K_opened == true
        dVdt = dVdt_close;
        K_opened = false;
        fprintf('potassium channel current close time: %g.mS\n',t(i))
    end
    
    % RK4 Simulations
    t(i+1) = t(i) + dt;
    
    dKin1 = dt * dK_indt(t(i));
    dKout1 = dt * dK_outdt(t(i));
    dNain1 = dt * dNa_indt(t(i));
    dNaout1 = dt * dNa_outdt(t(i));
    dn1 = dt * dndt(t(i), n(i), V(i));
    dm1 = dt * dmdt(t(i), m(i), V(i));
    dh1 = dt * dhdt(t(i), h(i), V(i));
    dV1 = dt * dVdt(t(i), V(i), n(i), m(i), h(i), I(i), ...
                    Na_ion_in(i), K_ion_out(i));

    dKin2 = dt * dK_indt(t(i)+0.5*dt);
    dKout2 = dt * dK_outdt(t(i)+0.5*dt);
    dNain2 = dt * dNa_indt(t(i)+0.5*dt);
    dNaout2 = dt * dNa_outdt(t(i)+0.5*dt);
    dn2 = dt * dndt(t(i)+0.5*dt, n(i)+0.5*dn1, V(i)+0.5*dV1);
    dm2 = dt * dmdt(t(i)+0.5*dt, m(i)+0.5*dm1, V(i)+0.5*dV1);
    dh2 = dt * dhdt(t(i)+0.5*dt, h(i)+0.5*dh1, V(i)+0.5*dV1);
    dV2 = dt * dVdt(t(i)+0.5*dt, V(i)+0.5*dV1, ...
                    n(i)+0.5*dn1, m(i)+0.5*dm1, h(i)+0.5*dh1, I(i),...
                    Na_ion_in(i)+0.5*dNain1, K_ion_out(i)+0.5*dKout1);

    dKin3 = dt * dK_indt(t(i)+0.5*dt);
    dKout3 = dt * dK_outdt(t(i)+0.5*dt);
    dNain3 = dt * dNa_indt(t(i)+0.5*dt);
    dNaout3 = dt * dNa_outdt(t(i)+0.5*dt);
    dn3 = dt * dndt(t(i)+0.5*dt, n(i)+0.5*dn2, V(i)+0.5*dV2);
    dm3 = dt * dmdt(t(i)+0.5*dt, m(i)+0.5*dm2, V(i)+0.5*dV2);
    dh3 = dt * dhdt(t(i)+0.5*dt, h(i)+0.5*dh2, V(i)+0.5*dV2);
    dV3 = dt * dVdt(t(i)+0.5*dt, V(i)+0.5*dV2, ...
                    n(i)+0.5*dn2, m(i)+0.5*dm2, h(i)+0.5*dh2, I(i),...
                    Na_ion_in(i)+0.5*dNain2, K_ion_out(i)+0.5*dKout2);

    dKin4 = dt * dK_indt(t(i)+dt);
    dKout4 = dt * dK_outdt(t(i)+dt);
    dNain4 = dt * dNa_indt(t(i)+dt);
    dNaout4 = dt * dNa_outdt(t(i)+dt);
    dn4 = dt * dndt(t(i)+dt, n(i)+dn3, V(i)+dV3);
    dm4 = dt * dmdt(t(i)+dt, m(i)+dm3, V(i)+dV3);
    dh4 = dt * dhdt(t(i)+dt, h(i)+dh3, V(i)+dV3);
    dV4 = dt * dVdt(t(i)+dt, V(i)+dV3, ...
                    n(i)+dn3, m(i)+dm3, h(i)+dh3, I(i), ...
                    Na_ion_in(i)+dNain3, K_ion_out(i)+0.5*dKout3);
                
    K_ion_in(i+1) = K_ion_in(i) + (dKin1 + 2*dKin2 + 2*dKin3 + dKin4)/6;
    K_ion_out(i+1) = K_ion_out(i) + (dKout1 + 2*dKout2 + 2*dKout3 + dKout4)/6;
    Na_ion_in(i+1) = Na_ion_in(i) + (dNain1 + 2*dNain2 + 2*dNain3 + dNain4)/6;
    Na_ion_out(i+1) = Na_ion_out(i) + (dn1 + 2*dNaout2 + 2*dNaout3 + dNaout4)/6; 
    
    if K_ion_in(i+1) < 0
        K_ion_out(i+1) = K_ion_out(i+1) + K_ion_in(i+1);
        K_ion_in(i+1) = 0;
    end
    
    if K_ion_out(i+1) < 0
        K_ion_in(i+1) = K_ion_in(i+1) + K_ion_out(i+1);
        K_ion_out(i+1) = 0;
    end
    
    if Na_ion_in(i+1) < 0
        Na_ion_out(i+1) = Na_ion_out(i+1) + Na_ion_in(i+1);
        Na_ion_in(i+1) = 0;
    end
    
    if Na_ion_out(i+1) < 0
        Na_ion_in(i+1) = Na_ion_in(i+1) + Na_ion_out(i+1);
        Na_ion_out(i+1) = 0;
    end
    
    n(i+1) = n(i) + (dn1 + 2*dn2 + 2*dn3 + dn4)/6;
    m(i+1) = m(i) + (dm1 + 2*dm2 + 2*dm3 + dm4)/6;
    h(i+1) = h(i) + (dh1 + 2*dh2 + 2*dh3 + dh4)/6;
    V(i+1) = V(i) + (dV1 + 2*dV2 + 2*dV3 + dV4)/6;
end

%% Plot
figure;
plot(t,V)
title('Membrane Potential vs. Time')
xlabel('Time (ms)')
ylabel('Voltage (mV)')

figure;
plot(t,n,'r-',t,m,'g-',t,h,'b-')
title('Activation Probabilities vs. Time')
xlabel('Time (ms)')
legend('n','m','h')

fprintf('max membrane potential: %g.mV\n',max(V))
% figure;
% plot(t,K_ion_out,t,Na_ion_in)