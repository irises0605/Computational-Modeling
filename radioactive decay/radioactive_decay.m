%{
    Iris Liu
    Spring 2020
    Matlab script that computes and plots the function for what fraction of 
    an original quantity Q0 of radium-226 remains, as a function of years
    To run: Type script name in command line
%}

% (a)

lambda = 0.000427869;       % decay rate
Q0 = 1.0;                   % Q0 inital at t = 0
t0 = 0;                     % beginning time
T = 10000;                  % ending time

N = 10000;                  % number of steps taken
hN = (T-t0)/N;              % time step

t = zeros(1,N);             % initializes time vector
dQ = zeros(1,N);            % initializes rate of change vector

t(1) = t0;                  % set initial time
dQ(1) = -lambda*Q0;         % set initial rate of change
Q(1) = Q0;                   % set initial percent mass (100%)

for k=1:N                   
    t(k+1)=t(k)+hN;         
    Q(k+1)=Q(k)+dQ(k)*hN;
    dQ(k+1)=-lambda*Q(k);
end

plot(t,Q/Q0)              
hold on

% Qt = Q0*exp(-lambda*t);
% plot(t,Qt/Q0)

title('Radioactive Decay Part (a)')
ylabel('Q(t)/Q0 Fraction')
xlabel('t')

% (b)

T = 500;                    % time at 500 years
hN = (T-t0)/N;

for k = 1:N                   
    t(k+1)=t(k)+hN;         
    Q(k+1)=Q(k)+dQ(k)*hN;
    dQ(k+1)=-lambda*Q(k);
end
fprintf('**Exercise 1 Radioactive Decay Part (b)**\n')
fprintf('After %.0f years, %.3f of original quantity is left\n',...
        T,Q(N+1))
    
T = 5000;                   % time at 5000 years
hN = (T-t0)/N;

for k = 1:N                   
    t(k+1)=t(k)+hN;         
    Q(k+1)=Q(k)+dQ(k)*hN;
    dQ(k+1)=-lambda*Q(k);
end
fprintf('After %.0f years, %.3f of original quantity is left\n\n',...
        T,Q(N+1))

% (c)

M = 0.60;                   % Qt/Q0 fraction
Y = 0;                      % initialize year approximation

T = 10000;                  % time of 10000 years
hN = (T-t0)/N;

for k = 1:N               
    t(k+1)=t(k)+hN;         
    Q(k+1)=Q(k)+dQ(k)*hN;
    dQ(k+1)=-lambda*Q(k);
    if abs(Q(k) - M) <= 0.0001
        Y = t0 + hN*k;
    end
end

fprintf('**Exercise 1 Radioactive Decay Part (c)**\n')
fprintf(['If %.2f of original quantity is left, the radium-266 '... 
        'is around %.0f years old \n\n'],M,Y)
    