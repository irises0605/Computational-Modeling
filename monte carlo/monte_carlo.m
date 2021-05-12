%{
    Iris Liu
    Spring 2020
    Matlab script that estimates the value of ? using Monte Carlo simulation
    techniques involving random number generation for the estimate
    To run: Type script name in command line
%}

% a)

n = 10000;

x = rand(n,1);
y = rand(n,1);

r = sqrt(x.^2+y.^2);

calcpi = 4*sum(r<=1)/n;

fprintf('**Estimate Pi)**\n')
fprintf(['The estimated area of pi using Monte Carlo simulation '... 
        'is %.4f\n\n'],calcpi)

% (b)

plot(x(r<=1),y(r<=1),'b.');
hold on
plot(x(r>1),y(r>1),'r.');
title('Estimated Pi Value')
hold off

% (c)

results = zeros(1,1000);
for i = 1:1000
    x1 = rand(1,n);
    y1 = rand(1,n);
    r1 = sqrt(x1.^2+y1.^2);
    results(i) = 4*sum(r1<=1)/n;
end

fprintf('**Estimate Pi**\n')
fprintf(['The minimum estimated pi value among 1000 simulations'...
        'is %.4f\n'],min(results))
fprintf(['The maximum estimated pi value among 1000 simulations'...
        'is %.4f\n'],max(results))
fprintf(['The mean estimated pi value among 1000 simulations'...
        'is %.4f\n\n'],mean(results))
    
