%{
    Iris Liu
    Spring 2020    
    To run: Type script name in command line
%}


% 1-D Random Walks
W = 20;                         % Number of walks
numSteps = 20;                  % Number of steps in walk
B = 5;                          % Boundary for how far they can go
steps = zeros(numSteps,W);      % Mtx with position down rows, W columns

collisionPosition = zeros(1,W); % Mtx records step collision occurs at

for i=1:W    
    test = 0;
    steps(1,i)=0;   % The start position for every walk is zero
    
    % Loop for how many steps they take during walk
    % The next position is equal to last position plus random number
    % if our test value has been set to 1, we've already exceeded boundary
    % otherwise test if new step is exceeding boundary
    for j=2:numSteps
        nextStep = steps(j-1,i)+randn;
        if nextStep<=-B || nextStep>=B
            % New position will stay at boundary if exceeding the boundary
            nextStep = -B*(nextStep<=-B) + B*(nextStep>=B);
            collisionPosition(1,i) = j-1;
            nextStep = steps(j-1,i);
        end
        steps(j,i) = nextStep;
    end 
end

% Plot the datapoints
figure
plot(steps)
title('Position Along Walks')
xlabel('timestep')
ylabel('position')

% Calculate avg steps before collision happens
numCollisions = sum(steps(numSteps,:)==B) + sum(steps(numSteps,:)==-B);
avgSteps = sum(collisionPosition)/numCollisions;
fprintf("Average number of steps before a walker collides: %f", avgSteps);

% Record number of nonfrozen walkers at each step  
nonfrozenPeople = zeros(numSteps,1);
for i=1:numSteps
   nonfrozenPeople(i,1) = W-(sum(steps(i,:)==B) + sum(steps(i,:)==-B));
end

% Plot number of nonfrozen walkers at each timestep
% figure
% plot(nonfrozenPeople)
% title('Number of Non-Frozen Walkers at Each Step')
% xlabel('timestep')
% ylabel('Non-Frozen Walkers')
