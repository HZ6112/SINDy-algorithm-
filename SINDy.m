
%% SINDy 
%% Lorenz function defined 
function dx = lorenz(t,x,Beta) 
dx = [Beta(1)*(x(2)-x(1)); x(1)*(Beta(2)-x(3))-x(2); x(1)*x(2)-Beta(3)*x(3);];
Beta = [10; 28; 8/3]; % Lorenz1s parameters (chaotic) 
n = 3;
x0=[-8; 8; 27]; % Initial condition 
tspan=.01:.01:50;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n)); 
[t,x]=ode45(@(t,x) lorenz(t,x,Beta),tspan,x0,options);
%% Compute Derivative
for i=1:length(x)
    dx(i,:) = lorenz(0,x(i,:),Beta); 
end
%% SparsityDynamics defined 
function Xi = sparsifyDynamics(Theta,dXdt,lambda,n)
% Compute Sparse regression: sequential least squares 
Xi = Theta\dXdt; % Initial guess: Least-squares
% Lambda is our sparsification knob.
for k=1:10
    smallinds = (abs(Xi)<lambda); % Find small coefficients Xi(smallinds)=0; % and threshold
    for ind = 1:n % n is state dimension
        biginds = ~smallinds(:,ind);
% Regress dynamics onto remaining terms to find sparse Xi
        Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind); 
    end
end
%% Build library and compute sparse regression
Theta = poolData(x,n,3); % up to third order polynomials 
lambda = 0.025; % lambda is our sparsification knob. 
Xi = sparsifyDynamics(Theta,dx,lambda,n)
