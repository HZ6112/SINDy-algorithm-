function dx = lorenz(t,x,Beta) 
dx = [Beta(1)*(x(2)-x(1)); x(1)*(Beta(2)-x(3))-x(2); x(1)*x(2)-Beta(3)*x(3);];
Beta = [10; 28; 8/3]; % Lorenz1s parameters (chaotic) 
x0=[0; 1; 20]; % Initial condition
dt = 0.001;
tspan=dt:dt:50;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));
[t,x]=ode45(@(t,x) lorenz(t,x,Beta),tspan,x0,options); 
plot3(x(:,1),x(:,2),x(:,3));
