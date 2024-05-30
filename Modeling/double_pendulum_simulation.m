% This MATLAB script simulates the motion of the double pendulum given by 
% the function DoublePendulumDynamics using the Runge-Kutta 4 (RK4) integration method.

clc; clear; close all;

% Initial state (position and velocity)
x0 = [0; 0; 0; 0];

% Simulation time
sim_time = 10;

% Discretization time
delta = 0.01;

% Number of simulation steps
N = sim_time/delta;

% Input
u = 1*ones(1,N);

x = x0;
t = 0;
for k=1:N
    
   % Runge-Kutta 4 integration
   k1 = DoublePendulumDynamics(x(:,k),         u(:,k));
   k2 = DoublePendulumDynamics(x(:,k)+delta/2*k1, u(:,k));
   k3 = DoublePendulumDynamics(x(:,k)+delta/2*k2, u(:,k));
   k4 = DoublePendulumDynamics(x(:,k)+delta*k3,   u(:,k));
   x(:,k+1) = x(:,k) + delta/6*(k1+2*k2+2*k3+k4);
   t(k+1) = t(k)+delta;

   % Draw the current pendulums configuration
   drawpendulum(x(1,k),x(2,k));

end
t = t(1:end-1); x = x(:,1:end-1);

% Plot the pendulums angle and angular velocities
figure
subplot(2,1,1)
plot(t,x([1,3],:))
legend({'$\theta_1$','$\dot \theta_1$'}, Interpreter="latex")
title('$\theta_1$')
subplot(2,1,2)
plot(t,x([2,4],:))
legend({'$\theta_2$','$\dot \theta_2$'}, Interpreter="latex")
xlabel('$t$');
title('$\theta_2$')


function drawpendulum(theta1,theta2)
% This function draws the double pendulum, at the current angels theta1 and theta2
    figure(3);    
    axis equal    
    L1 = 0.5;
    L2 = 0.75;
    x1 = L1*sin(theta1); y1 = L1*cos(theta1);
    x2 = x1 + L2*sin(theta2); y2 = y1 + L2*cos(theta2);
    plot([0,x1,x2],[0,y1,y2],'o-','LineWidth',2.5,'color',[0 .447 .741], 'MarkerFaceColor',[0 .447 .741]);
    xlim([-1,1])
    ylim([-1.5,1.5])
    
end
