% This MATLAB script simulates the motion of the double pendulum given by 
% the function DoublePendulumDynamics using the Runge-Kutta 4 (RK4) integration method.

clc; clear; close all;

% Initial state (position and velocity)
x0 = zeros(6,1);

% Simulation time
sim_time = 4;

% Discretization time
delta = 0.01;

% Number of simulation steps
N = sim_time/delta;

% Input
u = [1*ones(1,N/2) -1*ones(1,N/2)];

x = x0;
t = 0;
for k=1:N
    
   % Runge-Kutta 4 integration
   k1 = DoublePendulumCartDynamics(x(:,k),         u(:,k));
   k2 = DoublePendulumCartDynamics(x(:,k)+delta/2*k1, u(:,k));
   k3 = DoublePendulumCartDynamics(x(:,k)+delta/2*k2, u(:,k));
   k4 = DoublePendulumCartDynamics(x(:,k)+delta*k3,   u(:,k));
   x(:,k+1) = x(:,k) + delta/6*(k1+2*k2+2*k3+k4);
   t(k+1) = t(k)+delta;

   % Draw the current pendulums configuration
   drawpendulum_on_cart(x(1,k), x(2,k),x(3,k));

end
t = t(1:end-1); x = x(:,1:end-1);

% Plot the pendulums angle and angular velocities
figure
subplot(3,1,1)
plot(t,x([1,4],:))
legend({'$\theta_0$','$\dot \theta_0$'}, Interpreter="latex")
title('$\theta_0$')

subplot(3,1,2)
plot(t,x([2,5],:))
legend({'$\theta_1$','$\dot \theta_1$'}, Interpreter="latex")
xlabel('$t$');
title('$\theta_1$')

subplot(3,1,3)
plot(t,x([3,6],:))
legend({'$\theta_2$','$\dot \theta_2$'}, Interpreter="latex")
xlabel('$t$');
title('$\theta_2$')



function drawpendulum_on_cart(x_cart, theta1, theta2)
    % This function draws the double pendulum on a cart, at the current 
    % position x_cart and angles theta1 and theta2
    figure(3);
    clf;  % Clear current figure
    hold on;  % Hold on to draw multiple items
    
    % Pendulum parameters
    L1 = 0.5;
    L2 = 0.75;
    
    % Cart parameters
    cart_width = 0.4;
    cart_height = 0.2;
    
    % Calculate pendulum positions
    x1 = x_cart + L1*sin(theta1); 
    y1 = L1*cos(theta1);
    x2 = x1 + L2*sin(theta2); 
    y2 = y1 + L2*cos(theta2);
    
    % Draw cart
    rectangle('Position', [x_cart - cart_width/2, -cart_height/2, cart_width, cart_height], 'Curvature', 0.1, 'FaceColor', [0.6 0.6 0.6]);
    
    % Draw pendulum
    plot([x_cart, x1, x2], [0, y1, y2], 'o-', 'LineWidth', 2.5, 'color', [0 .447 .741], 'MarkerFaceColor', [0 .447 .741]);
    
    % Set plot limits
    xlim([-2, 2]);
    ylim([-1.5, 1.5]);
    
    % Add grid and labels
    grid on;
    xlabel('X Position');
    ylabel('Y Position');
    
    hold off;  % Release the hold
end


function drawpendulum_on_cart2(x_cart, theta1, theta2)
    % This function draws the double pendulum on a cart, at the current 
    % position x_cart and angles theta1 and theta2
    % It keeps the cart always in the center of the figure to avoid moving
    % out of the limits

    figure(3);
    clf;  % Clear current figure
    hold on;  % Hold on to draw multiple items
    
    % Pendulum parameters
    L1 = 0.5;
    L2 = 0.75;
    
    % Cart parameters
    cart_width = 0.4;
    cart_height = 0.2;
    
    % Calculate pendulum positions
    x1 = x_cart + L1*sin(theta1); 
    y1 = L1*cos(theta1);
    x2 = x1 + L2*sin(theta2); 
    y2 = y1 + L2*cos(theta2);
    
    % Draw cart
    rectangle('Position', [x_cart - cart_width/2, -cart_height/2, cart_width, cart_height], 'Curvature', 0.1, 'FaceColor', [0.6 0.6 0.6]);
    
    % Draw pendulum
    plot([x_cart, x1, x2], [0, y1, y2], 'o-', 'LineWidth', 2.5, 'color', [0 .447 .741], 'MarkerFaceColor', [0 .447 .741]);
    
    % Set dynamic plot limits
    x_buffer = 1;  % Buffer space around the cart's position
    xlim([x_cart - x_buffer, x_cart + x_buffer]);
    ylim([-1.5, 1.5]);
    
    % Add grid and labels
    grid on;
    xlabel('X Position');
    ylabel('Y Position');
    
    hold off;  % Release the hold
end