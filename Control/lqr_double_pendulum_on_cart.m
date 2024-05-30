% MATLAB code for Linear–quadratic regulator (LQR) of a double inverted
% pendulum on a cart

clear; close all; clc;

% Add path to folder containing the nonlinear dynamcs function
addpath('../Modeling/')

% Set the random number generator seed for reproducibility
rng("default")

% Import CasADi library for symbolic computations
import casadi.*

% Parameters
dt = 0.01;              % Sampling time [s]
sim_time = 3;           % Simulation time

alpha = 0.4;          % Deviation of initial condition to origin
% (for 0.2 the origin can be stabilized, for 0.3 not)

% System dynamics
states = SX.sym('states',6,1);     % System states [theta_1; theta_2; theta_1_dot; theta_2_dot]
n_states = length(states);         % Number of states
controls = SX.sym('controls');      % System inputs
n_controls = length(controls);      % Number of inputs

% Create a function handle for the dynamics
f = Function('f',{states,controls},{DoublePendulumCartDynamics(states,controls)});

% System parameters
m0 = 1.5;       % Mass of the cart
m1 = 0.75;      % Mass of first link
m2 = 1;         % Mass of second link
g = 9.81;       % Gravity
L1 = 0.75;      % Length 1
L2 = 1;         % Length 2
l1 = 1/2*L1;    % Half-length 1
l2 = 1/2*L2;    % Half-length 2
J1 = m1*L1^2/12;    % Inertia 1
J2 = m2*L2^2/12;    % Inertia 2
r1 = 0.02; r2 = 0.04;   % Frictrion parameters

p = (J2*L1^2*m2^2 + J2*l1^2*m1^2 + J1*l2^2*m2^2 + 2*J1*J2*m0 + 2*J1*J2*m1 + 2*J1*J2*m2 + L1^2*l2^2*m0*m2^2 + L1^2*l2^2*m1*m2^2 + 2*J2*L1^2*m0*m2 + 2*J2*L1^2*m1*m2 + l1^2*l2^2*m1*m2^2 + l1^2*l2^2*m1^2*m2 + 2*J2*l1^2*m0*m1 + 2*J1*l2^2*m0*m2 + 2*J1*l2^2*m1*m2 + 2*J2*l1^2*m1*m2 - J2*L1^2*m2^2 - J2*l1^2*m1^2 - J1*l2^2*m2^2 - L1*l1*l2^2*m1*m2^2 - 2*J2*L1*l1*m1*m2 + 2*l1^2*l2^2*m0*m1*m2 - l1^2*l2^2*m1^2*m2 - l1^2*l2^2*m1*m2^2 - L1^2*l2^2*m0*m2^2 - L1^2*l2^2*m1*m2^2 - L1*l1*l2^2*m1*m2^2 + L1*l1*l2^2*m1*m2^2 - 2*J2*L1*l1*m1*m2 + L1*l1*l2^2*m1*m2^2);

A= [0,                                                                         0,                                                                                                         0, 1,                                                           0,                                                                                                     0;
    0,                                                                         0,                                                                                                         0, 0,                                                           1,                                                                                                     0;
    0,                                                                         0,                                                                                                         0, 0,                                                           0,                                                                                                     1;
    0,            -(2*g*(L1*m2 + l1*m1)*(l1*m1*m2*l2^2 + J2*L1*m2 + J2*l1*m1))/p,                                                              -(2*g*l2^2*m2^2*(m1*l1^2 - L1*m1*l1 + J1))/p, 0,              (2*r1*(l1*m1*m2*l2^2 + J2*L1*m2 + J2*l1*m1))/p,                                                              (2*l2*m2*r2*(m1*l1^2 - L1*m1*l1 + J1))/p;
    0, (2*g*(L1*m2 + l1*m1)*(J2*m0 + J2*m1 + J2*m2 + l2^2*m0*m2 + l2^2*m1*m2))/p,                                                                -(2*g*l2^2*m2^2*(L1*m0 + L1*m1 - l1*m1))/p, 0, -(2*r1*(J2*m0 + J2*m1 + J2*m2 + l2^2*m0*m2 + l2^2*m1*m2))/p,                                                                (2*l2*m2*r2*(L1*m0 + L1*m1 - l1*m1))/p;
    0,                    -(2*g*l2*m2*(L1*m2 + l1*m1)*(L1*m0 + L1*m1 - l1*m1))/p, (2*g*l2*m2*(J1*m0 + J1*m1 + J1*m2 + l1^2*m0*m1 + l1^2*m1*m2 + L1^2*m0*m2 + L1^2*m1*m2 - 2*L1*l1*m1*m2))/p, 0,                      (2*l2*m2*r1*(L1*m0 + L1*m1 - l1*m1))/p, -(2*r2*(J1*m0 + J1*m1 + J1*m2 + l1^2*m0*m1 + l1^2*m1*m2 + L1^2*m0*m2 + L1^2*m1*m2 - 2*L1*l1*m1*m2))/p];
 
B= [                                                                       0;
                                                                           0;
                                                                           0;
(2*J2*m2*L1^2 + 2*m1*m2*l1^2*l2^2 + 2*J2*m1*l1^2 + 2*J1*m2*l2^2 + 2*J1*J2)/p;
                              -(2*l1*m1*m2*l2^2 + 2*J2*L1*m2 + 2*J2*l1*m1)/p;
                                      -(2*l2*m2*(m1*l1^2 - L1*m1*l1 + J1))/p];


% Output matrices
C = [eye(3) zeros(3)];


% Weighting matrices for LQR
Q = blkdiag(50,10*eye(2),eye(3));
R = 0.1;    % Control weight

% Compute LQR gain matrix
K = lqr(A,B,Q,R);

% Generate initial conditions with some deviation to the origin
x0 = zeros(n_states,1) + alpha*rand(n_states, 1);
x0(1) = -1.5;

t = 0:dt:sim_time;

% Closed loop simulation
x_cl(:,1) = x0;
for i = 1:sim_time/dt

    % Compute control input using LQR
    u_cl(:,i) = -K*x_cl(:,i);

    % Runge-Kutta integration for closed loop
    k1 = f(x_cl(:,i),  u_cl(:,i));   
    k2 = f(x_cl(:,i) + dt/2*k1, u_cl(:,i)); 
    k3 = f(x_cl(:,i) + dt/2*k2, u_cl(:,i));
    k4 = f(x_cl(:,i) + dt*k3, u_cl(:,i)); 
    x_cl(:,i+1) = full(x_cl(:,i) + dt/6*(k1 + 2*k2 + 2*k3 + k4));    
    drawpendulum_on_cart(x_cl(1,i+1),x_cl(2,i+1),x_cl(3,i+1)); % Draw pendulum at each step
end

% Closed-Loop
figure
subplot(4,1,1)
plot(t(1:end-1),x_cl([1,4],1:end-1))
legend({'$\theta_1$','$\dot \theta_1$'}, Interpreter="latex")
subplot(4,1,2)
plot(t(1:end-1),x_cl([2,5],1:end-1))
legend({'$\theta_2$','$\dot \theta_2$'}, Interpreter="latex")
subplot(4,1,3)
plot(t(1:end-1),x_cl([3,6],1:end-1))
legend({'$\theta_3$','$\dot \theta_3$'}, Interpreter="latex")
subplot(4,1,4)
plot(t(1:end-1),u_cl)
xlabel('$t$'); ylabel('$u$')
sgtitle('Closed-loop trajectories: Double Pendulum')



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