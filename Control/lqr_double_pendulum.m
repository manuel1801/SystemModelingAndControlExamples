% MATLAB code for Linear–quadratic regulator (LQR) of a double inverted pendulum

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

alpha = 0.2;          % Deviation of initial condition to origin
% (for 0.2 the origin can be stabilized, for 0.3 not)

% System dynamics
states = SX.sym('states',4,1);     % System states [theta_1; theta_2; theta_1_dot; theta_2_dot]
n_states = length(states);         % Number of states
controls = SX.sym('controls');      % System inputs
n_controls = length(controls);      % Number of inputs

% Create a function handle for the dynamics
f = Function('f',{states,controls},{DoublePendulumDynamics(states,controls)});

% System parameters
m1 = 0.75;      % Mass 1
m2 = 1;         % Mass 2
g = 9.81;       % Gravity
L1 = 0.75;      % Length 1
L2 = 1;         % Length 2
l1 = 1/2*L1;    % Half-length 1
l2 = 1/2*L2;    % Half-length 2
J1 = m1*L1^2/12;    % Inertia 1
J2 = m2*L2^2/12;    % Inertia 2
r1 = 0.02; r2 = 0.04;   % Frictrion parameters

% Combined parameter
p = J2*m2*L1^2 + m1*m2*l1^2*l2^2 + J2*m1*l1^2 + J1*m2*l2^2 + J1*J2;
 
% System matrices
A=[
                                   0,                                    0,                      1,                                0
                                   0,                                    0,                      0,                                1
(g*(m2*l2^2 + J2)*(L1*m2 + l1*m1))/p,                  -(L1*g*l2^2*m2^2)/p, -(r1*(m2*l2^2 + J2))/p,                  (L1*l2*m2*r2)/p
     -(L1*g*l2*m2*(L1*m2 + l1*m1))/p, (g*l2*m2*(m2*L1^2 + m1*l1^2 + J1))/p,        (L1*l2*m2*r1)/p, -(r2*(m2*L1^2 + m1*l1^2 + J1))/p];
 
B=[
               0
               0
(m2*l2^2 + J2)/p
   -(L1*l2*m2)/p];


% Output matrices
C = [1 0 0 0;
     0 1 0 0];

D = 0;

% Weighting matrices for LQR
Q = zeros(n_states);
Q(1,1) = 10; Q(2,2) = 10; Q(3,3) = 1; Q(4,4) = 1; % State weights
R = 0.1;    % Control weight

% Compute LQR gain matrix
K = lqr(A,B,Q,R);

% Generate initial conditions with some deviation to the origin
x0 = zeros(n_states,1) + alpha*rand(n_states, 1);

t = 0:dt:sim_time;


% % Open loop simulation
% x_ol(:,1) = x0;
% for i = 1:sim_time/dt
% 
%     % Runge-Kutta integration for open loop
%     k1 = f(x_ol(:,i), 0);   
%     k2 = f(x_ol(:,i) + dt/2*k1, 0); 
%     k3 = f(x_ol(:,i) + dt/2*k2, 0);
%     k4 = f(x_ol(:,i) + dt*k3, 0); 
%     x_ol(:,i+1) = full(x_ol(:,i) + dt/6*(k1 + 2*k2 + 2*k3 + k4));    
%     drawpendulum(x_ol(1,i+1),x_ol(2,i+1)); % Draw pendulum at each step
% end
% figure
% subplot(2,1,1)
% plot(t(1:end-1),x_ol([1,3],1:end-1))
% legend({'$\theta_1$','$\dot \theta_1$'}, Interpreter="latex")
% subplot(2,1,2)
% plot(t(1:end-1),x_ol([2,4],1:end-1))
% legend({'$\theta_2$','$\dot \theta_2$'}, Interpreter="latex")
% sgtitle('Open-loop trajectories: Double Pendulum')


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
    drawpendulum(x_cl(1,i+1),x_cl(2,i+1)); % Draw pendulum at each step
end

% Closed-Loop
figure
subplot(3,1,1)
plot(t(1:end-1),x_cl([1,3],1:end-1))
legend({'$\theta_1$','$\dot \theta_1$'}, Interpreter="latex")
subplot(3,1,2)
plot(t(1:end-1),x_cl([2,4],1:end-1))
legend({'$\theta_2$','$\dot \theta_2$'}, Interpreter="latex")
subplot(3,1,3)
plot(t(1:end-1),u_cl)
xlabel('$t$'); ylabel('$u$')
sgtitle('Closed-loop trajectories: Double Pendulum')


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