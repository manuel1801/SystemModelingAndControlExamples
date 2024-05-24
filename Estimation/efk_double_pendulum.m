% MATLAB code for Extended Kalman Filter (EKF) of a double inverted pendulum
% We add uniformly distributed reandom measurement noise
% As process noise we use the friction parameters r1 and r2

clear; close all; clc;

% Add path to folder containing the nonlinear dynamcs function
addpath('../Modeling/')

% Set the random number generator seed for reproducibility
rng("default")

% Import CasADi library for symbolic computations
import casadi.*

% Parameters
dt = 0.05;              % Sampling time [s]
sim_time = 3;           % Simulation time
N_sim = sim_time/dt;

v_max = 2e-2; v_min = -v_max;
w_max = 1e-4; w_min = -w_max;



alpha = 0.2;          % Deviation of initial condition to origin

% System dynamics
states = SX.sym('states',4,1);     % System states [theta_1; theta_2; theta_1_dot; theta_2_dot]
n_states = length(states);         % Number of states
controls = SX.sym('controls');      % System inputs
n_controls = length(controls);      % Number of inputs
n_outputs = 2;
process_noise = SX.sym('process_noise',4,1);
n_process_noise = length(process_noise);

% Create a function handle for the dynamics
% f = Function('f',{states,controls},{DoublePendulumDynamics(states,controls)});
f = Function('f',{states,controls,process_noise},{DoublePendulumDynamics(states,controls)+process_noise});

% Output function
h = Function('h',{states},{states(1:2)});

% Generate random noise sequences
v = v_max*(2*rand(n_outputs,N_sim)-1);
w = w_max*(2*rand(n_process_noise,N_sim)-1);
w0 = zeros(n_process_noise,1);

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


% Closed loop simulation
x_cl(:,1) = x0;
x_cl_noisy(:,1) = x0 + w_max*(2*rand(n_process_noise,1)-1);
for i = 1:N_sim

    % Compute control input using LQR
    u_cl(:,i) = -K*x_cl_noisy(:,i);
    
    y_cl(:,i) = full(h(x_cl_noisy(:,i))) + v(:,i);

    % Runge-Kutta integration for closed loop
    k1 = f(x_cl(:,i),  u_cl(:,i), w0);   
    k2 = f(x_cl(:,i) + dt/2*k1, u_cl(:,i), w0); 
    k3 = f(x_cl(:,i) + dt/2*k2, u_cl(:,i), w0);
    k4 = f(x_cl(:,i) + dt*k3, u_cl(:,i), w0); 
    x_cl(:,i+1) = full(x_cl(:,i) + dt/6*(k1 + 2*k2 + 2*k3 + k4));    

    k1 = f(x_cl_noisy(:,i),  u_cl(:,i), w(:,i));   
    k2 = f(x_cl_noisy(:,i) + dt/2*k1, u_cl(:,i), w(:,i)); 
    k3 = f(x_cl_noisy(:,i) + dt/2*k2, u_cl(:,i), w(:,i));
    k4 = f(x_cl_noisy(:,i) + dt*k3, u_cl(:,i), w(:,i)); 
    x_cl_noisy(:,i+1) = full(x_cl_noisy(:,i) + dt/6*(k1 + 2*k2 + 2*k3 + k4));    
    
%     drawpendulum(x_cl_noisy(1,i+1),x_cl_noisy(2,i+1)); % Draw pendulum at each step
end

t = t(1:end-1);
x_cl = x_cl(:,1:end-1);
x_cl_noisy = x_cl_noisy(:,1:end-1);

% Plotting
figure
subplot(3,1,1)
plot(t,x_cl_noisy([1,2],:)); hold on
plot(t,y_cl([1,2],:)); 
legend({'$\theta_1$','$\theta_2$', '$y_1$', '$y_2$'}, Interpreter="latex")
subplot(3,1,2)
plot(t,x_cl_noisy([3,4],:));
legend({'$\dot \theta_1$','$\dot \theta_2$'}, Interpreter="latex")
subplot(3,1,3)
plot(t,u_cl) 
xlabel('$t$'); ylabel('$u$')
sgtitle('Closed-loop trajectories: Double Pendulum')



% State estimation using EKF
x_hat_prior = + 0.1*(2*rand(n_states, 1)-1);
P_prior = (x_hat_prior - x0)*(x_hat_prior - x0)'; %  + 20*(2*rand(n_states,n_states)-0.5);

R = 1/3*(v_max)^2 * eye(n_outputs);      % variance of v
Q = 1/3*(v_max)^2 * eye(n_states);      % variance of v
R = 1e-2;

V = eye(n_outputs);
W = eye(n_states);

for i = 1:N_sim

    % Measurement Update (“Correct”)

    % Compute the Kalman gain
    K = P_prior*C'*inv(C*P_prior*C' + V*R*V');

    % Update estimate with measurement y    
    x_hat(:,i) = x_hat_prior + K*(y_cl(:,i) - C*x_hat_prior);

    % Update the error covariance
    P = (eye(n_states) - K*C)*P_prior;


    % Time Update ("Predict")

    % Predicte the next state
    k1 = f(x_hat(:,i), u_cl(:,i), w0);   
    k2 = f(x_hat(:,i) + dt/2*k1, u_cl(:,i), w0);
    k3 = f(x_hat(:,i) + dt/2*k2, u_cl(:,i), w0);
    k4 = f(x_hat(:,i) + dt*k3, u_cl(:,i), w0); 
    x_hat_prior = x_hat(:,i) + full(dt/6*(k1 +2*k2 +2*k3 +k4)); 

    % Predict the next error covariance

    % add linearization of A that can be evaluated at any x (not at origin)

    A = getA(x_hat(:,i),u_cl(:,i));
    P_prior = A*P*A' + W*Q*W';

    if mod(i,10)==0
        disp(['Iteration ', num2str(i), ' / ', num2str(N_sim)])
    end



end

% Plotting
figure
subplot(3,1,1)
plot(t,x_cl_noisy([1,2],:),'--'); hold on
plot(t,y_cl([1,2],:)); hold on
plot(t,x_hat([1,2],:));
legend({'$\theta_1$','$\theta_2$', '$y_1$', '$y_2$', '$\hat \theta_1$','$\hat \theta_2$'}, Interpreter="latex")
subplot(3,1,2)
plot(t,x_cl_noisy([3,4],:),'--');hold on
plot(t,x_hat([1,2],:));
legend({'$\dot \theta_1$','$\dot \theta_2$','$\hat{\dot \theta_1}$','$\hat{\dot \theta_2}$'}, Interpreter="latex")
subplot(3,1,3)
plot(t,u_cl) 
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