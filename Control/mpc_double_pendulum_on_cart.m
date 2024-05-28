% MATLAB code for Model Predictive Control (MPC) of a double inverted
% pendulum on a cart.

clear; close all; clc;

% Add path to folder containing the nonlinear dynamcs function
addpath('../Modeling/')

% Import CasADi library for symbolic computations
import casadi.*

% Parameters
dt = 0.05;              % Sampling time [s]
N_MPC = 25;                 % Prediction horizon
sim_time = 8;           % Maximum simulation time
N_sim = sim_time/dt;

x_lim = 5;              % constraints for the carts position

% System dynamics
states = SX.sym('states',6,1);     % System states [theta_1; theta_2; theta_1_dot; theta_2_dot]
n_states = length(states);         % Number of states
controls = SX.sym('controls');      % System inputs
n_controls = length(controls);      % Number of inputs

% Create a function handle for the dynamics
f = Function('f',{states,controls},{DoublePendulumCartDynamics(states,controls)});

% Decision variables
U = SX.sym('U',n_controls,N_MPC);       % Control variables
X_ref = SX.sym('XS',n_states);      % Reference state
Xt = SX.sym('Xt', n_states);        % Initial condition

X = SX.sym('X',n_states,(N_MPC+1));     % States over the optimization problem

% Objective and constraints
obj = 0;    % Objective function
g = [];     % Constraints vector

% Control limits
u_ub = 100; u_lb = -100;

% Weighting matrices
Q = blkdiag(200,1e3*eye(2),eye(3));

R = 0.1;    % Control weight

% Initial condition constraint
g = [g;X(:,1)-Xt];

% Loop over prediction horizon
for k = 1:N_MPC
    % Objective function (tracking error and control effort)
    obj = obj + (X(:,k)-X_ref)'*Q*(X(:,k)-X_ref) + U(:,k)'*R*U(:,k);
    
    % Runge-Kutta integration
    k1 = f(X(:,k), U(:,k));   
    k2 = f(X(:,k) + dt/2*k1, U(:,k)); 
    k3 = f(X(:,k) + dt/2*k2, U(:,k));
    k4 = f(X(:,k) + dt*k3, U(:,k)); 
    st_next = X(:,k) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    
    % State constraint
    g = [g;X(:,k+1)-st_next]; 
end

% Decision variables
OPT_variables = [X(:);U(:)];

% Store all parameters in one vector
params = [Xt;X_ref];

% NLP problem
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', params);

% Solver options
opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level = 0;         % Set to 0 for minimal output
opts.print_time = 0;
opts.ipopt.acceptable_tol = 1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

% Initialize solver
solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

% Initialize constraints and bounds
args = struct;
args.lbg(1:n_states*(N_MPC+1)) = 0;     % Equality constraints
args.ubg(1:n_states*(N_MPC+1)) = 0;     % Equality constraints
args.lbx(1:n_states:n_states*(N_MPC+1),1) = -x_lim;    % State lower bounds
args.ubx(1:n_states:n_states*(N_MPC+1),1) = x_lim;     % State upper bounds
args.lbx(2:n_states:n_states*(N_MPC+1),1) = -inf;    % State lower bounds
args.ubx(2:n_states:n_states*(N_MPC+1),1) = inf;     % State upper bounds
args.lbx(3:n_states:n_states*(N_MPC+1),1) = -inf;    % State lower bounds
args.ubx(3:n_states:n_states*(N_MPC+1),1) = inf;     % State upper bounds
args.lbx(4:n_states:n_states*(N_MPC+1),1) = -inf;    % State lower bounds
args.ubx(4:n_states:n_states*(N_MPC+1),1) = inf;     % State upper bounds
args.lbx(5:n_states:n_states*(N_MPC+1),1) = -inf;    % State lower bounds
args.ubx(5:n_states:n_states*(N_MPC+1),1) = inf;     % State upper bounds
args.lbx(6:n_states:n_states*(N_MPC+1),1) = -inf;    % State lower bounds
args.ubx(6:n_states:n_states*(N_MPC+1),1) = inf;     % State upper bounds
args.lbx(n_states*(N_MPC+1)+1:n_controls:n_states*(N_MPC+1)+n_controls*N_MPC,1) = u_lb;  % Control lower bounds
args.ubx(n_states*(N_MPC+1)+1:n_controls:n_states*(N_MPC+1)+n_controls*N_MPC,1) = u_ub;   % Control upper bounds

% Initial conditions
x0 = [0; pi; pi; 0; 0; 0];    % Initial condition
xs = [0; 0; 0; 0; 0; 0];      % Reference posture

x_cl(:,1) = x0;
t(1) = 0;
U0 = zeros(n_controls,N_MPC);
X0 = repmat(x0,1,N_MPC+1);

% MPC loop
xt = x0;
u_cl=[];
for k  = 1 : N_sim
    args.p   = [xt;xs]; % Set the parameter vector
    
    % Initial value of the optimization variables
    args.x0  = [X0(:);U0(:)]; 

    % Solve optimization problem
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);

    % Extract optimal states and controls
    x_ol = reshape(full(sol.x(1:n_states*(N_MPC+1))), n_states, N_MPC+1);
    u_ol = reshape(full(sol.x(n_states*(N_MPC+1)+1:end)),n_controls,N_MPC);
    u_cl(:,k) = u_ol(:,1);
        
    % Apply the input and get new state measurement
    k1 = f(xt, u_cl(:,k));   
    k2 = f(xt + dt/2*k1, u_cl(:,k)); 
    k3 = f(xt + dt/2*k2, u_cl(:,k));
    k4 = f(xt + dt*k3, u_cl(:,k)); 
    xt = xt + dt/6*(k1 + 2*k2 + 2*k3 + k4);

    xt = full(xt);
    x_cl(:,k+1) = xt;
    U0 = [u_ol(:,2:N_MPC) u_ol(:,N_MPC)];
    X0 = [x_ol(:,2:end) x_ol(:,end)];
    t(k+1) = t(k)+dt;
    drawpendulum_on_cart(xt(1), xt(2), xt(3));
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
    xlim([-5, 5]);
    ylim([-1.5, 1.5]);
    
    % Add grid and labels
    grid on;
    xlabel('X Position');
    ylabel('Y Position');
    
    hold off;  % Release the hold
end


