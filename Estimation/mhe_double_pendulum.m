% MATLAB code for Moving Horizon Estimation (MHE) of a double inverted pendulum

clear; close all; clc;

% Add path to folder containing the nonlinear dynamcs function
addpath('../Modeling/')

% Import CasADi library for symbolic computations
import casadi.*

% Parameters
dt = 0.05;              % Sampling time [s]
N_MCP = 20;                 % Prediction horizon
sim_time = 5;           % Maximum simulation time
N_sim = sim_time/dt;
v_max = 1e-1; v_min = -v_max;

w_true = [0.02; 0.04];

% System dynamics
states = SX.sym('states',4,1);     % System states [theta_1; theta_2; theta_1_dot; theta_2_dot]
n_states = length(states);         % Number of states
controls = SX.sym('controls');      % System inputs
n_controls = length(controls);      % Number of inputs
n_outputs = 2;
process_noise = SX.sym('process_noise',2,1);

% Create a function handle for the dynamics
f = Function('f',{states,controls},{DoublePendulumDynamics(states,controls)});
f_noisy = Function('f',{states,controls,process_noise},{DoublePendulumDynamicsNoise(states,controls,process_noise)});

% Decision variables
U = SX.sym('U',n_controls,N_MCP);       % Control variables
X_ref = SX.sym('XS',n_states);      % Reference state
Xt = SX.sym('Xt', n_states);        % Initial condition

X = SX.sym('X',n_states,(N_MCP+1));     % States over the optimization problem

% Objective and constraints
obj = 0;    % Objective function
g = [];     % Constraints vector

% Weighting matrices
Q = zeros(n_states);
Q(1,1) = 10; Q(2,2) = 10; Q(3,3) = 1; Q(4,4) = 1; % State weights
R = 0.1;    % Control weight

% Initial condition constraint
g = [g;X(:,1)-Xt];

% Loop over prediction horizon
for k = 1:N_MCP
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

% NLP problem
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', [Xt;X_ref]);

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
args.lbg(1:n_states*(N_MCP+1)) = 0;     % Equality constraints
args.ubg(1:n_states*(N_MCP+1)) = 0;     % Equality constraints
args.lbx(1:n_states:n_states*(N_MCP+1),1) = -inf;    % State lower bounds
args.ubx(1:n_states:n_states*(N_MCP+1),1) = inf;     % State upper bounds
args.lbx(2:n_states:n_states*(N_MCP+1),1) = -inf;    % State lower bounds
args.ubx(2:n_states:n_states*(N_MCP+1),1) = inf;     % State upper bounds
args.lbx(3:n_states:n_states*(N_MCP+1),1) = -inf;    % State lower bounds
args.ubx(3:n_states:n_states*(N_MCP+1),1) = inf;     % State upper bounds
args.lbx(4:n_states:n_states*(N_MCP+1),1) = -inf;    % State lower bounds
args.ubx(4:n_states:n_states*(N_MCP+1),1) = inf;     % State upper bounds
args.lbx(n_states*(N_MCP+1)+1:n_controls:n_states*(N_MCP+1)+n_controls*N_MCP,1) = -10;  % Control lower bounds
args.ubx(n_states*(N_MCP+1)+1:n_controls:n_states*(N_MCP+1)+n_controls*N_MCP,1) = 10;   % Control upper bounds

% Initial conditions
x0 = [pi; pi; 0; 0];    % Initial condition
xs = [0; 0; 0; 0];      % Reference posture

v = v_max*(2*rand(n_outputs,N_sim+1)-1);

x_cl(:,1) = x0;
y_cl(:,1) = x0(1:n_outputs) + v(:,1);
t(1) = 0;
U0 = zeros(n_controls,N_MCP);
X0 = repmat(x0,1,N_MCP+1);

% MPC loop
mpciter = 0;
u_cl = [];
while mpciter < N_sim
    args.p   = [x0;xs]; % Set the parameter vector
    
    % Initial value of the optimization variables
    args.x0  = [X0(:);U0(:)]; 

    % Solve optimization problem
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);

    % Extract optimal states and controls
    x_ol = reshape(full(sol.x(1:n_states*(N_MCP+1))), n_states, N_MCP+1);
    u_ol = reshape(full(sol.x(n_states*(N_MCP+1)+1:end)),n_controls,N_MCP);
    u_cl(:,mpciter+1) = u_ol(:,1);
    
    % Update states
    x0 = x0 + dt*f(x0,u_cl(:,mpciter+1));
    x0 = full(x0);
    x_cl(:,mpciter+2) = x0;
    y_cl(:,mpciter+2) = x0(1:n_outputs) + v(:,mpciter+2);
    U0 = [u_ol(:,2:N_MCP) u_ol(:,N_MCP)];
    X0 = [x_ol(:,2:end) x_ol(:,end)];
    t(mpciter+2) = t(mpciter+1)+dt;    
    mpciter = mpciter + 1;
    % drawpendulum(x0(1), x0(2));
end

% Plotting
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

figure
subplot(2,1,1)
plot(t(1:end-1),x_cl(1,1:end-1)); hold on
plot(t(1:end-1),y_cl(1,1:end-1)); 
legend({'$\theta_1$','$\tilde \theta_1$'}, Interpreter="latex")
subplot(2,1,2)
plot(t(1:end-1),x_cl(2,1:end-1));hold on
plot(t(1:end-1),y_cl(2,1:end-1));
legend({'$\theta_2$','$\tilde \theta_2$'}, Interpreter="latex")
sgtitle('Closed-loop output trajectories: Double Pendulum')



% Apply MHE using noisy memasurements and a model where the friction terms
% are unknown and treated as process noise

N_MHE = 20; % or 15 (minimum) % prediction horizon

% w_max = 0.05; w_min = 0;   % Bound on process disturbance
w_max = w_true(2)+0.0001; w_min = w_true(1)-0.0001;   % Bound on process disturbance (unknown!!!)

x_min = -10;   x_max = 10;      % Bound on the state
y_min = -10; y_max = 10;  % Bound on the output


process_noise = SX.sym('process_noise',2,1);     
n_process_noise = length(process_noise);         


x0 = [pi; pi; 0; 0] ; % + 0.001*rand(n_states,1);    % Initial guess of the state
x0_hat = x0; % Initial value for \hat x(t-N)



% Weighting matrix for the estimated process disturbance
Q = 1 / (1/12*(w_max - w_min)^2) * eye(n_process_noise);
Q = diag([50 50]);

% Weighting matrix for the estimated measurement noise
R = 1 / (1/12*(v_max - v_min)^2) * eye(n_outputs);
R = 1e3;

% Weighting matrix for prior weighting
P = 1*eye(n_states);


% Parameters
X0_hat = SX.sym('X0_hat', n_states,1);   % To store \hat x(t-N)

% To store estimated state-, and output sequences 
x_est_cl = zeros(n_states,N_sim);
y_est_cl = zeros(n_outputs,N_sim);
w_est_cl = zeros(n_process_noise,N_sim);

% Initial guess for the MHE optimization problem
X0_est_init = x0_hat; % initial state guess
W_est_init = zeros(n_process_noise,0); % process noise sequence guess


% MHE loop
for k  = 1 : N_sim

    N_MHE_t = min(k-1,N_MHE);   % Estimation horizon based on current time k

    if k-1 <= N_MHE % Rebuild optimization problem as long as estimation horizon not reached

        X = SX.sym('X', n_states,N_MHE_t+1);        % Estimated state sequence
        W = SX.sym('W', n_process_noise,N_MHE_t);   % Estimated process noise sequence
        Y_meas = SX.sym('Y', n_outputs,N_MHE_t);    % Output measurements
        U = SX.sym('U',n_controls,N_MHE_t);         % Control variables

                    
        g_x = [];
        g_y = [];        % Output constraint vector
        obj = 0;         % Objective function         

        for i = 1:N_MHE_t
            
            % Compute output symbolically
            y_hat = X(1:n_outputs,i);
        
            % Add w_k'Qw_k + v_k'Rv_k to the objective function (with v=y_meas-y_est)
            obj = obj + W(:,i)'*Q*W(:,i) + (y_hat - Y_meas(:,i))'*R*(y_hat - Y_meas(:,i));
            
            % Add output to output constraint vector
            g_y = [g_y;y_hat];
            
            % Compute next state symbolically
            x_next = f_noisy(X(:,i), U(:,i), W(:,i)); 
        
            % Add next state to state constraint vector
            g_x = [g_x; X(:,i+1) - x_next];
                
        end
        
        % Store state and output constraints in nonlinear constraint vector
        g = [g_x; g_y]; 
                
        % Add prior weighting to the cost function
        obj = obj + (X(:,1) - X0_hat)'*P*(X(:,1) - X0_hat);   % x - \hat x(t-N)
        
        % Store all optimization variables in one vector
        OPT_variables = [X(:);W(:)];
        
        % Store all parameters in one vector
        params = [X0_hat; Y_meas(:); U(:)];
        
        % Create struct object for the NLP
        nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', params);
                
        % Create solver object
        solver = nlpsol('solver', 'ipopt', nlp_prob, opts);
                
        % Specify lower and upper bounds on the nonlinear constraint vector g
        args = struct;
        args.lbg(1:n_states*N_MHE_t) = 0; %  Upper bound on estimated state
        args.ubg(1:n_states*N_MHE_t) = 0; %  Lower bound on estimated state
        args.lbg(n_states*N_MHE_t+1:n_states*N_MHE_t+n_outputs*N_MHE_t) = y_min; % Upper bound on estimated output
        args.ubg(n_states*N_MHE_t+1:n_states*N_MHE_t+n_outputs*N_MHE_t) = y_max; % Upper bound on estimated output
                
        % Specify lower and upper bounds on the optimization variables
        args.lbx(1:n_states*(N_MHE_t+1)) = x_min;
        args.ubx(1:n_states*(N_MHE_t+1)) = x_max;        
        args.lbx(n_states*(N_MHE_t+1)+1:n_states*(N_MHE_t+1)+n_process_noise*N_MHE_t) = w_min;
        args.ubx(n_states*(N_MHE_t+1)+1:n_states*(N_MHE_t+1)+n_process_noise*N_MHE_t) = w_max;        
           
    end
      
    % Write \hat x(t-N) and past N measurements into param vector
    args.p = [x0_hat; reshape(y_cl(:,k-N_MHE_t:k-1),n_outputs*N_MHE_t,1); reshape(u_cl(:,k-N_MHE_t:k-1),n_controls*N_MHE_t,1)];

    % Set initial guess for MHE optimization problem
    args.x0 = [X0_est_init(:); W_est_init(:)];

    % Solve the MHE optimization problem
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx, 'lbg', args.lbg, 'ubg', args.ubg, 'p',args.p);

    % Get estimated initial state
    X_est = reshape(sol.x(1:n_states*(N_MHE_t+1)), n_states, N_MHE_t+1);

    % Get estimated process disturbance sequence
    W_est = reshape(sol.x(n_states*(N_MHE_t+1)+1:n_states*(N_MHE_t+1)+n_process_noise*N_MHE_t), n_process_noise, N_MHE_t);
   
    % Estimate of current state x_t and process noise w_t
    x_est_cl(:,k) = full(X_est(:,end));
    if N_MHE_t > 0
        w_est_cl(:,k) = full(W_est(:,end)); 
    end
   
    % \hat x(t-N) and initial guess for next iteration
    if k > N_MHE % N_MHE reached 
        
        % Use estimate from time t-N_MHE_T
        x0_hat = x_est_cl(:,k-N_MHE_t+1);     
        
        % Use estimated sequence shifted by one        
        X0_est_init = [X_est(:,2:end) X_est(:,end)];
        W_est_init = [W_est(:,2:end) W_est(:,end)];

    else % N_MHE not reached 
        
        % Use first estimate
        x0_hat = x_est_cl(:,1);

        % Append last estimated to estimated sequence
        X0_est_init = [X_est X_est(:,end)];
        W_est_init = [W_est zeros(n_process_noise,1)];
    end

    
end


% Plotting
figure
subplot(3,1,1)
plot(t(1:end-1),x_cl([1,3],1:end-1)); hold on
plot(t(1:end-1),x_est_cl([1,3],:)); 
legend({'$\theta_1$','$\dot \theta_1$','$\hat \theta_1$','$\hat{\dot \theta_1}$'}, Interpreter="latex")
subplot(3,1,2)
plot(t(1:end-1),x_cl([2,4],1:end-1)); hold on
plot(t(1:end-1),x_est_cl([2,4],:)); 
legend({'$\theta_2$','$\dot \theta_2$','$\hat \theta_2$','$\hat{\dot \theta_2}$'}, Interpreter="latex")
subplot(3,1,3)
plot(t(1:end-1),w_est_cl(1,:)); hold on
plot(t(1:end-1),w_est_cl(2,:)); hold on
plot(t(1:end-1),w_true(1)*ones(1,N_sim), '--'); hold on
plot(t(1:end-1),w_true(2)*ones(1,N_sim), '--');
legend({'$\hat w_1$','$\hat w_2$', '$w_1$, $w_2$'}, Interpreter="latex")




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