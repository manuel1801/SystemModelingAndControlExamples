% MATLAB code for Moving Horizon Estimation (MHE) of a double inverted pendulum

clear; close all; clc;

% Add path to folder containing the nonlinear dynamcs function
addpath('../Modeling/')

% Import CasADi library for symbolic computations
import casadi.*

% Parameters
dt = 0.05;              % Sampling time [s]
sim_time = 5;           % Maximum simulation time
N_sim = sim_time/dt;

v_max = 5e-2; v_min = -v_max;
w_max = 1e-3; w_min = -w_max;

N_MHE = 30;
N_MPC = 20;                 % Prediction horizon


% System dynamics
states = SX.sym('states',4,1);     % System states [theta_1; theta_2; theta_1_dot; theta_2_dot]
n_states = length(states);         % Number of states
controls = SX.sym('controls');      % System inputs
n_controls = length(controls);      % Number of inputs
n_outputs = 2;
process_noise = SX.sym('process_noise',4,1);
n_process_noise = length(process_noise);

% Create a function handle for the dynamics
f = Function('f',{states,controls,process_noise},{DoublePendulumDynamics(states,controls)+process_noise});

% Output function
h = Function('h',{states},{states(1:2)});


v = v_max*(2*rand(n_outputs,N_sim)-1);
w = w_max*(2*rand(n_process_noise,N_sim+N_MPC)-1);
w0 = zeros(n_process_noise,1);



% Decision variables
U = SX.sym('U',n_controls,N_MPC);       % Control variables
X_ref = SX.sym('XS',n_states);      % Reference state
Xt = SX.sym('Xt', n_states);        % Initial condition

X = SX.sym('X',n_states,(N_MPC+1));     % States over the optimization problem
W = SX.sym('W',n_process_noise,N_MPC);     % process noise over the optimization problem

% Objective and constraints
obj = 0;    % Objective function
g = [];     % Constraints vector

% Control limits
u_ub = 20; u_lb = -20;

% Weighting matrices
Q = zeros(n_states);
Q(1,1) = 10; Q(2,2) = 10; Q(3,3) = 1; Q(4,4) = 1; % State weights
R = 0.1;    % Control weight

% Initial condition constraint
g = [g;X(:,1)-Xt];

% Loop over prediction horizon
for k = 1:N_MPC
    % Objective function (tracking error and control effort)
    obj = obj + (X(:,k)-X_ref)'*Q*(X(:,k)-X_ref) + U(:,k)'*R*U(:,k);
    
    % Runge-Kutta integration
    k1 = f(X(:,k), U(:,k), W(:,k));   
    k2 = f(X(:,k) + dt/2*k1, U(:,k), W(:,k)); 
    k3 = f(X(:,k) + dt/2*k2, U(:,k), W(:,k));
    k4 = f(X(:,k) + dt*k3, U(:,k), W(:,k)); 
    st_next = X(:,k) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    
    % State constraint
    g = [g;X(:,k+1)-st_next]; 
end

% Decision variables
OPT_variables = [X(:);U(:)];

% Store all parameters in one vector
params = [Xt;X_ref;W(:)];

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
args.lbx(1:n_states:n_states*(N_MPC+1),1) = -inf;    % State lower bounds
args.ubx(1:n_states:n_states*(N_MPC+1),1) = inf;     % State upper bounds
args.lbx(2:n_states:n_states*(N_MPC+1),1) = -inf;    % State lower bounds
args.ubx(2:n_states:n_states*(N_MPC+1),1) = inf;     % State upper bounds
args.lbx(3:n_states:n_states*(N_MPC+1),1) = -inf;    % State lower bounds
args.ubx(3:n_states:n_states*(N_MPC+1),1) = inf;     % State upper bounds
args.lbx(4:n_states:n_states*(N_MPC+1),1) = -inf;    % State lower bounds
args.ubx(4:n_states:n_states*(N_MPC+1),1) = inf;     % State upper bounds
args.lbx(n_states*(N_MPC+1)+1:n_controls:n_states*(N_MPC+1)+n_controls*N_MPC,1) = u_lb;  % Control lower bounds
args.ubx(n_states*(N_MPC+1)+1:n_controls:n_states*(N_MPC+1)+n_controls*N_MPC,1) = u_ub;   % Control upper bounds

% Initial conditions
x0 = [pi; pi; 0; 0];    % Initial condition
xs = [0; 0; 0; 0];      % Reference posture


t(1) = 0;
U0 = zeros(n_controls,N_MPC);
X0 = repmat(x0,1,N_MPC+1);

% MPC loop
xt = x0;
xt_noisy = xt + w_max*(2*rand(n_process_noise,1)-1);

x_cl(:,1) = xt;
x_cl_noisy(:,1) = xt_noisy;

u_cl=[];
for k  = 1 : N_sim
    args.p   = [xt_noisy;xs; reshape(w(:,k:k+N_MPC-1), n_process_noise*N_MPC,1)]; % Set the parameter vector
    
    % Initial value of the optimization variables
    args.x0  = [X0(:);U0(:)]; 

    % Solve optimization problem
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);

    % Extract optimal states and controls
    x_ol = reshape(full(sol.x(1:n_states*(N_MPC+1))), n_states, N_MPC+1);
    u_ol = reshape(full(sol.x(n_states*(N_MPC+1)+1:end)),n_controls,N_MPC);
    u_cl(:,k) = u_ol(:,1);

    y_cl(:,k) = full(h(xt_noisy)) + v(:,k);
        
    % Apply the input and get new state measurement
    k1 = f(xt, u_cl(:,k), w0);   
    k2 = f(xt + dt/2*k1, u_cl(:,k), w0); 
    k3 = f(xt + dt/2*k2, u_cl(:,k), w0);
    k4 = f(xt + dt*k3, u_cl(:,k), w0); 
    xt = xt + dt/6*(k1 + 2*k2 + 2*k3 + k4);

    k1 = f(xt_noisy, u_cl(:,k), w(:,k));   
    k2 = f(xt_noisy + dt/2*k1, u_cl(:,k), w(:,k)); 
    k3 = f(xt_noisy + dt/2*k2, u_cl(:,k), w(:,k));
    k4 = f(xt_noisy + dt*k3, u_cl(:,k), w(:,k)); 
    xt_noisy = xt_noisy + dt/6*(k1 + 2*k2 + 2*k3 + k4);

    xt = full(xt);
    xt_noisy = full(xt_noisy);

    x_cl(:,k+1) = xt;
    x_cl_noisy(:,k+1) = xt_noisy;

    U0 = [u_ol(:,2:N_MPC) u_ol(:,N_MPC)];
    X0 = [x_ol(:,2:end) x_ol(:,end)];
    t(k+1) = t(k)+dt;
    drawpendulum(xt_noisy(1), xt_noisy(2));
end

t = t(1:end-1);
x_cl = x_cl(:,1:end-1);
x_cl_noisy = x_cl_noisy(:,1:end-1);


% Plotting
figure
subplot(3,1,1)
plot(t,x_cl_noisy([1,2],:))
legend({'$\theta_1$','$\theta_3$'}, Interpreter="latex")
subplot(3,1,2)
plot(t,x_cl_noisy([3,4],:))
legend({'$\dot \theta_1$','$\dot \theta_2$'}, Interpreter="latex")
subplot(3,1,3)
plot(t,u_cl)
xlabel('$t$'); ylabel('$u$')
sgtitle('Closed-loop trajectories: Double Pendulum')


x_max = 2*ceil(max(max(x_cl_noisy))); x_min = -x_max;
y_max = 2*ceil(max(max(y_cl))); y_min = -y_max;

      

x0_hat = [pi; pi; 0; 0]; % Initial guess for \hat x(t-N)


% Weighting matrix for the estimated process disturbance
Q = 1 / (1/12*(w_max - w_min)^2) * eye(n_process_noise);
Q = 100*eye(n_process_noise);

% Weighting matrix for the estimated measurement noise
R = 1 / (1/12*(v_max - v_min)^2) * eye(n_outputs);
R = 10;

% Weighting matrix for prior weighting
P = 1*eye(n_states);


% Parameters
X0_hat = SX.sym('X0_hat', n_states,1);   % To store \hat x(t-N)

% To store estimated state-, and output sequences 
x_est_cl = zeros(n_states,N_sim);
y_est_cl = zeros(n_outputs,N_sim);

% Initial guess for the MHE optimization problem
X0_est_init = x0_hat; % initial state guess
W_est_init = zeros(n_process_noise,0); % process noise sequence guess


% MHE loop
for k  = 1 : N_sim

    N_MHE_t = min(k-1,N_MHE);   % Estimation horizon based on current time k

    if k-1 <= N_MHE % Rebuild optimization problem as long as estimation horizon not reached

        X = SX.sym('X', n_states,N_MHE_t+1);        % Estimated state sequence
        W = SX.sym('W', n_process_noise,N_MHE_t);   % Estimated process noise sequence
        U = SX.sym('W', n_controls,N_MHE_t);   % External control sequence sequence
        Y_meas = SX.sym('Y', n_outputs,N_MHE_t);    % Output measurements

                    
        g_x = [];
        g_y = [];        % Output constraint vector
        obj = 0;         % Objective function         

        for i = 1:N_MHE_t
            
            % Compute output symbolically
            y_hat = h(X(:,i));
        
            % Add w_k'Qw_k + v_k'Rv_k to the objective function (with v=y_meas-y_est)
            obj = obj + W(:,i)'*Q*W(:,i) + (y_hat - Y_meas(:,i))'*R*(y_hat - Y_meas(:,i));
            
            % Add output to output constraint vector
            g_y = [g_y;y_hat];
            
            % Compute next state symbolically
            k1 = f(X(:,i), U(:,i), W(:,i));   
            k2 = f(X(:,i) + dt/2*k1, U(:,i), W(:,i));
            k3 = f(X(:,i) + dt/2*k2, U(:,i), W(:,i));
            k4 = f(X(:,i) + dt*k3, U(:,i), W(:,i)); 
            x_next = X(:,i) +dt/6*(k1 +2*k2 +2*k3 +k4); 

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
        args.lbg(1:n_states*N_MHE_t) = 0; 
        args.ubg(1:n_states*N_MHE_t) = 0; 
        args.lbg(n_states*N_MHE_t+1:n_states*N_MHE_t+n_outputs*N_MHE_t) = y_min;
        args.ubg(n_states*N_MHE_t+1:n_states*N_MHE_t+n_outputs*N_MHE_t) = y_max;
                
        % Specify lower and upper bounds on the optimization variables
        args.lbx(1:n_states*(N_MHE_t+1)) = x_min;
        args.ubx(1:n_states*(N_MHE_t+1)) = x_max;        
        args.lbx(n_states*(N_MHE_t+1)+1:n_states*(N_MHE_t+1)+n_process_noise*N_MHE_t) = w_min;
        args.ubx(n_states*(N_MHE_t+1)+1:n_states*(N_MHE_t+1)+n_process_noise*N_MHE_t) = w_max;        
           
    end
      
    % Write \hat x(t-N) and past N measurements into param vector
    args.p = [x0_hat; reshape(y_cl(:,k-N_MHE_t:k-1),n_outputs*N_MHE_t,1);reshape(u_cl(:,k-N_MHE_t:k-1),n_controls*N_MHE_t,1)];

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
%     if N_MHE_t > 0
%         w_est_cl(:,k) = full(W_est(:,end)); 
%     end
   
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

    if mod(k,10)==0
        disp(['Iteration ', num2str(k), ' / ', num2str(N_sim)])
    end

    
end


% Plotting
figure
subplot(2,1,1)
plot(t,x_cl_noisy([1,2],:), '--'); hold on
plot(t,x_est_cl([1,2],:));
legend({'$\theta_1$','$\theta_2$','$\hat \theta_1$','$\hat \theta_2$'}, Interpreter="latex")
subplot(2,1,2)
plot(t,x_cl_noisy([3,4],:), '--'); hold on
plot(t,x_est_cl([3,4],:));
legend({'$\dot \theta_1$','$\dot \theta_2$','$\hat{\dot \theta_1}$','$\hat{\dot \theta_2}$'}, Interpreter="latex")



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