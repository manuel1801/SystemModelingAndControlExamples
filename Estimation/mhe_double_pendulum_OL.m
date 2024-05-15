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
v_max = 1e-1; v_min = -v_max;
w_max = 1e-4; w_min = -w_max;
x_min = -10;   x_max = 10;      % Bound on the state
y_min = -10; y_max = 10;  % Bound on the output

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



% Solver options
opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level = 0;         % Set to 0 for minimal output
opts.print_time = 0;
opts.ipopt.acceptable_tol = 1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;


% Initial conditions
x0 = [0.01; 0.01; 0; 0];    % Initial condition

v = v_max*(2*rand(n_outputs,N_sim+1)-1);
w = w_max*(2*rand(n_states,N_sim+1)-1);

x_OL(:,1) = x0;
x_OL_noisy(:,1) = x0;
y_OL(:,1) = x0(1:n_outputs) + v(:,1);
t(1) = 0;
u = zeros(n_controls,1);
for k=1:N_sim
    
   % Runge-Kutta 4 integration
%    k1 = f(x_OL(:,k),         u);
%    k2 = f(x_OL(:,k)+dt/2*k1, u);
%    k3 = f(x_OL(:,k)+dt/2*k2, u);
%    k4 = f(x_OL(:,k)+dt*k3,   u);
%    x_OL(:,k+1) = x_OL(:,k) + full(dt/6*(k1+2*k2+2*k3+k4));

   x_OL(:,k+1) = x_OL(:,k) + full(dt*f(x_OL(:,k),u));
   x_OL_noisy(:,k+1) = x_OL_noisy(:,k) + full(dt*f(x_OL_noisy(:,k),u)) + w(:,k);
   y_OL(:,k+1) = x_OL_noisy(1:n_outputs,k) + v(:,k);

   t(k+1) = t(k)+dt;

   % Draw the current pendulums configuration
%    drawpendulum(x_OL(1,k),x_OL(2,k));

end

% Plotting
figure
subplot(3,1,1)
plot(t,x_OL([1,2],:)); hold on
plot(t,x_OL_noisy([1,2],:)); 
legend({'$\theta_1$','$\theta_2$'}, Interpreter="latex")
subplot(3,1,2)
plot(t,x_OL([3,4],:));hold on
plot(t,x_OL_noisy([3,4],:)); 
legend({'$\dot \theta_1$','$\dot \theta_2$'}, Interpreter="latex")
subplot(3,1,3)
plot(t,x_OL([1,2],:)); hold on
plot(t,y_OL); 
legend({'$\theta_1$','$\theta_2$'}, Interpreter="latex")




N_MHE = 20;

process_noise = SX.sym('process_noise',n_states,1);     
n_process_noise = length(process_noise);         

meas_noise = SX.sym('process_noise',n_states,1);     
n_meas_noise = length(meas_noise);         

x0 = [0; 0; 0; 0];
x0_hat = x0; % Initial value for \hat x(t-N)



% Weighting matrix for the estimated process disturbance
Q = 1 / (1/12*(w_max - w_min)^2) * eye(n_process_noise);
Q = 1e3*eye(n_states);

% Weighting matrix for the estimated measurement noise
R = 1 / (1/12*(v_max - v_min)^2) * eye(n_outputs);
R = 1e2;

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
            x_next = f(X(:,i), 0) + W(:,i); 
        
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
        params = [X0_hat; Y_meas(:)];
        
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
    args.p = [x0_hat; reshape(y_OL(:,k-N_MHE_t:k-1),n_outputs*N_MHE_t,1)];

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

    
end


% Plotting
figure
subplot(2,1,1)
plot(t(1:end-1),x_OL_noisy([1,2],1:end-1)); hold on
plot(t(1:end-1),x_est_cl([1,2],:));
legend({'$\theta_1$','$\theta_2$'}, Interpreter="latex")
subplot(2,1,2)
plot(t(1:end-1),x_OL_noisy([3,4],1:end-1)); hold on
plot(t(1:end-1),x_est_cl([3,4],:));
legend({'$\dot \theta_1$','$\dot \theta_2$'}, Interpreter="latex")




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