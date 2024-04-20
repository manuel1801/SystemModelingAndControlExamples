% point stabilization + Multiple shooting + Runge Kutta
clear; close all; clc;

import casadi.*

dt = 0.05;   %[s]
N = 20;     % prediction horizon
sim_time = 5; % Maximum simulation time

states = SX.sym('states',4,1); n_states = length(states); %  [theta_1; theta_2; theta_1_dot; theta_2_dot]
controls = SX.sym('controls'); n_controls = length(controls);

f = Function('f',{states,controls},{DoublePendulumDynamics(states,controls)}); % nonlinear mapping function f(x,u)

U = SX.sym('U',n_controls,N);   % Decision variables (controls)

X_ref = SX.sym('XS',n_states);  % reference value
Xt = SX.sym('Xt', n_states);    % initial condition

X = SX.sym('X',n_states,(N+1)); % A vector that represents the states over the optimization problem.

obj = 0; % Objective function
g = [];  % constraints vector

Q = zeros(n_states); Q(1,1) = 10; Q(2,2) = 10; Q(3,3) = 1; Q(4,4) = 1; % weighing matrices (states)
R = 0.1; % weighing matrices (controls)

g = [g;X(:,1)-Xt]; % initial condition constraints
for k = 1:N
    obj = obj+(X(:,k)-X_ref)'*Q*(X(:,k)-X_ref) + U(:,k)'*R*U(:,k); % calculate obj    
    k1 = f(X(:,k), U(:,k));   
    k2 = f(X(:,k) + dt/2*k1, U(:,k)); 
    k3 = f(X(:,k) + dt/2*k2, U(:,k));
    k4 = f(X(:,k) + dt*k3, U(:,k)); 
    st_next=X(:,k) +dt/6*(k1 +2*k2 +2*k3 +k4);
    g = [g;X(:,k+1)-st_next]; 
end

% make the decision variable one column  vector
OPT_variables = [X(:);U(:)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', [Xt;X_ref]);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;

args.lbg(1:n_states*(N+1)) = 0; % Equality constraints
args.ubg(1:n_states*(N+1)) = 0; % Equality constraints

args.lbx(1:n_states:n_states*(N+1),1) = -inf; %state 1 lower bound
args.ubx(1:n_states:n_states*(N+1),1) = inf; %state 1 upper bound
args.lbx(2:n_states:n_states*(N+1),1) = -inf; %state 2 lower bound
args.ubx(2:n_states:n_states*(N+1),1) = inf; %state 2 upper bound
args.lbx(3:n_states:n_states*(N+1),1) = -inf; %state 3 lower bound
args.ubx(3:n_states:n_states*(N+1),1) = inf; %state 3 upper bound
args.lbx(4:n_states:n_states*(N+1),1) = -inf; %state 4 lower bound
args.ubx(4:n_states:n_states*(N+1),1) = inf; %state 4 upper bound

args.lbx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = -10; % input lower bound
args.ubx(n_states*(N+1)+1:n_controls:n_states*(N+1)+n_controls*N,1) = 10; %v input upper bound

x0 = [pi; pi; 0; 0];    % initial condition.
xs = [0; 0; 0; 0];     % Reference posture.

x_cl(:,1) = x0;

t(1) = 0;

U0 = zeros(n_controls,N);        % two control inputs for each robot
X0 = repmat(x0,1,N+1);          % initialization of the states decision variables



% Start MPC
mpciter = 0;
u_cl=[];
u_disturbed = [];

main_loop = tic;
while mpciter < sim_time / dt
    args.p   = [x0;xs]; % set the values of the parameters vector

    % initial value of the optimization variables
    args.x0  = [X0(:);U0(:)]; 

    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);

    x_ol = reshape(full(sol.x(1:n_states*(N+1))), n_states, N+1);
    u_ol = reshape(full(sol.x(n_states*(N+1)+1:end)),n_controls,N); % get controls only from the solution

    u_cl(:,mpciter+1) = u_ol(:,1);
    
    x0 = x0 + dt*f(x0,u_cl(:,mpciter+1));
    x0 = full(x0);

    x_cl(:,mpciter+2) = x0;

    U0 = [u_ol(:,2:N) u_ol(:,N)];
    X0 = [x_ol(:,2:end) x_ol(:,end)];


    t(mpciter+2) = t(mpciter+1)+dt;
    mpciter = mpciter + 1;

    drawpendulum(x0(1), x0(2));
end


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