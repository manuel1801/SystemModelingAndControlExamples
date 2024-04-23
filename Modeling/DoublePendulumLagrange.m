% This MATLAB script uses the Euler-Lagrange formalism to derive the equations
% of motion for a double pendulum with a fixed base. It transforms the equations
% into a nonlinear state space representation and subsequently linearizes them
% around the upright position of the pendulum.

clc; clear; close all;

% Symbolic variables for the double pendulumn's parameters
syms m1 m2 L1 L2 l1 l2 J1 J2 g r1 r2 p 
params_sym = {m1 m2 L1 L2 l1 l2 J1 J2 g r1 r2 p};

% Symbolic variables for the state q=[q1,q2], its 1st and 2nd time
% derivatives q_dot=[q1_dot,q2_dot] and q_ddot=[q1_ddot,q2_ddot]
% and the input
syms q1 q1_dot q2 q2_dot q1_ddot q2_ddot u

% Numeric values for the parameters
m1_num  = 0.5;
m2_num  = 0.75;
g_num  = 9.81;
L1_num = 0.5;
L2_num = 0.75;
l1_num   = 1/2*L1_num;
l2_num   = 1/2*L2_num;
J1_num = m1_num*L1_num^2/12;
J2_num = m2_num*L2_num^2/12;
r1_num = 0.1;
r2_num = 0.1;
p_num = J2_num*m2_num*L1_num^2 + m1_num*m2_num*l1_num^2*l2_num^2 + J2_num*m1_num*l1_num^2 + J1_num*m2_num*l2_num^2 + J1_num*J2_num;
params_num = {m1_num m2_num L1_num L2_num l1_num l2_num J1_num J2_num g_num r1_num r2_num p_num};

% Lagrange Function for the double pendulum
L = 1/2*(m1*l1^2 + m2*L1^2+J1)*q1_dot^2 ...
    + 1/2*(m2*l2^2 + J2)*q2_dot^2 ...
    + m2*L1*l2*cos(q1 - q2)*q1_dot*q2_dot ...
    - (m1*l1 + m2*L1)*g*cos(q1) ...
    - m2*l2*g*cos(q2);

% Partial derivatives w.r.t. q
dLdq1 = diff(L,q1);
dLdq2 = diff(L,q2);

% Partial derivatives w.r.t. q_dot
dLdq1_dot = diff(L,q1_dot);
dLdq2_dot = diff(L,q2_dot);

% Define symbolic functions for q1 and q2
syms q1_f(t) q2_f(t)

% Get time derivative of q1 and q2 (also as symbolic function)
q1_dot_f(t)= diff(q1_f,t);
q2_dot_f(t)= diff(q2_f,t);

% Substitute symbolic variables with symbolic functions
dLdq1_f = subs(dLdq1,{q1 q2 q1_dot q2_dot},{q1_f q2_f q1_dot_f q2_dot_f});
dLdq2_f = subs(dLdq2,{q1 q2 q1_dot q2_dot},{q1_f q2_f q1_dot_f q2_dot_f});
dLdq1_dot_f = subs(dLdq1_dot,{q1 q2 q1_dot q2_dot},{q1_f q2_f q1_dot_f q2_dot_f});
dLdq2_dot_f = subs(dLdq2_dot,{q1 q2 q1_dot q2_dot},{q1_f q2_f q1_dot_f q2_dot_f});

% Now, take the time derivative of dLdq_dot
dLdq1_ddot_f = diff(dLdq1_dot_f,t);
dLdq2_ddot_f = diff(dLdq2_dot_f,t);

% The final euler lagrange equation
euler_lagrange_eq_1 = dLdq1_ddot_f - dLdq1_f - u;
euler_lagrange_eq_2 = dLdq2_ddot_f - dLdq2_f;
euler_lagrange_eq_1 = simplify(euler_lagrange_eq_1);
euler_lagrange_eq_2 = simplify(euler_lagrange_eq_2);
euler_lagrange_eq = [euler_lagrange_eq_1;euler_lagrange_eq_2];

% Backsubstitute the symbolic variables for the symbolic functions
euler_lagrange_eq = subs(euler_lagrange_eq,diff(q1_f(t), t, t),q1_ddot);
euler_lagrange_eq = subs(euler_lagrange_eq,diff(q2_f(t), t, t),q2_ddot);
euler_lagrange_eq = subs(euler_lagrange_eq,diff(q1_f(t), t),q1_dot);
euler_lagrange_eq = subs(euler_lagrange_eq,diff(q2_f(t), t),q2_dot);
euler_lagrange_eq = subs(euler_lagrange_eq,q1_f(t),q1);
euler_lagrange_eq = subs(euler_lagrange_eq,q2_f(t),q2);

disp('The Equation of motion for the double pendulum are:')
disp(euler_lagrange_eq )

% Bring them into the form M(q)*q_ddot + D(q,q_dot)+C(q) + R*q_dot = H*u

% Positive definite intertia matrix
M = jacobian(euler_lagrange_eq,[q1_ddot;q2_ddot]);

% Gravitational forces (vecor)
C = subs(euler_lagrange_eq,{q1_dot q2_dot q1_ddot q2_ddot u},{0 0 0 0 0});

% Input matrix
H = [1;0];

% Centripetal and coriolis forces (vector)
D = euler_lagrange_eq -M*[q1_ddot;q2_ddot] - C + H*u;
D = simplify(D);

% Friction terms
R = [r1 0;0 r2];

q = [q1;q2]; q_dot = [q1_dot;q2_dot];
x = [q1;q2;q1_dot;q2_dot];

% Inverse of M
M_inv = inv(M);

% Bring the equation of motion into a state space form (1st Order ODE)
% x_dot = f(x,u) with the state x = [q;q_dot]
double_pendulum_ode = [q_dot;-M_inv*D-M_inv*C-M_inv*R*q_dot+M_inv*H*u];
double_pendulum_ode = simplify(double_pendulum_ode);

% Linearization of f around the origin (up right position of the pendulum)
% Compute jacobian of f w.r.t. x and u
A = jacobian(double_pendulum_ode,x);
B = jacobian(double_pendulum_ode,u);

% Insert the origin x = [0;0;0;0]
A = subs(A,{q1 q2 q1_dot q2_dot},{0 0 0 0});
A = simplify(A);

B = subs(B,{q1 q2 q1_dot q2_dot},{0 0 0 0});
B = simplify(B);

% For better readability use p := J2*m2*L1^2 + m1*m2*l1^2*l2^2 + J2*m1*l1^2 + J1*m2*l2^2 + J1*J2
A = subs(A,J2*m2*L1^2 + m1*m2*l1^2*l2^2 + J2*m1*l1^2 + J1*m2*l2^2 + J1*J2,p);
B = subs(B,J2*m2*L1^2 + m1*m2*l1^2*l2^2 + J2*m1*l1^2 + J1*m2*l2^2 + J1*J2,p);

% A and B matrix of the linearized model
disp('A='); disp(A);
disp('B='); disp(B);

% And with the numeric values for the parameters
A_num = double(subs(A,params_sym, params_num));
B_num = double(subs(B,params_sym, params_num));
disp('A='); disp(A_num);
disp('B='); disp(B_num);