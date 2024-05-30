function [A,B] = getDoublePendulumAB(x,u)
%GETDOUBLEPENDULUMAB Summary of this function goes here
%   Detailed explanation goes here

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

% Get the angles
q1 = x(1);
q2 = x(2);

% and the angular velocities
q1_dot = x(3);
q2_dot = x(4);

p = 2*J1*J2 + L1^2*l2^2*m2^2 + 2*J2*L1^2*m2 + 2*J2*l1^2*m1 + 2*J1*l2^2*m2 - L1^2*l2^2*m2^2*cos(2*q1 - 2*q2) + 2*l1^2*l2^2*m1*m2;

A= [0,                                                                                                                                                                                                                                                                                                                                                     0,                                                                                                                           1,                                                                                            0;
    0,                                                                                                                                                                                                                                                                                                                                                     0,                                                                                                                           0,                                                                                            1;
    (L1*g*l2^2*m2^2*cos(q1) - 2*L1*l2^3*m2^2*q2_dot^2*cos(q1 - q2) + 2*J2*L1*g*m2*cos(q1) + L1*g*l2^2*m2^2*cos(q1 - 2*q2) + 2*J2*g*l1*m1*cos(q1) - 2*L1^2*l2^2*m2^2*q1_dot^2*cos(2*q1 - 2*q2) - 2*L1*l2*m2*q2_dot*r2*sin(q1 - q2) - 2*J2*L1*l2*m2*q2_dot^2*cos(q1 - q2) + 2*g*l1*l2^2*m1*m2*cos(q1))/p,                                                                                                                                                                                      (2*L1*l2*m2*(m2*cos(q1 - q2)*l2^2*q2_dot^2 + L1*m2*cos(2*q1 - 2*q2)*l2*q1_dot^2 - g*m2*cos(q1 - 2*q2)*l2 + J2*cos(q1 - q2)*q2_dot^2 + r2*sin(q1 - q2)*q2_dot))/p,                                                      -(2*q1_dot*sin(2*q1 - 2*q2)*L1^2*l2^2*m2^2 + 2*r1*l2^2*m2 + 2*J2*r1)/p, -(2*L1*l2*m2*(2*m2*q2_dot*sin(q1 - q2)*l2^2 - r2*cos(q1 - q2) + 2*J2*q2_dot*sin(q1 - q2)))/p;
    (2*L1*l2*m2*(m2*cos(q1 - q2)*L1^2*q1_dot^2 + l2*m2*cos(2*q1 - 2*q2)*L1*q2_dot^2 - g*m2*cos(2*q1 - q2)*L1 + m1*cos(q1 - q2)*l1^2*q1_dot^2 - g*m1*cos(2*q1 - q2)*l1 + J1*cos(q1 - q2)*q1_dot^2 - r1*sin(q1 - q2)*q1_dot + u*sin(q1 - q2)))/p, -(l2*m2*(2*L1*u*sin(q1 - q2) - 2*J1*g*cos(q2) - L1^2*g*m2*cos(2*q1 - q2) + 2*L1^3*m2*q1_dot^2*cos(q1 - q2) - 2*L1*q1_dot*r1*sin(q1 - q2) - L1^2*g*m2*cos(q2) + 2*J1*L1*q1_dot^2*cos(q1 - q2) - 2*g*l1^2*m1*cos(q2) + 2*L1^2*l2*m2*q2_dot^2*cos(2*q1 - 2*q2) - L1*g*l1*m1*cos(2*q1 - q2) + 2*L1*l1^2*m1*q1_dot^2*cos(q1 - q2) + L1*g*l1*m1*cos(q2)))/p, (2*L1*l2*m2*(2*m2*q1_dot*sin(q1 - q2)*L1^2 + 2*m1*q1_dot*sin(q1 - q2)*l1^2 + r1*cos(q1 - q2) + 2*J1*q1_dot*sin(q1 - q2)))/p,      -(- 2*q2_dot*sin(2*q1 - 2*q2)*L1^2*l2^2*m2^2 + 2*r2*L1^2*m2 + 2*m1*r2*l1^2 + 2*J1*r2)/p];
 
B = [0;
     0;
     (2*m2*l2^2 + 2*J2)/p;
    -(2*L1*l2*m2*cos(q1 - q2))/p];

end

