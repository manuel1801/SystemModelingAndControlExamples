function [A,B] = getDoublePendulumCartAB(x,u)
%GETDOUBLEPENDULUMAB Summary of this function goes here
%   Detailed explanation goes here


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

% Get the angles
q1 = x(1);
q2 = x(2);

% and the angular velocities
q1_dot = x(3);
q2_dot = x(4);



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


end

