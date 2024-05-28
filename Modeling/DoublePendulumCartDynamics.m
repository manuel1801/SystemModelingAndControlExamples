% This function implements the nonlinear continuous-time dynamics of the
% inverted double pendulum on a cart xdot(t) = f( x(t), u(t) ), derived using
% DoublePendulumCartLagrange.m Given the current state x and the current input u, 
% DoublePendulumCartDynamics(x,u) computes the current time derivative xdot of the
% pendulum's states. 

function xdot = DoublePendulumCartDynamics(x,u)

    % System paramter
    m0 = 1.5;
    m1 = 0.75;
    m2 = 1;
    g = 9.81;
    L1 = 0.75;
    L2 = 1;     
    l1 = 1/2*L1;
    l2 = 1/2*L2;
    J1 = m1*L1^2/12;
    J2 = m2*L2^2/12;
    r1 = 0.02; r2 = 0.04;

    % Get the angles
    q0 = x(1);
    q1 = x(2);
    q2 = x(3);

    % and the angular velocities
    q0_dot = x(4);
    q1_dot = x(5);
    q2_dot = x(6);

    p = (J2*L1^2*m2^2 + J2*l1^2*m1^2 + J1*l2^2*m2^2 + 2*J1*J2*m0 + 2*J1*J2*m1 + 2*J1*J2*m2 + L1^2*l2^2*m0*m2^2 + L1^2*l2^2*m1*m2^2 + 2*J2*L1^2*m0*m2 + 2*J2*L1^2*m1*m2 + l1^2*l2^2*m1*m2^2 + l1^2*l2^2*m1^2*m2 + 2*J2*l1^2*m0*m1 + 2*J1*l2^2*m0*m2 + 2*J1*l2^2*m1*m2 + 2*J2*l1^2*m1*m2 - J2*L1^2*m2^2*cos(2*q1) - J2*l1^2*m1^2*cos(2*q1) - J1*l2^2*m2^2*cos(2*q2) - L1*l1*l2^2*m1*m2^2 - 2*J2*L1*l1*m1*m2 + 2*l1^2*l2^2*m0*m1*m2 - l1^2*l2^2*m1^2*m2*cos(2*q1) - l1^2*l2^2*m1*m2^2*cos(2*q2) - L1^2*l2^2*m0*m2^2*cos(2*q1 - 2*q2) - L1^2*l2^2*m1*m2^2*cos(2*q1 - 2*q2) - L1*l1*l2^2*m1*m2^2*cos(2*q1) + L1*l1*l2^2*m1*m2^2*cos(2*q2) - 2*J2*L1*l1*m1*m2*cos(2*q1) + L1*l1*l2^2*m1*m2^2*cos(2*q1 - 2*q2));


    xdot = [q0_dot;
            q1_dot;
            q2_dot;
            (2*J1*J2*u + L1^2*l2^2*m2^2*u + 2*J2*L1^2*m2*u + 2*J2*l1^2*m1*u + 2*J1*l2^2*m2*u + 2*l1^2*l2^2*m1*m2*u - L1^2*l2^2*m2^2*u*cos(2*q1 - 2*q2) - J2*L1^2*g*m2^2*sin(2*q1) + 2*J2*L1^3*m2^2*q1_dot^2*sin(q1) - J2*g*l1^2*m1^2*sin(2*q1) - J1*g*l2^2*m2^2*sin(2*q2) + 2*J2*l1^3*m1^2*q1_dot^2*sin(q1) + 2*J1*l2^3*m2^2*q2_dot^2*sin(q2) + 2*J2*L1*m2*q1_dot*r1*cos(q1) - L1*l2^2*m2^2*q1_dot*r1*cos(q1 - 2*q2) + 2*J2*l1*m1*q1_dot*r1*cos(q1) + 2*J1*l2*m2*q2_dot*r2*cos(q2) + J1*L1*l2^2*m2^2*q1_dot^2*sin(q1) + J2*L1^2*l2*m2^2*q2_dot^2*sin(q2) + 2*J1*J2*L1*m2*q1_dot^2*sin(q1) - J1*L1*l2^2*m2^2*q1_dot^2*sin(q1 - 2*q2) + 2*J1*J2*l1*m1*q1_dot^2*sin(q1) + 2*J1*J2*l2*m2*q2_dot^2*sin(q2) - L1^2*l2*m2^2*q2_dot*r2*cos(2*q1 - q2) + J2*L1^2*l2*m2^2*q2_dot^2*sin(2*q1 - q2) - g*l1^2*l2^2*m1^2*m2*sin(2*q1) - g*l1^2*l2^2*m1*m2^2*sin(2*q2) + 2*l1^3*l2^2*m1^2*m2*q1_dot^2*sin(q1) + 2*l1^2*l2^3*m1*m2^2*q2_dot^2*sin(q2) + L1*l2^2*m2^2*q1_dot*r1*cos(q1) + L1^2*l2*m2^2*q2_dot*r2*cos(q2) + 2*l1*l2^2*m1*m2*q1_dot*r1*cos(q1) + 2*l1^2*l2*m1*m2*q2_dot*r2*cos(q2) + L1*l1^2*l2^2*m1*m2^2*q1_dot^2*sin(q1) + L1^2*l1*l2^2*m1*m2^2*q1_dot^2*sin(q1) + 2*J2*L1*l1^2*m1*m2*q1_dot^2*sin(q1) + 2*J2*L1^2*l1*m1*m2*q1_dot^2*sin(q1) + L1*l1*l2^3*m1*m2^2*q2_dot^2*sin(2*q1 - q2) - L1*l1^2*l2^2*m1*m2^2*q1_dot^2*sin(q1 - 2*q2) + L1^2*l1*l2^2*m1*m2^2*q1_dot^2*sin(q1 - 2*q2) + 2*J1*l1*l2^2*m1*m2*q1_dot^2*sin(q1) + 2*J2*l1^2*l2*m1*m2*q2_dot^2*sin(q2) - L1*g*l1*l2^2*m1*m2^2*sin(2*q1) + L1*g*l1*l2^2*m1*m2^2*sin(2*q2) - L1*l1*l2^3*m1*m2^2*q2_dot^2*sin(q2) - 2*J2*L1*g*l1*m1*m2*sin(2*q1) - L1*l1*l2*m1*m2*q2_dot*r2*cos(q2) - J2*L1*l1*l2*m1*m2*q2_dot^2*sin(q2) - L1*l1*l2*m1*m2*q2_dot*r2*cos(2*q1 - q2) + J2*L1*l1*l2*m1*m2*q2_dot^2*sin(2*q1 - q2))/p;
            -(l2^2*m2^2*q1_dot*r1 + 2*J2*m0*q1_dot*r1 + 2*J2*m1*q1_dot*r1 + 2*J2*m2*q1_dot*r1 - 2*J2*L1*g*m2^2*sin(q1) - l2^2*m2^2*q1_dot*r1*cos(2*q2) - 2*J2*g*l1*m1^2*sin(q1) + J2*L1^2*m2^2*q1_dot^2*sin(2*q1) + J2*l1^2*m1^2*q1_dot^2*sin(2*q1) + L1*l2^2*m2^2*u*cos(q1) + 2*J2*L1*m2*u*cos(q1) - L1*l2^2*m2^2*u*cos(q1 - 2*q2) + 2*J2*l1*m1*u*cos(q1) + 2*l2^2*m0*m2*q1_dot*r1 + 2*l2^2*m1*m2*q1_dot*r1 + J2*L1*l2*m2^2*q2_dot^2*sin(q1 - q2) - 2*J2*L1*g*m0*m2*sin(q1) - 2*J2*L1*g*m1*m2*sin(q1) - g*l1*l2^2*m1*m2^2*sin(q1) - 2*g*l1*l2^2*m1^2*m2*sin(q1) - L1*g*l2^2*m0*m2^2*sin(q1 - 2*q2) - L1*g*l2^2*m1*m2^2*sin(q1 - 2*q2) - 2*J2*g*l1*m0*m1*sin(q1) - 2*J2*g*l1*m1*m2*sin(q1) + l1^2*l2^2*m1^2*m2*q1_dot^2*sin(2*q1) + g*l1*l2^2*m1*m2^2*sin(q1 - 2*q2) + L1^2*l2^2*m0*m2^2*q1_dot^2*sin(2*q1 - 2*q2) + L1^2*l2^2*m1*m2^2*q1_dot^2*sin(2*q1 - 2*q2) + l1*l2^3*m1*m2^2*q2_dot^2*sin(q1 + q2) + L1*l2*m2^2*q2_dot*r2*cos(q1 + q2) + 2*L1*l2^3*m0*m2^2*q2_dot^2*sin(q1 - q2) + 2*L1*l2^3*m1*m2^2*q2_dot^2*sin(q1 - q2) + J2*L1*l2*m2^2*q2_dot^2*sin(q1 + q2) - l1*l2^3*m1*m2^2*q2_dot^2*sin(q1 - q2) + 2*l1*l2^2*m1*m2*u*cos(q1) - L1*l2*m2^2*q2_dot*r2*cos(q1 - q2) - L1*g*l2^2*m0*m2^2*sin(q1) - L1*g*l2^2*m1*m2^2*sin(q1) + 2*J2*L1*l2*m0*m2*q2_dot^2*sin(q1 - q2) + 2*J2*L1*l2*m1*m2*q2_dot^2*sin(q1 - q2) - 2*g*l1*l2^2*m0*m1*m2*sin(q1) + L1*l1*l2^2*m1*m2^2*q1_dot^2*sin(2*q1) - J2*l1*l2*m1*m2*q2_dot^2*sin(q1 - q2) + 2*J2*L1*l1*m1*m2*q1_dot^2*sin(2*q1) - L1*l1*l2^2*m1*m2^2*q1_dot^2*sin(2*q1 - 2*q2) + l1*l2*m1*m2*q2_dot*r2*cos(q1 + q2) - 2*L1*l2*m0*m2*q2_dot*r2*cos(q1 - q2) - 2*L1*l2*m1*m2*q2_dot*r2*cos(q1 - q2) + J2*l1*l2*m1*m2*q2_dot^2*sin(q1 + q2) + l1*l2*m1*m2*q2_dot*r2*cos(q1 - q2))/p;
            (l1^2*m1^2*q2_dot*r2*cos(2*q1) - l1^2*m1^2*q2_dot*r2 - 2*J1*m0*q2_dot*r2 - 2*J1*m1*q2_dot*r2 - 2*J1*m2*q2_dot*r2 - L1^2*m2^2*q2_dot*r2 + L1^2*l2*m2^2*u*cos(2*q1 - q2) + 2*J1*g*l2*m2^2*sin(q2) - J1*l2^2*m2^2*q2_dot^2*sin(2*q2) - L1^2*l2*m2^2*u*cos(q2) - 2*L1^2*m0*m2*q2_dot*r2 - 2*L1^2*m1*m2*q2_dot*r2 - 2*J1*l2*m2*u*cos(q2) - 2*l1^2*m0*m1*q2_dot*r2 - 2*l1^2*m1*m2*q2_dot*r2 + L1^2*m2^2*q2_dot*r2*cos(2*q1) + J1*L1*l2*m2^2*q1_dot^2*sin(q1 - q2) + 2*g*l1^2*l2*m1*m2^2*sin(q2) + g*l1^2*l2*m1^2*m2*sin(q2) + 2*J1*g*l2*m0*m2*sin(q2) + 2*J1*g*l2*m1*m2*sin(q2) - l1^2*l2^2*m1*m2^2*q2_dot^2*sin(2*q2) + L1^2*l2^2*m0*m2^2*q2_dot^2*sin(2*q1 - 2*q2) + L1^2*l2^2*m1*m2^2*q2_dot^2*sin(2*q1 - 2*q2) - l1^3*l2*m1^2*m2*q1_dot^2*sin(q1 + q2) - L1*l2*m2^2*q1_dot*r1*cos(q1 + q2) - L1^2*g*l2*m0*m2^2*sin(2*q1 - q2) - L1^2*g*l2*m1*m2^2*sin(2*q1 - q2) + 2*L1^3*l2*m0*m2^2*q1_dot^2*sin(q1 - q2) + 2*L1^3*l2*m1*m2^2*q1_dot^2*sin(q1 - q2) + g*l1^2*l2*m1^2*m2*sin(2*q1 - q2) - J1*L1*l2*m2^2*q1_dot^2*sin(q1 + q2) - l1^3*l2*m1^2*m2*q1_dot^2*sin(q1 - q2) - 2*l1^2*l2*m1*m2*u*cos(q2) + L1*l2*m2^2*q1_dot*r1*cos(q1 - q2) + 2*L1*l1*m1*m2*q2_dot*r2 + L1^2*g*l2*m0*m2^2*sin(q2) + L1^2*g*l2*m1*m2^2*sin(q2) + L1*l1^2*l2*m1*m2^2*q1_dot^2*sin(q1 - q2) + L1*l1^2*l2*m1^2*m2*q1_dot^2*sin(q1 - q2) - 3*L1^2*l1*l2*m1*m2^2*q1_dot^2*sin(q1 - q2) - 3*L1*g*l1*l2*m1*m2^2*sin(q2) - L1*g*l1*l2*m1^2*m2*sin(q2) + 2*L1*l1*m1*m2*q2_dot*r2*cos(2*q1) + 2*J1*L1*l2*m0*m2*q1_dot^2*sin(q1 - q2) + 2*J1*L1*l2*m1*m2*q1_dot^2*sin(q1 - q2) + 2*g*l1^2*l2*m0*m1*m2*sin(q2) + L1*l1*l2^2*m1*m2^2*q2_dot^2*sin(2*q2) + L1*l1*l2*m1*m2*u*cos(2*q1 - q2) - J1*l1*l2*m1*m2*q1_dot^2*sin(q1 - q2) - L1*l1*l2^2*m1*m2^2*q2_dot^2*sin(2*q1 - 2*q2) + L1*g*l1*l2*m1*m2^2*sin(2*q1 - q2) - L1*g*l1*l2*m1^2*m2*sin(2*q1 - q2) - l1*l2*m1*m2*q1_dot*r1*cos(q1 + q2) - L1*l1^2*l2*m1*m2^2*q1_dot^2*sin(q1 + q2) + L1*l1^2*l2*m1^2*m2*q1_dot^2*sin(q1 + q2) + L1^2*l1*l2*m1*m2^2*q1_dot^2*sin(q1 + q2) + L1*l1*l2*m1*m2*u*cos(q2) + 2*L1*l2*m0*m2*q1_dot*r1*cos(q1 - q2) + 2*L1*l2*m1*m2*q1_dot*r1*cos(q1 - q2) - J1*l1*l2*m1*m2*q1_dot^2*sin(q1 + q2) - l1*l2*m1*m2*q1_dot*r1*cos(q1 - q2) + 2*L1*l1^2*l2*m0*m1*m2*q1_dot^2*sin(q1 - q2) - L1*g*l1*l2*m0*m1*m2*sin(q2) - L1*g*l1*l2*m0*m1*m2*sin(2*q1 - q2))/p];



end