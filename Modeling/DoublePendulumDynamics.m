% This function implements the nonlinear continuous-time dynamics of the
% inverted double pendulum xdot(t) = f( x(t), u(t) ), derived using
% DoublePendulumLagrange.m Given the current state x and the current input u, 
% DoublePendulumDynamics(x,u) computes the current time derivative xdot of the inverted double 
% pendulum's states. 

function xdot = DoublePendulumDynamics(x,u)

    % System paramter
    m1 = 0.75;
    m2 = 1;
    g = 9.81;
    L1 = 0.75;
    L2 = 1;     
    l1 = 1/2*L1;
    l2 = 1/2*L2;
    J1 = m1*L1^2/12;
    J2 = m2*L2^2/12;
    r1 = 0.5; r2 = 0.5;

    % Get the angles
    q1 = x(1);
    q2 = x(2);

    % and the angular velocities
    q1_dot = x(3);
    q2_dot = x(4);

    % Evaluate the nonlinear dynamics function
    xdot = [q1_dot
            q2_dot
            (2*J2*u + 2*l2^2*m2*u - 2*J2*q1_dot*r1 - 2*l2^2*m2*q1_dot*r1 - 2*L1*l2^3*m2^2*q2_dot^2*sin(q1 - q2) + L1*g*l2^2*m2^2*sin(q1) + 2*J2*L1*g*m2*sin(q1) + L1*g*l2^2*m2^2*sin(q1 - 2*q2) + 2*J2*g*l1*m1*sin(q1) - L1^2*l2^2*m2^2*q1_dot^2*sin(2*q1 - 2*q2) + 2*L1*l2*m2*q2_dot*r2*cos(q1 - q2) - 2*J2*L1*l2*m2*q2_dot^2*sin(q1 - q2) + 2*g*l1*l2^2*m1*m2*sin(q1))/(2*J1*J2 + L1^2*l2^2*m2^2 + 2*J2*L1^2*m2 + 2*J2*l1^2*m1 + 2*J1*l2^2*m2 - L1^2*l2^2*m2^2*cos(2*q1 - 2*q2) + 2*l1^2*l2^2*m1*m2)
            (2*L1^3*l2*m2^2*q1_dot^2*sin(q1 - q2) - 2*L1^2*m2*q2_dot*r2 - 2*l1^2*m1*q2_dot*r2 - L1^2*g*l2*m2^2*sin(2*q1 - q2) - 2*J1*q2_dot*r2 + L1^2*g*l2*m2^2*sin(q2) + 2*J1*g*l2*m2*sin(q2) + L1^2*l2^2*m2^2*q2_dot^2*sin(2*q1 - 2*q2) - 2*L1*l2*m2*u*cos(q1 - q2) + 2*L1*l2*m2*q1_dot*r1*cos(q1 - q2) + 2*J1*L1*l2*m2*q1_dot^2*sin(q1 - q2) + 2*g*l1^2*l2*m1*m2*sin(q2) - L1*g*l1*l2*m1*m2*sin(2*q1 - q2) + 2*L1*l1^2*l2*m1*m2*q1_dot^2*sin(q1 - q2) - L1*g*l1*l2*m1*m2*sin(q2))/(2*J1*J2 + L1^2*l2^2*m2^2 + 2*J2*L1^2*m2 + 2*J2*l1^2*m1 + 2*J1*l2^2*m2 - L1^2*l2^2*m2^2*cos(2*q1 - 2*q2) + 2*l1^2*l2^2*m1*m2)];

end

