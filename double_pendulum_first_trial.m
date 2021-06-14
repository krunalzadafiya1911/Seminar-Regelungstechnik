clc;
clear;
close all;

addpath('E:\casadi');
import casadi.*

T = 10; % Time horizon
N = 100; % number of control intervals
T1 = 100; % damping coefficient
T2 = 300;
m1 = 1; % mass
m2 = 2;
g = 9.81; % gravitational acceleration
l1 = 0.5; % length of pendulum
l2 = 0.3;

%final position of pendulum
x_final = [ pi; 0; pi; 0];

% Declare model variables
x1 = SX.sym('x1');
y1 = SX.sym('y1');
x2 = SX.sym('x2');
y2 = SX.sym('y2');

x1_dot = SX.sym('x1_dot');
y1_dot = SX.sym('y1_dot');
x2_dot = SX.sym('x2_dot');
y2_dot = SX.sym('y2_dot');

x = [x1; x1_dot; y1; y1_dot; x2; x2_dot; y2; y2_dot];

F1 = SX.sym('F1');
F2 = SX.sym('F2');

u = [F1; F2];
% Model equations
xdot = [x1_dot;
    1/m1*(F1- T1*x1/l1);
    y1_dot;
    1/m2*(T1-m1*g-T1*y1/l1);
    x2_dot;
    1/m2*(F2-T2*(x2-x1)/l2);
    y2_dot;
    1/m2*(T2-m2*g-T2*(x2-x1)/l2)];

% Objective term
Q = [1 0.5 1 0.5 1 0.5 1 0.5];
q = diag(Q);

R = [0.0001 0.00001];
r = diag(R);

L = x'*q*x + u'*r*u;
%L = 10*(x1-x_final(1))^2 + 10*(x2-x_final(2))^2 + (x3- x_final(3))^2 + 0.5*(x4- x_final(4))^2 + 1e-1*u1^2 + 1e-1*u2^2 ;

% Formulate discrete time dynamics
if false
   % CVODES from the SUNDIALS suite
   dae = struct('x',x,'p',u,'ode',xdot,'quad',L);
   opts = struct('tf',T/N);
   F = integrator('F', 'cvodes', dae, opts);
else
   % Fixed step Runge-Kutta 4 integrator
   M = 4; % RK4 steps per interval
   DT = T/N/M;
   f = Function('f', {x, u}, {xdot, L});
   X0 = MX.sym('X0', 8);
   U = MX.sym('U');
   X = X0;
   Q = 0;
   for j=1:M
       [k1, k1_q] = f(X, U);
       [k2, k2_q] = f(X + DT/2 * k1, U);
       [k3, k3_q] = f(X + DT/2 * k2, U);
       [k4, k4_q] = f(X + DT * k3, U);
       X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
       Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
    end
    F = Function('F', {X0, U}, {X, Q}, {'x0','p'}, {'xf', 'qf'});
end

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

% Formulate the NLP
Xk = [l1*sin(x_final(1)); 0; -l1*cos(x_final(1)); 0; l1*sin(x_final(1)) + l2* sin(x_final(3)); 0; -l1*cos(x_final(1))-l2*cos(x_final(3)); 0];

for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)]);
    w = {w{:}, Uk};
    lbw = [lbw, -1000];           %range for input para
    ubw = [ubw,  1000];
    w0 = [w0,  0];

    % Integrate till the end of the interval
    Fk = F('x0',Xk,'p', Uk);
    Xk = Fk.xf;
    J=J+Fk.qf;

    % Add inequality constraint
    g = {g{:}, Xk(1)};
    lbg = [lbg; -1000];
    ubg = [ubg;  1000];
end

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, ...
             'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

% Plot the solution
u_opt = w_opt;
x_opt = [pi/2;0; 0; 0];

x_opt1 = [l1*sin(x_opt(1)); 0; -l1*cos(x_opt(1)); 0; l1*sin(x_opt(1)) + l2* sin(x_opt(3)); 0; -l1*cos(x_opt(1))-l2*cos(x_opt(3)); 0];
for k=0:N-1
    Fk = F('x0', x_opt1(:,end), 'p', u_opt(k+1));
    x_opt1 = [x_opt1, full(Fk.xf)];
end
x1_opt = x_opt1(1,:);
x2_opt = x_opt1(2,:);
tgrid = linspace(0, T, N+1);
clf;
hold on
plot(tgrid, x1_opt, '--')
plot(tgrid, x2_opt, '-')
stairs(tgrid, [u_opt; nan], '-.')
xlabel('t')
legend('x1','x2','u')

%%
t=0:1:(length(x1_opt)-1);
X1 = x_opt1(1,:);
Y1 = x_opt1(3,:);
X2 = x_opt1(5,:);
Y2 = x_opt1(7,:);

figure(2)
for ind = 1:length(t)
  plot([0,X1(ind)],[0,Y1(ind)], [X1(ind) , X2(ind)],[Y1(ind), Y2(ind)], X1(ind),Y1(ind),'r.', X2(ind),Y2(ind),'b.','MarkerSize',30)
  %plot([X1(ind),X2(ind)],[Y1(ind),Y2(ind)], X2(ind),Y2(ind),'b.', 'MarkerSize',50)
  axis([-l1-l2-0.1 l1++l2+0.1 -l1-l2-0.1 l1+l2+0.1])
  text(-0.2, 0.5,"Timer : " + num2str(t(ind)*0.01))
  grid on
  drawnow
  pause(0.1)
end