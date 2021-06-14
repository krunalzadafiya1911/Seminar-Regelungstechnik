clc;
clear;
close all;

addpath('E:\casadi');
import casadi.*

T = 10; % Time horizon
N = 100; % number of control intervals
b = 1; % damping coefficient
m = 1; % mass
g = 9.81; % gravitational acceleration
l = 0.5; % length of pendulum

%final position of pendulum
x_final = [ pi/2; 0];

% Declare model variables
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1; x2];
u = SX.sym('u');

% Model equations
xdot = [x2; u/(m*l) - b*x2/(m*l) - g*sin(x1)/l];

% Objective term
L = (x1)^2 + 0.1*(x2)^2 + 1e-6*u^2;

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
   X0 = MX.sym('X0', 2);
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
Xk = [pi; 0];

for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)]);
    w = {w{:}, Uk};
    lbw = [lbw, -10];           %range for input para
    ubw = [ubw,  10];
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
x_opt = [pi;0];
for k=0:N-1
    Fk = F('x0', x_opt(:,end), 'p', u_opt(k+1));
    x_opt = [x_opt, full(Fk.xf)];
end
x1_opt = x_opt(1,:);
x2_opt = x_opt(2,:);
tgrid = linspace(0, T, N+1);
clf;
hold on
plot(tgrid, x1_opt, '--')
plot(tgrid, x2_opt, '-')
stairs(tgrid, [u_opt; nan], '-.')
xlabel('t')
legend('x1','x2','u')


t=0:1:(length(x1_opt)-1);
X = l*sin(x_opt(1,:));
Y = -l*cos(x_opt(1,:));

figure(2)
for ind = 1:length(t)
  plot([0,X(ind)],[0,Y(ind)], X(ind),Y(ind),'r.', 'MarkerSize',30) 
  axis([-l-0.1 l+0.1 -l-0.1 l+0.1])
  text(-0.2, 0.5,"Timer : " + num2str(t(ind)*0.01))
  grid on
  drawnow
  pause(0.1)
end