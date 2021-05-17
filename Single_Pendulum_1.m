%% Seminar Control

%Problem: stabilize the single pendulum using LQR

clc;
clear;
clear all;

%% state space model repersentation

% Initial conditions
X0 = [pi;0];

% system installation
syms m l g phi(t) f(t) gamma b

dphi(t) = diff(phi(t));         %Angle velocity
ddphi(t) = diff(dphi(t));       %Angular accekeration

eqn = m*l*ddphi(t) + m*g*sin(phi) + b*dphi == f(t);

eq_isolateDDPhi = isolate(eqn,ddphi(t))

eq_expand = expand(eq_isolateDDPhi);

rhs_eq = rhs(eq_expand);

syms x_1 x_2 u

subphi = phi(t)+gamma;
fsub = u;
gsub = 9.81;
lsub = 0.5;
msub = 1;

replace_1 = [m,l,g,f,phi];
replaced_by = [msub,lsub,gsub,fsub,subphi];
eq_subs = subs(eq_expand,replace_1,replaced_by)

gammasub = 0;
bsub = 0;
eq_subs1 = subs(eq_subs,[gamma,b],[gammasub,bsub]);

%system of state variables
subphi = x_1;
subdphi = x_2;

clear u 
syms u
f(t) = u;           %controller input
 
y(t) = x_1;         %function output
 
xd_1 = x_2;         %state space equations
xd_2 = subs(rhs(eq_subs),[phi,dphi],[subphi,subdphi]);

X = [x_1; x_2];     %state vector

F = [xd_1; xd_2];   %x'

% System dynamics A,B,C,D
A_1 = jacobian(F, X);

xa_1 = X0(1);
xa_2 = X0(2);
damping = 1;
gamma_a = 0;
A = subs(A_1,[gamma x_1 x_2 b], [gamma_a,xa_1,xa_2,damping]);
B = jacobian(F, u);
C = jacobian(y, X);
D = jacobian(y, u);

A = [0 1; 981/50 -2];
B = [0;2];
C = [1 0];
D = 0;
% Control parameters
Q = [1 0;0 1];
R = 1;
K = lqr(A, B, Q, R);

% Closed Loop system
sys = ss((A - B*K), B, C, D);

% Run response to intialize the condition
t = 0:0.005:30;
[y,t,x] = initial(sys, X0, t);

plot(t,y)
 %% animation of single double pendulum

% sol.x = x;
% sol.y = y;
% 
% x_1 = @(t) lsub*sin(sol.y);
% y_1 = @(t) -lsub*cos(sol.y);
% %x_2 = @(t) L_1*sin(deval(sols,t,3))+L_2*sin(deval(sols,t,1));
% %y_2 = @(t) -L_1*cos(deval(sols,t,3))-L_2*cos(deval(sols,t,1));
% 
% fanimator(@(t) plot(x_1(t),y_1(t),'ro','MarkerSize',msub*5,'MarkerFaceColor','r'));
% axis equal;
% 
% hold on;
% fanimator(@(t) plot([0 x_1(t)],[0 y_1(t)],'r-'));
% %fanimator(@(t) plot(x_2(t),y_2(t),'go','MarkerSize',m_2*10,'MarkerFaceColor','g'));
% %fanimator(@(t) plot([x_1(t) x_2(t)],[y_1(t) y_2(t)],'g-'));
% 
% fanimator(@(t) text(-0.3,0.3,"Timer: "+num2str(t,2)));
% hold off;
% 
% 
% playAnimation