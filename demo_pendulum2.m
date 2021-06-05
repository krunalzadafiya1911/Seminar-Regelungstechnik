function demo_pendulum2
% A demo of iLQG/DDP with Single pendulum dynamics
clc;
close all;

fprintf(['\nA demonstration of the iLQG algorithm '...
'with single pendulum dynamics.\n'...
'\"Control-Limited Differential Dynamic Programming\"\n'])

% Set full_DDP=true to compute 2nd order derivatives of the 
% dynamics. This will make iterations more expensive, but
% final convergence will be much faster (quadratic)
full_DDP = false;

% set up the optimization problem
DYNCST  = @(x,u,i) pendulum_dyn_cst(x,u,full_DDP);
T       = 500;              % horizon
x0      = [0;0;0;0];            % initial state
u0      = 10*randn(2,T);    % initial controls
Op.lims  = [-1000000 1000000];  % external force force (N)
Op.plot = -1;               % plot the derivatives as well

% prepare the visualization window and graphics callback
figure(9);
set(gcf,'name','pendulum','Menu','none','NumberT','off')
set(gca,'xlim',[-10 10],'ylim',[-12 12],'DataAspectRatio',[1 1 1])
grid on
box on

% prepare and install trajectory visualization callback
line_handle = line([0 0],[0 0],'color','b','linewidth',2);
plotFn = @(x) set(line_handle,'Xdata',x(3,:),'Ydata',x(4,:));
Op.plotFn = plotFn;

% === run the optimization!
[x,u]= iLQG(DYNCST, x0, u0, Op);
save output_2_pendulum.mat x u

pendulum_plot(x, u);

function y = pendulum_dynamics(x,u)

% === states and controls:
% x = [phi dphi]' = [angle; angular velocity]
% u = [F]'     = [force applied]

% constants
L_1  = 0.5;      % l = length of pendulum
L_2 = 0.3;
h  = 0.01;     % h = timestep (seconds)
m_1 = 1;         % m = mass in Kg
m_2 = 2;         % m2 = mass
%b = 1;         % damping coefficient

g = 9.81;      % gravitational acceleration in m/s^2

% controls
Q_1  = u(1,:,:);      % F = force applied tangetially
Q_2 = u(2,:,:); 

theta_1  = x(1,:,:);    % phi = angel w.r.t. verticle axis
theta_1_dot = x(2,:,:);   
theta_2 = x(3,:,:);
theta_2_dot = x(4,:,:);

theta_1_dot_dot = (- L_1*L_2*m_2*sin(theta_1 - theta_2).*theta_2_dot.^2 + Q_1 + L_1*sin(theta_1).*((g*m_1) + (981*m_2)/100) - (L_1*cos(theta_1 - theta_2).*(m_2*sin(theta_1 - theta_2)*L_2^2.*theta_1_dot.^2 - (g*m_2*sin(theta_2)*L_2) + Q_2))/L_2)./(L_1^2*(m_1 + m_2) - L_1^2*m_2*cos(theta_1 - theta_2).^2);
theta_2_dot_dot = (Q_2 - (g*L_2*m_2*sin(theta_2)) + L_2^2*m_2*(theta_1_dot.^2).*sin(theta_1 - theta_2) - (L_2*m_2*cos(theta_1 - theta_2).*(- L_1*L_2*m_2*sin(theta_1 - theta_2).*theta_2_dot.^2 + Q_1 + L_1*sin(theta_1).*((g*m_1) + (g*m_2))))/(L_1*(m_1 + m_2)))./(L_2^2*m_2 - (L_2^2*m_2^2*cos(theta_1 - theta_2).^2)/(m_1 + m_2));
 

z  = [Q_1* 0 + theta_1_dot;
    theta_1_dot_dot;
    Q_2*0 + theta_2_dot;
    theta_2_dot_dot];

dy = h*z;           % change in state
y  = x + dy;        % new state


function c = pendulum_cost(x, u)
% cost function for single pendulum problem
% cost function is given by LQR 
% Q MATRIX
% R MATRIX

final = isnan(u(1,:));
u(:,final)  = 0;

Q = diag([100 1 100 1]);           %cost Q matrix
R = diag([1e-2 1e-2]);                 %cost R matrix

x_final = [pi;0;pi;0];            %final state

n = size(u);
for i =1:n(2)
    c(i) = (x(:,i)-x_final)'*Q*(x(:,i)-x_final) + u(:,i)'*R*u(:,i);    %final cost
end


function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = pendulum_dyn_cst(x,u,full_DDP)
% combine pendulum dynamics and cost
% use helper function finite_difference() to compute derivatives

if nargout == 2
    f = pendulum_dynamics(x,u);
    c = pendulum_cost(x,u);
else
    % state and control indices
    ix = 1:4;
    iu = 5:6;
    
    % dynamics first derivatives
    xu_dyn  = @(xu) pendulum_dynamics(xu(ix,:),xu(iu,:));
    J       = finite_difference(xu_dyn, [x; u]);
    fx      = J(:,ix,:);
    fu      = J(:,iu,:);
    
    
%     % dynamics second derivatives
%     if full_DDP
%         xu_Jcst = @(xu) finite_difference(xu_dyn, xu);
%         JJ = finite_difference(xu_Jcst, [x; u]);
%         JJ = reshape(JJ, [4 6 size(J,2) size(J,3)]);
%         JJ = 0.5*(JJ + permute(JJ,[1 3 2 4])); %symmetrize
%         fxx = JJ(:,ix,ix,:);
%         fxu = JJ(:,ix,iu,:);
%         fuu = JJ(:,iu,iu,:);
%     else
%         [fxx,fxu,fuu] = deal([]);
%     end


%for single derivative   for full_DDP == false
    if full_DDP
        xu_Jcst = @(xu) finite_difference(xu_dyn, xu);
        JJ      = finite_difference(xu_Jcst, [x; u]);
        JJ      = reshape(JJ, [2 3 size(J)]);
        JJ      = 0.5*(JJ + permute(JJ,[1 3 2 4])); %symmetrize
        fxx     = JJ(:,ix,ix,:);
        fxu     = JJ(:,ix,iu,:);
        fuu     = JJ(:,iu,iu,:);
    else
        [fxx,fxu,fuu] = deal([]);
    end
    
    % cost first derivatives
    xu_cost = @(xu) pendulum_cost(xu(ix,:),xu(iu,:));
    J       = squeeze(finite_difference(xu_cost, [x; u]));
    cx      = J(ix,:);
    cu      = J(iu,:);
    
    % cost second derivatives
    xu_Jcst = @(xu) squeeze(finite_difference(xu_cost, xu));
    JJ      = finite_difference(xu_Jcst, [x; u]);
    JJ      = 0.5*(JJ + permute(JJ,[2 1 3]));      %symmetrize
    cxx     = JJ(ix,ix,:);
    cxu     = JJ(ix,iu,:);
    cuu     = JJ(iu,iu,:);
    
    [f,c] = deal([]);
end

function J = finite_difference(fun, x, h)
% simple finite-difference derivatives
% assumes the function fun() is vectorized

if nargin < 3
    h = 2^-17;
end

[n, K]  = size(x);
H       = [zeros(n,1) h*eye(n)];
H       = permute(H, [1 3 2]);
X       = pp(x, H);
X       = reshape(X, n, K*(n+1));
Y       = fun(X);
m       = numel(Y)/(K*(n+1));
Y       = reshape(Y, m, K, n+1);
J       = pp(Y(:,:,2:end), -Y(:,:,1)) / h;
J       = permute(J, [1 3 2]);


% ======== graphics functions ========
function pendulum_plot(x, u)
% pedulum plot function animate the pendulum system

L_1 = 0.5;
L_2 = 0.3;

t=0:1:(length(x)-1);
X1 = L_1*sin(x(1,:));
Y1 = -L_1*cos(x(1,:));
X2 = L_1*sin(x(1,:)) + L_2*sin(x(3,:));
Y2 = -L_1*cos(x(1,:)) - L_2*cos(x(3,:));


theta_1 = x(1,:);
theta_1_dot = x(2,:);
theta_2 = x(3,:);
theta_2_dot = x(4,:);

tgrid = linspace(0, max(size(u)), 501);
figure(9)
clf;
subplot(1,2,1)
plot(tgrid, theta_1, '--')
hold on
plot(tgrid, theta_1_dot, '-')
hold on
stairs(tgrid, [u(1,:)'; nan], '-.')
xlabel('t')
legend('theta_1','theta_1_dot','u1')

subplot(2,2,1)
plot(tgrid, theta_2, '--')
hold on
plot(tgrid, theta_2_dot, '-')
hold on
stairs(tgrid, [u(2,:)'; nan], '-.')
hold on
xlabel('t')
legend('theta_2','theta_2_dot','u2')

figure(10)
for ind = 1:length(t)
  plot([0,X1(ind)],[0,Y1(ind)],[X1(ind), X2(ind)],[Y1(ind), Y2(ind)], X1(ind),Y1(ind),'r.', X2(ind),Y2(ind),'b.', 'MarkerSize',50) 
  axis([-L_1-L_2-0.1 L_1+L_2+0.3 -L_1-L_2-0.1 L_1+L_2+0.1])
  text(-0.1, L_1+L_2+0.2, "Timer : " + num2str(t(ind)*0.01))
  grid on
  drawnow
  pause(0.1)
end

% utility functions: singleton-expanded addition and multiplication
function c = pp(a,b)
c = bsxfun(@plus,a,b);