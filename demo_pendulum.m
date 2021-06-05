function demo_pendulum
% A demo of iLQG/DDP with Single pendulum dynamics
clc;
close all;

fprintf(['\nA demonstration of the iLQG algorithm '...
'with single pendulum dynamics.\n'...
'\"Control-Limited Differential Dynamic Programming\"\n'])

% Set full_DDP=true to compute 2nd order derivatives of the 
% dynamics. This will make iterations more expensive, but
% final convergence will be much faster (quadratic)
full_DDP = true;

% set up the optimization problem
DYNCST  = @(x,u,i) pendulum_dyn_cst(x,u,full_DDP);
T       = 500;              % horizon
x0      = [pi/2;0];            % initial state
u0      = 10*randn(1,T);    % initial controls
Op.lims  = [-10000 10000];  % external force force (N)
Op.plot = -1;               % plot the derivatives as well

% prepare the visualization window and graphics callback
figure(9);
set(gcf,'name','pendulum','Menu','none','NumberT','off')
set(gca,'xlim',[-10 10],'ylim',[-12 12],'DataAspectRatio',[1 1 1])
grid on
box on

% prepare and install trajectory visualization callback
line_handle = line([0 0],[0 0],'color','b','linewidth',2);
plotFn = @(x) set(line_handle,'Xdata',x(1,:),'Ydata',x(2,:));
Op.plotFn = plotFn;

% === run the optimization!
[x,u]= iLQG(DYNCST, x0, u0, Op);
save output_single_pendulum.mat x u

pendulum_plot(x);

function y = pendulum_dynamics(x,u)

% === states and controls:
% x = [phi dphi]' = [angle; angular velocity]
% u = [F]'     = [force applied]

% constants
l  = 0.5;      % l = length of pendulum
h  = 0.01;     % h = timestep (seconds)
m = 1;         % m = mass in Kg
b = 1;         % damping coefficient
g = 9.81;      % gravitational acceleration in m/s^2

% controls
F  = u(1,:,:);      % F = force applied tangetially

phi  = x(1,:,:);    % phi = angel w.r.t. verticle axis
dphi = x(2,:,:);    % dphi = angular velocity

z  = [u* 0 + dphi; F/(m*l) - b*dphi/(m*l) - g*sin(phi)/l]; %linear system

dy = h*z;           % change in state
y  = x + dy;        % new state


function c = pendulum_cost(x, u)
% cost function for single pendulum problem
% cost function is given by LQR 
% Q MATRIX
% R MATRIX

final = isnan(u(1,:));
u(:,final)  = 0;

Q = [1 0; 0 0.01];           %cost Q matrix
R = 0.00001;                 %cost R matrix

x_final = [260*pi/180;0];            %final state

n = size(u);
for i =1:max(n(:))
    c(i) = (x(:,i)-x_final)'*Q*(x(:,i)-x_final) + u(i)'*R*u(i);    %final cost
end


function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = pendulum_dyn_cst(x,u,full_DDP)
% combine pendulum dynamics and cost
% use helper function finite_difference() to compute derivatives

if nargout == 2
    f = pendulum_dynamics(x,u);
    c = pendulum_cost(x,u);
else
    % state and control indices
    ix = 1:2;
    iu = 3;
    
    % dynamics first derivatives
    xu_dyn  = @(xu) pendulum_dynamics(xu(ix,:),xu(iu,:));
    J       = finite_difference(xu_dyn, [x; u]);
    fx      = J(:,ix,:);
    fu      = J(:,iu,:);
    
    
    % dynamics second derivatives
    if full_DDP
        xu_Jcst = @(xu) finite_difference(xu_dyn, xu);
        JJ = finite_difference(xu_Jcst, [x; u]);
        JJ = reshape(JJ, [2 3 size(J,2) size(J,3)]);
        JJ = 0.5*(JJ + permute(JJ,[1 3 2 4])); %symmetrize
        fxx = JJ(:,ix,ix,:);
        fxu = JJ(:,ix,iu,:);
        fuu = JJ(:,iu,iu,:);
    else
        [fxx,fxu,fuu] = deal([]);
    end


% for single derivative   for full_DDP == false
%     if full_DDP
%         xu_Jcst = @(xu) finite_difference(xu_dyn, xu);
%         JJ      = finite_difference(xu_Jcst, [x; u]);
%         JJ      = reshape(JJ, [2 3 size(J)]);
%         JJ      = 0.5*(JJ + permute(JJ,[1 3 2 4])); %symmetrize
%         fxx     = JJ(:,ix,ix,:);
%         fxu     = JJ(:,ix,iu,:);
%         fuu     = JJ(:,iu,iu,:);
%     else
%         [fxx,fxu,fuu] = deal([]);
%     end
%     
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
function pendulum_plot(x)
% pedulum plot function animate the pendulum system

t=0:1:(length(x)-1);
X = 0.5*sin(x(1,:));
Y = -0.5*cos(x(1,:));

figure(10)
for ind = 1:length(t)
  plot([0,X(ind)],[0,Y(ind)], X(ind),Y(ind),'r.', 'MarkerSize',30) 
  axis([-0.6 0.6 -0.6 0.6])
  text(-0.2, 0.5,"Timer : " + num2str(t(ind)*0.01))
  grid on
  drawnow
  pause(0.1)
end

% utility functions: singleton-expanded addition and multiplication
function c = pp(a,b)
c = bsxfun(@plus,a,b);