addpath('E:\casadi');
import casadi.*

T = 0.2;            %sampling time
N = 3;              %prediction horizon
rob_diam = 0.3;

v_max = 0.6;
v_min = -v_max;
omega_max = pi/4;
omega_min = -omega_max;

%define the state variables 
x = SX.sym('x');
x_dot = SX.sym('x_dot');
y = SX.sym('y');
y_dot = SX.sym('y_dot');

states = [x; x_dot; y; y_dot];
n_states = length(states);

%define the input parameters
Fx = SX.sym('Fx');
Fy = SX.sym('Fy');
%omega= SX.sym('omega');

controls = [Fx, Fy];
n_controls = length(controls);
rhs = [v*cos(theta);v*sin(theta);omega];    %system

f = Function('f',{states, controls}, {rhs});    %nonliner mapping function
U = SX.sym('U', n_controls, N);                 %decision variables (controls)
P = SX.sym('P', n_states+ n_states);

X= SX.sym('X', n_states,(N+1));

%compute system symbolically
X(:,1) = P(1:3);        %intial state

for k = 1:N
    st = X(:,k);
    con = U(:,k);
    f_value = f(st,con);
    st_next = st + (T*f_value);
    X(:,k+1) = st_next;
end


%this function is to get the optimal tranjectory

ff = Function('ff',{U,P}, {X});

obj = 0;            %obejctive function
g = [];             %constraints vector

Q = zeros(3,3); Q(1,1) = 1; Q(2,2) = 5; Q(3,3) = 0.1; % weighting matrix (state)
R = zeros(2,2); R(1,1) = 0.5; R(2,2) = 0.05; % weighting matrix (Input)

%compute objective

for i = 1:N
    st = X(:,k); con = U(:,k);
    obj = obj + (st-P(4:6))'*Q*(st-P(4:6)) + con'*R*con;        %calculate obj
end

%compute constraints
for k = 1:N+1
   g = [g; X(1,k)];         %state x
   g = [g; X(2,k)];         %state Y
end

%NON linear programmin structre
OPT_variables = reshape(U, 2*N, 1);
nlp_prob = struct('f',obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level = 0;
opts.print_time = 0;
opts.ipopt.acceptable_tol = 1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

%arguments

args = struct;
%inequality constrains
args.lbg = -2;          %lower bound of states
args.ubg = 2;           %upper bound of states

%input constraints
args.lbx(1:2:2*N-1,1)= v_min;    
args.lbx(2:2:2*N,1)= omega_min;
args.ubx(1:2:2*N-1,1)= v_max;
args.ubx(2:2:2*N,1)= omega_max;


%--------------------------
%SIMULATION LOOP SHOULD STARTS FROM HERE
%------------------------------------------

t0 = 0;
x0 = [0;0;0.0];         %initial conditions
xs = [1.5;1.5;0.0];     %reference posture
xx(:,1) = x0;           %xx contains the history of states
t(1) = t0;
u0 = zeros(N,2);        %   two control inputs
sim_tim = 20;           %maximum simulation time

%start MPC
mpciter = 0;
xx1= [];
u_cl = [];

while (norm((x0-xs),2) > 1e-2 && mpciter < sim_tim/T)
   
    args.p = [x0;xs];
    args.x0 = reshape(u0',2*N, 1);
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx',args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p);
    u = reshape(full(sol.x)', 2, N)';
    ff_value = ff(u', args.p);
    xx1(:,1:3,mpciter+1) = full(ff_value)';
    u_cl = [u_cl ; u(1,:)];
    
    t(mpciter+1) =t0;
    [t0,x0,u0] = shift(T, t0, x0, u,f);
    xx(:,mpciter+2) = x0;
    mpciter = mpciter +1;
end

function [t0,x0,u0] = shift(T, t0, x0, u,f)
st = x0;
con = u(1,:)';
f_value =   f(st, con);
st = st + (T+f_value);
x0 = full(st);

t0 = t0 + T;
u0 = [u(2:size(u,1),:) ; u(size(u,1),:)];
end
