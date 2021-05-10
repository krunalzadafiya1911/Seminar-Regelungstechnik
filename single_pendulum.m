function dx =  single_pendulum(t, x)
mass = 1;
l = 0.5;
g = 9.8; %gravitational acceleration
    
%define the state space model
dx = [x(2,1); -(g/l)*sin(x(1,1))];

end