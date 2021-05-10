%% Seminar Control

%Problem: stabilize the single pendulum using LQR

clc;
clear;
clear all;

%% state space model repersentation

% Initial conditions
X0 = [];

% System dynamics A,B,C,D
A = [];
B = [];
C = [];
D = [];

% Control parameters
Q = [];
R = 1;
K = lqr(A, B, Q, R);

% Closed Loop system
sys = ss((A - B*K), B, C, D);

% Run response to intialize the condition
t = 0:0.005:30;
[y,t,x] = initial(sys, x0, t);