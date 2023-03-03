% This program computes the pseudospectra of a matrix using two different
% curve tracing methods.
%
% PF(lamda0,epsilon,A,tol2,tk,flag)
% Cobra2(lamda0,epsilon,A,tol2,h1,H,flag)

%% ---------- 1st Example ----------
clear all;
clc;
A = gallery('redheff',64);

% epsilon = 1
PF(6,1,A,1e-2,9e-2,3);
PF(2,1,A,3e-2,9e-2,3);
PF(-1,1,A,5e-2,9e-2,3);

% epsilon = 1.23
PF(1,1.23,A,3e-1,5e-1,3);
PF(1,1.23,A,1e-2,5e-2,3);
PF(1,1.23,A,9e-3,1e-2,3);
PF(8,1.23,A,1e-2,9e-2,3);

% epsilon = 2.9
PF(8,2.9,A,9e-2,1e-1,3);

%% ---------- 2nd Example ----------
clear all;
clc;
A = gallery('kahan',50);

% epsilon = 1e-2
PF(-1,1e-2,A,4e-2,1e-1,3);
PF(-1,1e-2,A,4e-2,8e-2,3);
PF(-1,1e-2,A,9e-3,2e-2,3);
Cobra(-1,1e-2,A,2e-2,5e-3,4.4,3); %210,1673
% Cobra(-1,1e-2,A,9e-3,5e-3,2.4,3); %385,3073
% Cobra(-1,1e-2,A,6e-3,3e-3,4.4,3); %350,2793

% epsilon = 1e-4
PF(-1,1e-4,A,4e-2,1e-1,3);
PF(-1,1e-4,A,4e-2,8e-2,3);
PF(-1,1e-4,A,9e-3,21e-3,3);
Cobra(-1,1e-4,A,9e-3,5e-3,3.1,3);

%% ---------- 3rd Example ----------
clear all;
clc;
A = gallery('grcar',64);

% epsilon = 1e-2
PF(1,1e-2,A,7e-2,1e-1,3);
PF(1,1e-2,A,1e-2,5e-2,3);
PF(1,1e-2,A,1e-2,25e-3,3);
Cobra(1,1e-2,A,1e-2,15e-3,0.9,3);

% epsilon = 1e-3
PF(1,1e-3,A,7e-2,1e-1,3);
PF(1,1e-3,A,1e-2,5e-2,3);
PF(1,1e-3,A,1e-2,25e-3,3);
Cobra(1,1e-3,A,5e-3,15e-3,0.4,3);
