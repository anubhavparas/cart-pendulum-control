%% Start
clear all;
clc;
close all;

%% Symbols
syms F M m1 m2 l1 l2 g;
syms x xd xdd th1 th1d th1dd th2 th2d thdd;

%% Defining state variables
q = [x xd th1 th1d th2 th2d];
u = [F];

%% Defining the eqns - Non-Linear system

% double derivative of x
xdd_num = F - m1*g*cos(th1)*sin(th1) - - m2*g*cos(th2)*sin(th2) - m1*l1*(th1d^2)*sin(th1) - m2*l2*(th2d^2)*sin(th2);
xdd_den = M + m1*(sin(th1))^2 + m2*(sin(th2))^2;
xdd = xdd_num/xdd_den;
syms xdd;
% double derivative of th1
th1dd = (xdd*cos(th1) - g*sin(th1))/l1;

% double derivative of th2
th2dd = (xdd*cos(th2) - g*sin(th2))/l2;

% M = 1000;
% m1 = 100;
% m2 = 100;
% l1 = 20;
% l2 = 10;
% g = 10;


vals = [1000 100 100 20 10 10];
xdd = subs(xdd, [M m1 m2 l1 l2 g], vals);
th1dd = subs(th1dd, [M m1 m2 l1 l2 g], vals);
th2dd = subs(th2dd, [M m1 m2 l1 l2 g], vals);
disp(th2dd);

u = [th1 th1d th2 th2d F];
f1 = -(2000*sin(u(1))*u(2)^2 + 1000*sin(u(3))*u(4)^2 - u(5) + 1000*cos(u(1))*sin(u(1)) - 1000*cos(u(3))*sin(u(3)))/(100*sin(u(1))^2 + 100*sin(u(3))^2 + 1000);
 
u = [xdd th1];
f2 = (u(1)*cos(u(2)))/20 - sin(u(2))/2;
u = [xdd th2];
f3 = (u(1)*cos(u(2)))/10 - sin(u(2));




