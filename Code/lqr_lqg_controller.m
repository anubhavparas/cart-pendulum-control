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

% double derivative of th1
th1dd = (xdd*cos(th1) - g*sin(th1))/l1;

% double derivative of th2
th2dd = (xdd*cos(th2) - g*sin(th2))/l2;

%% Non-linear function
f1 = xd;
f2 = xdd;
f3 = th1d;
f4 = th1dd;
f5 = th2d;
f6 = th2dd;

F = [f1 f2 f3 f4 f5 f6];

%% Linearization
%% 
% Jacobian
J = jacobian(F, q);

% Equilibrium points:
q_e = [0 0 0 0 0 0];
J = subs(J, q, q_e);
disp(J);

%% State-space model
%%
A = J;

B = jacobian(F, u);
B = subs(B, q, q_e);

C = eye(size(A));

D = 0;



%% Checking controllability
%%
CTRB = [B A*B (A^2)*B (A^3)*B (A^4)*B (A^5)*B];
det_ctrb = det(CTRB);
pretty(det_ctrb);


%% Substituting the values of constants
%%
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 10;

A = double(subs(A));
B = double(subs(B));
disp(A);
disp(B);

%% Check contrallability with values
%%
%disp(rank(crtb(A,B)));


%% Check stability
%%
poles_open = eigs(A);
disp('Open loop poles')
disp(poles_open);


%% LQR Design
%%
Q = zeros(6);
Q(1,1) = 1000;
Q(2,2) = 10;
Q(3,3) = 1000000;
Q(4,4) = 10;
Q(5,5) = 10000000;
Q(6,6) = 10;

R = 0.0001;

Kr = lqr(A,B,Q,R);
disp(Kr);


%% Closed Loop State-space
Ac = A-B*Kr;

%% Check closed loop stability

poles_cl = eigs(Ac);
disp(poles_cl);

ss_open = ss(A, B, C, D);
%% Simulation:
x_o = 0.5; th1_o = deg2rad(10); th2_o = deg2rad(10);
init_state = [x_o, 0, th1_o, 0, th2_o, 0];
ss_cl = ss(Ac, B, C, D);

[y_op,t_op,x_op] = initial(ss_open, init_state);
[y,t,x] = initial(ss_cl, init_state);

% hold on
% plot(t(:,1),y(:,1),t(:,1),y(:,3), t(:,1),y(:,5));
% legend x theta1 theta2
% title 'Unit step response of LQR Controller Qx=1000, Qth1=1000000, Qth2=10000000 , R=0.0001'
% hold off

% hold on
% plot(t(:,1),y(:,1));
% legend x
% title 'Unit step response of LQR Controller Q=1000, R=0.0001'
% hold off

% hold on
% plot(t(:,1),y(:,3));
% legend theta1
% title 'Unit step response of LQR Controller Q=1000000, R=0.0001'
% hold off
% 
% 
% hold on
% plot(t(:,1),y(:,5));
% legend theta2
% title 'Unit step response of LQR Controller Q=10000000, R=0.0001'
% hold off
%% Observability Check
%% Case 1: x(t)
C1 = [1 0 0 0 0 0];
obsv_mat1 = obsv(A,C1);
disp(rank(obsv_mat1));  % rank=6 - observable


%% Case 2: th1(t), th2(t) 
C2 = [0 0 1 0 0 0; 0 0 0 0 1 0];
obsv_mat2 = obsv(A,C2);
disp(rank(obsv_mat2)); % rank=4 - not-observable


%% Case 3: x(t), th2(t) 
C3 = [1 0 0 0 0 0; 0 0 0 0 1 0];
obsv_mat3 = obsv(A,C3);
disp(rank(obsv_mat3)); % rank=6 - observable


%% Case 4: x(t), th1(1), th2(t) 
C4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
obsv_mat4 = obsv(A,C4);
disp(rank(obsv_mat4)); % rank=6 - observable


%% Luenberger Observer
%%
% We can select the poles of the observer control gain to be less than
% those of the full-state feedback controller so that the error dynamics
% converge to stability faster than the systems state.
% leftmost pole for the Kr is -3.9
%poles_obsv = [-4.0 -4.05, -4.1, -4.15, -4.2, -4.25];
poles_obsv = [-1.6 -1.75, -1.8, -2.15, -2.2, -2.5];

%% Case1
L1 = place(A', C1', poles_obsv)';
Aob1 = A-L1*C1;
Bob1 = [B L1];
Cob1 = eye(size(A));
Dob1 = zeros(6,2);


ss_ob1 = ss(Aob1, Bob1, Cob1, Dob1);
ss_op1 = ss(A, B, C1, D);


dt = 0.01;
t_sim = dt:dt:50;
inp_u = ones(size(t_sim));
[y_f_op, t_] = lsim(ss_open, inp_u, t_sim);
[y_op, t_] = lsim(ss_op1, inp_u, t_sim);
[x_est, t_] = lsim(ss_ob1, [inp_u; y_op'], t_sim);


% hold on
% subplot(3,1,1);
% plot(t_(:,1), y_f_op(:,1),'k', t_(:,1), x_est(:,1), 'r--', 'Linewidth',2);
% legend x x-estimate
% title 'Observer tracking - x'
% hold off
% 
% hold on
% subplot(3,1,2)
% plot(t_(:,1), y_f_op(:,3),'k', t_(:,1), x_est(:,3), 'r--', 'Linewidth',2);
% legend th1 th1-estimate
% title 'Observer tracking - th1'
% hold off
% 
% hold on
% subplot(3,1,3)
% plot(t_(:,1), y_f_op(:,5),'k', t_(:,1), x_est(:,5), 'r--', 'Linewidth',2);
% legend th2 th2-estimate
% title 'Observer tracking - th2'
% hold off
 

%% Case3
L3 = place(A', C3', poles_obsv)';
Aob3 = A-L3*C3;
Bob3 = [B L3];
Cob3 = eye(size(A));
Dob3 = zeros(6,3);


ss_ob3 = ss(Aob3, Bob3, Cob3, Dob3);
ss_op3 = ss(A, B, C3, D);


dt = 0.01;
t_sim = dt:dt:50;
inp_u = ones(size(t_sim));
[y_f_op, t_] = lsim(ss_open, inp_u, t_sim);
[y_op, t_] = lsim(ss_op3, inp_u, t_sim);
[x_est, t_] = lsim(ss_ob3, [inp_u; y_op'], t_sim);


% hold on
% subplot(3,1,1);
% plot(t_(:,1), y_f_op(:,1),'k', t_(:,1), x_est(:,1), 'r--', 'Linewidth',2);
% legend x x-estimate
% title 'Observer tracking - x'
% hold off
% 
% hold on
% subplot(3,1,2)
% plot(t_(:,1), y_f_op(:,3),'k', t_(:,1), x_est(:,3), 'r--', 'Linewidth',2);
% legend th1 th1-estimate
% title 'Observer tracking - th1'
% hold off
% 
% hold on
% subplot(3,1,3)
% plot(t_(:,1), y_f_op(:,5),'k', t_(:,1), x_est(:,5), 'r--', 'Linewidth',2);
% legend th2 th2-estimate
% title 'Observer tracking - th2'
% hold off


%% Case4
L4 = place(A', C4', poles_obsv)';
Aob4 = A-L4*C4;
Bob4 = [B L4];
disp(size(Bob4));
Cob4 = eye(size(A));
Dob4 = zeros(6,4);


ss_ob4 = ss(Aob4, Bob4, Cob4, Dob4);
ss_op4 = ss(A, B, C4, D);


dt = 0.01;
t_sim = dt:dt:50;
inp_u = ones(size(t_sim));
[y_f_op, t_] = lsim(ss_open, inp_u, t_sim);
[y_op, t_] = lsim(ss_op4, inp_u, t_sim);
[x_est, t_] = lsim(ss_ob4, [inp_u; y_op'], t_sim);


% hold on
% subplot(3,1,1);
% plot(t_(:,1), y_f_op(:,1),'k', t_(:,1), x_est(:,1), 'r--', 'Linewidth',2);
% legend x x-estimate
% title 'Observer tracking - x'
% hold off
% 
% hold on
% subplot(3,1,2)
% plot(t_(:,1), y_f_op(:,3),'k', t_(:,1), x_est(:,3), 'r--', 'Linewidth',2);
% legend th1 th1-estimate
% title 'Observer tracking - th1'
% hold off
% 
% hold on
% subplot(3,1,3)
% plot(t_(:,1), y_f_op(:,5),'k', t_(:,1), x_est(:,5), 'r--', 'Linewidth',2);
% legend th2 th2-estimate
% title 'Observer tracking - th2'
% hold off



%% LQG Output state feedback controller design
%%
%% Optimal state estimator design:
C_out = [1 0 0 0 0 0]; %% smallest output vector - observable

D_out = zeros(size(C,1), size(B,2));

Vd = 0.001*eye(size(A));  % Process disturbance
Vn = 0.0001; % measurement noise

Kf = lqe(A, Vd, C_out, Vd, Vn);  % optimal Kalman filter gain
A_kf = A - Kf*C_out;
B_kf = [B Kf];
C_kf = eye(size(A));
D_kf = 0*[B Kf];


sysKF = ss(A_kf, B_kf, C_kf, D_kf);

%% LQG System = LQR + Kalman Filter
%
A_lqg = [(A-B*Kr), B*Kr;
         zeros(size(A)), (A*Kf*C_out)   
        ];

B_lqg = [B; zeros(size(B))];

C_lqg = [C_out zeros(size(C_out))];

sys_lqg = ss(A_lqg, B_lqg, C_lqg, 0);

[y_lqg, t, x_lqg] = lsim(sys_lqg, inp_u, t_sim, [init_state init_state]);

num_states = 6;
out_state = x_lqg(:, 1:num_states);
err = x_lqg(:, num_states+1:end);

state_est = out_state - err;

x_res = out_state(:,1);
th1_res = out_state(:,3);
th2_res = out_state(:,5);

x_est = state_est(:,1);
th1_est = state_est(:,3);
th2_est = state_est(:,5);

hold on
subplot(3,1,1);
plot(t,x_res,'-r', t,x_est,':r')
legend x-response x-est
title 'LQG Controller reponse';

subplot(3,1,2);
plot(t,th1_res,'-b', t,th1_est,':b')
legend th1-reponse th1-est
title 'LQG Controller reponse';

subplot(3,1,3);
plot(t,th2_res,'-g', t,th2_est,':g');
legend th2-response th2-est;
title 'LQG Controller reponse';
