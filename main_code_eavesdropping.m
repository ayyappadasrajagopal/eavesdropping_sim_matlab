%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HIL based State Estimation and Control with Quanser Qube servo 3 with
% AWGN channel and BPSK modulated - single sensor and fusion based
% -------------------------------------------------------------------------
% Coded by:
% ---------
% Ayyappadas R
% Research Scholar
% Dept. of EE4 
% IIT Palakkad
%
% Descirption:
% ------------
% This is the main code that simulates the HIL system for brief 
% -> Quanser cube servo as plant
% -> Estimation and control, closed loop operation
% -> Theoretical and simulated BPSK - BER and PER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clc;
close all;
clear all;

%% Initializations
%%-------------------------------------------------------------------------
N = 100;           % number of iterations / number of samples / number of packets
N_A = 2;           % Count for averaging
n = 2;             % order of system / no of states
prob = 0:0.1:1;    % probability of arrival/ avg arrival rate - for measurement channel, 
                   % - we could add seprate control channel 
                   % note: min probability of 0.02 is given to avoid bad conditioning error as
                   % error tends to infinity
                   % f = waitbar(0, 'Starting'); % wait bar initialization

%% Quanser cube Speed/Position control Modeling
%%-------------------------------------------------------------------------
qube2_param; % Load parameters
% State-space matrices (with states x = [theta; theta_dot])
A = [0 1;
    0 -kt*km/Jeq/Rm];
B = [0; kt/Jeq/Rm];
C = eye(2);
D = zeros(2,1);
motor_ss=ss(A,B,C,D);

%% LQR Design
%%-------------------------------------------------------------------------
% % Bryson's Rule
u_max = 5; % max current based on 10V max voltage
x1_max = 1; % max angle (rad)
x2_max = u_max/km; % max angular rate (rad/s)
r11 = 1/u_max^2;
q11 = 1/x1_max^2;
q22 = 1/x2_max^2;

% LQR
Q=diag([q11 q22]); 
R=r11; 

% Convert to discrete, where dt is your discrete time-step (in seconds)
dt = 0.1; 
d_sys = c2d(motor_ss,dt,'tustin');

A_d  = d_sys.A;
B_d  = d_sys.B;
C_d  = d_sys.C;
D_d  = d_sys.D;

Q_c = 0.001*[1 0;
              0 1]; % Controller's weight matrix for the state
R_c = 0.001;       % Controller's weight matrix for the input

% Q_c = Q;
% R_c = R;

[K_d,P_d,E_d] = dlqr(A_d,B_d,Q_c,R_c);

% Construct Discrete Time State-Space Matrices
A = A_d;
B = B_d;
C = [1 0;
    0   1];
C_a = [1 0];
C_rem = C; % if homogenous
D = zeros(2,1);

%% Controllability and observability
%%-------------------------------------------------------------------------
Cab = ctrb(A,B);
Ctrb_rank = rank(Cab);
Oac = obsv(A,C);
Obs_rank = rank(Oac);

poles_discrete = eig(A);       % Open Loop Poles (Discrete Time)
[K_c,S,e] = dlqr(A,B,Q_c,R_c); % Find LQR gain Matrix K and new poles e 

Q = 0.0001*eye(n); 
R = 0.0001;
R_rem = 0.01;          % if homogenous
R_f = 1/((1/R)+(1/R_rem));

sqrtQ = sqrtm(Q);   % std.dev of process noise
sqrtR = sqrt(R);
sqrtR_rem = sqrt(R_rem);

%% Closed Loop Discrete Time State-Space Model
%%-------------------------------------------------------------------------
m1 = size(C,1);                   % number of outputs for sensor 1
m2 = size(C,2)
m1_rem = size(C_rem,1);
Sigma = eye(n);                   % initial error covaraince value
x(:,1)= sqrtm(Sigma)*randn(n,1);  % initializing vector for system dyanmics
y=zeros(m1,N);                    % output vector initialized to all zero values

%% SubSystem Dynamics
%%-------------------------------------------------------------------------

x1(:,1)= zeros(n,1);   % initializing vector for system dyanmics
y1a=zeros(m1,N+1);     % output vector initialized to all zero values

open_mse1 = zeros(n,n);
open_mse1_f = zeros(n,n);
open_mse1_c = zeros(n,n);
open_empiricalmse1 = zeros(n,n);
open_empiricalmse1_f = zeros(n,n);
open_empiricalmse1_c = zeros(n,n);

%% Working loops
%%-------------------------------------------------------------------------

for kkb = 1:length(prob)
p = prob(kkb);              % probability/rate of arrival is inherent by
                            % the SNR                         
open_mse1_iter = zeros(n);  % initailize error covaraince matrices for averaging loop
open_empiricalmse1_iter = zeros(n); 
open_empiricalmse1_c_iter = zeros(n);
open_mse1_c_iter = zeros(n);
x_des = [45 0];
x_des = x_des*(pi/180);

for kk=1:N_A % System dynamics averaging 

%%%%%%%%%%%%% For the first candidate only one time %%%%%%%%%%%%%%%%%%%%%%%
y(:,1) = C*x(:,1) + sqrtR*randn(m1,1);
y_a(:,1) = C_a*x(:,1) + sqrtR*randn(m2,1);
y_rem(:,1) = C_rem*x(:,1) + sqrtR_rem*randn(m1_rem,1); % Homogenous sensor 
y_f(:,1) = y(:,1)+ (R/(R+R_rem))*(y_rem(:,1)-y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of error covarinace matrices
P1 = eye(n); % Error covarinace matrices initialized to identity maytrices.
P1_f = eye(n);
P1_c = eye(n);

zeta = 1;  % varibale for packet arrival initilaized with one

open_empiricalmse = zeros(1,n); %average error*error
open_mse = zeros(n);            %average P_k
open_mse_c = zeros(n);
open_empiricalmse_c = zeros(n);

hatx1 = zeros(n,1);
hatx1_c = zeros(n,1);
hatx1_f = zeros(n,1);
hatx1_rem = zeros(n,1);

open_error1_vec = zeros(n,N);

for k = 1:N

K_opt= dlqr(A,B,Q_c,R_c);
u = K_opt*(x_des'-hatx1);
ub(k+1) = u;

% System dynamics
wk = real(sqrtQ*randn(n,1));     % process noise
x(:,k+1) = A*x(:,k) + B*u+ wk;   % system/state dynamics

% Subsystem 1
vk = sqrtR*randn(m1,1);          % measurement noise
vk_rem = sqrtR_rem*randn(m1,1);  % measurement noise of remote sensor

y(:,k+1) = C*x(:,k+1) + vk;
y_rem(:,k+1) = C_rem*x(:,k+1) + vk_rem;
y_f(:,k+1) = y(:,k+1)+ (R/(R+R_rem))*(y_rem(:,k+1)-y(:,k+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Local State Estimation 
%%------------------------

    zeta = double(rand(1,1) <= p);
  
    %prediction step
	hatx1 = A*hatx1 +B*u; 
    P1 = A * P1 * A' + Q;
    % Correction step
    L1 = P1 * C'/( C*P1*C' + R);
    if zeta~=0 % on if no local hardwired sensors are available
    hatx1 = (eye(n) - L1*C)*hatx1 + L1*y_rem(:,k);
    P1 = P1 - L1*C*P1;
    seq(1,k) = zeta;
    end
    hatx1_ar(:,k+1)= hatx1;

    open_mse1 = open_mse1 + P1;
    open_error1 = x(:,k) - hatx1;
    open_empiricalmse1 = open_empiricalmse1 + open_error1 * open_error1';
    open_error1_vec(:,k) = open_error1;   

%%  2. Fusion based State Estimation 
%%------------------------------------ 

    %prediction step
	hatx1_f = A*hatx1_f +B*u; 
    P1_f = A * P1_f * A' + Q;
    hatx1_ar_f(:,k+1)= hatx1_f;
    
    % Correction step
    L1_f = P1_f * C'/( C*P1_f*C' + R_f);

    y_f(:,k+1) = y(:,k+1)+ (R/(R+R_rem))*(y_rem(:,k+1)-y(:,k+1));
    if zeta~=0 
    hatx1_f = (eye(n) - L1_f*C)*hatx1_f + L1_f*y_f(:,k);
    P1_f = P1_f - L1_f*C*P1_f;
    seq_f(1,k) = zeta;
    end

    open_mse1_f = open_mse1_f + P1_f;
    open_error1_f = x(:,k) - hatx1_f;
    open_empiricalmse1_f = open_empiricalmse1_f + open_error1_f * open_error1_f';

end

open_empiricalmse1 = open_empiricalmse1/N;   % Empirical MSE - single remote sensor
open_mse1 = open_mse1/N;                     % Theoretical MSE - single remote sensor
tr_empmse1 = trace(open_empiricalmse1);      % trace of empirical MSE
tr_mse1 =  trace(open_mse1);                 % trace of theoretical MSE
tr_empmse1_arr(kk) = tr_empmse1;             % array tr MSE
tr_mse1_arr(kk) = tr_mse1;                   % array tr MSE

open_empiricalmse1_f = open_empiricalmse1_f/N; % Empirical MSE - fusion
open_mse1_f = open_mse1_f/N;                   % Theoretical MSE - fusion
tr_empmse1_f = trace(open_empiricalmse1_f);    % trace of empirical MSE
tr_mse1_f =  trace(open_mse1_f);               % trace of theoretical MSE
tr_empmse1_arr_f(kk) = tr_empmse1_f;           % array tr MSE
tr_mse1_arr_f(kk) = tr_mse1_f;                 % array tr MSE

end

% average the trace of empirical MSE over N_A number of iterations
tr_empmse1_arr_a = sum(tr_empmse1_arr)/N_A;    % Single remote sensor based
tr_mse1_arr_a = sum(tr_mse1_arr)/N_A;
% save the trace values for probability loop  
tr_empmse1_arr_ab(kkb) = tr_empmse1_arr_a;     % corresponding to system 1 - local estimator 
tr_mse1_arr_ab(kkb) = tr_mse1_arr_a; 

tr_empmse1_arr_af = sum(tr_empmse1_arr_f)/N_A; % Fusion based 
tr_mse1_arr_af = sum(tr_mse1_arr_f)/N_A;
tr_empmse1_arr_abf(kkb) = tr_empmse1_arr_af;   % corresponding to system 1 - fusion
tr_mse1_arr_abf(kkb) = tr_mse1_arr_af; 

end

t1 = 0:1:N;

figure;
subplot(2,1,1)
plot(t1,x(1,:)*(180/pi),'b',LineWidth=1)
hold on;
plot(t1,hatx1_ar(1,:)*(180/pi),'r',LineWidth=1)
xlabel('Time instants $k$','Interpreter','latex')
ylabel('$\theta$ and $\hat{\theta}$ in degrees','Interpreter','latex')
legend('Actual State $\theta$','Estimated  State $\hat{\theta}$','Interpreter','latex');
title('Position')
subplot(2,1,2)
plot(t1,x(2,:)*(180/pi),'b',LineWidth=1)
hold on;
plot(t1,hatx1_ar(2,:)*(180/pi),'r',LineWidth=1)
xlabel('Time instants $k$','Interpreter','latex')
ylabel('$\omega$ and $\hat{\omega}$ in degrees/sec','Interpreter','latex')
legend('Actual State - $\omega$','Estimated  State - $\hat{\omega}$','Interpreter','latex');
title('Angle')
sgtitle('Actual and Estimated states')

figure;
subplot(2,1,1)
plot(t1,x(1,:)*(180/pi)-hatx1_ar(1,:)*(180/pi),'b',LineWidth=1)
xlabel('Time instants $k$','Interpreter','latex')
ylabel('$\theta-\hat{\theta}$ in degrees','Interpreter','latex')
title('Error in Position')
subplot(2,1,2)
plot(t1,x(2,:)*(180/pi)-hatx1_ar(2,:)*(180/pi),'b',LineWidth=1)
xlabel('Time instants $k$','Interpreter','latex')
ylabel('Error in estimates for $\omega-\hat{\omega}$ in degrees/sec','Interpreter','latex')
title('Error in Angle')
sgtitle('Error in Estimated states')

figure;
plot(t1,ub,'b',LineWidth=1)
xlabel('Time instants $k$','Interpreter','latex')
ylabel('Control - Voltage in Volts')
title('Control')

figure;
sgtitle("Tr(MSE) vs. Avg arrival rate for System 1, System 4 and System 4 with consensus");
% subplot(2,1,1)
plot(prob, smooth(smooth(tr_empmse1_arr_ab*(180/pi))))
hold on;
%plot(prob, smooth(tr_mse1_arr_ab))
plot(prob, smooth(smooth(tr_empmse1_arr_abf*(180/pi))))
%plot(prob, smooth(tr_mse1_arr_abf))
legend('Empirical - single measurement','Theoretical - single measurement', 'Empirical - Fusion based', 'Theoretical - Fusion based')
xlabel('Avg arrival rate');
ylabel('Tr(MSE)');








% Fusion based Estimation

% Consensus based Estimation