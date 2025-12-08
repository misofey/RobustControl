load('Assignment_Data_SC42145_2025.mat');

complete_plant = ss(A, B, C, D);
G = tf(complete_plant(1:2, 1:2));

% States: omega, z1_dot, z1, z2_dot, z2
% Inputs: beta, tau(generator torque), wind_speed
ops =  bodeoptions;
s =tf('s');
mimo = ss(A, B, C, D);

mimo = tf(mimo);
mimo = minreal(-mimo);

mimo = FWT(1:2, 1:2);
%%% RGA %%%
bode(mimo)

omega = 0.3*2*pi;
G_crf = freqresp(mimo, omega);
RGA_crf = G_crf.*pinv(G_crf).';

omega = 0;
G_0 = evalfr(mimo, omega);
RGA_0 = G_0.*pinv(G_0).';

%%% Controller design %%%
% Low frequency attunation
a = 10e-4; 
% H infinity norm for sensitivity
M = 3;
% Cut off frequency
omega_b = 0.3*2*pi;

% Controller weights
W_p1 = inv(tf(makeweight(10^-4, [omega_b, 0.708], 3)));

%W_p1 = minreal(10^(-4)*(s*omega_b*10^(3.5)+1)/(s*omega_b*10^(-0.627)+1)); %(s*M^(-1) + omega_b)/(s + omega_b*a);
%W_p1 = minreal(10^(-4)*(s/(-1.775^(-4))-1)/(s/-5.325-1))
W_p2 = 0.05;
W_p = [W_p1 0; 0 W_p2];

W_u1 = 0.005;
W_u2 = (5e-3*s^2 + 7e-4*s + 5e-5)/(s^2 + 14e-4*s + 1e-6); 
W_u = [W_u1 0    0;
       0    W_u2 0 ;
       0    0    0];

%% Mixed-sensitivity design %%
[K,CL,GAM,INFO] = mixsyn(mimo,W_p,W_u,[]);
L1 = tf(minreal(mimo*K));
I = eye(size(L));
S = feedback(I,L); 
T= I-S;
ops.Title.String = "Sensitivity";
bode(S, T, L)
sigma(CL, ss(GAM))
%step(feedback(L, I))

%% HINF DESIGN %%
P = augw(mimo, W_p, W_u, []);
[K,CL,GAM,INFO] = hinfsyn(P);

L2 = tf(minreal(mimo*K));
I = eye(size(L));
%sigma(CL, ss(GAM))
%step(minreal(CL))
%step(feedback(L1, I), feedback(L2, I))



grid