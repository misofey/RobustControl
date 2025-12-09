load('Assignment_Data_SC42145_2025.mat');
% States: omega, z1_dot, z1, z2_dot, z2
% Inputs: beta, tau(generator torque), wind_speed
ops =  bodeoptions;
s =tf('s');

mimo = FWT(1:2, 1:2);
G = tf(mimo);
Gd = FWT(1:2, 3);
%%% RGA %%%
% bode(mimo)

omega = 0.3*2*pi;
G_crf = freqresp(mimo, omega);
RGA_crf = G_crf.*pinv(G_crf).';

omega = 0;
G_0 = evalfr(mimo, omega);
RGA_0 = G_0.*pinv(G_0).';

%%% Performance weight %%%
% Low frequency attunation
a = 10e-4; 
% H infinity norm for sensitivity
M = 3;
% Cut off frequency
omega_b = 0.3 * 2*pi;
    
% W_p1 = inv(makeweight(a, omega_b, M));

W_p1 = (s/M + omega_b) / (s + omega_b*a);

W_p2 = 0.1;
W_p = [W_p1 0; 0 W_p2];

W_u1 = 0.01;
W_u2 = (5e-3*s^2 + 7e-4*s + 5e-5)/(s^2 + 14e-4*s + 1e-6); 
W_u = [W_u1 0;
       0    W_u2];

%% Mixed-sensitivity design %%
% [K,CL,gamma,INFO] = mixsyn(mimo,W_p,W_u,[]);
% L = tf(minreal(mimo*K));
% I = eye(size(L));
% S = feedback(I,L); 
% T = I-S;
% ops.Title.String = "Sensitivity";
% bode(S, T, L)
% sigma(CL, gamma)
% legend();
%step(feedback(L, I))

% % HINF DESIGN %%
% P = minreal(augw(G, W_p, W_u, []));
% [K,CL,gamma,INFO] = hinfsyn(P, 2, 2);
[K,CL,gamma,INFO] = mixsyn(G, W_p, W_u, []);

% sigma(CL, ss(gamma))

L=G*K;
S = feedback(eye(2),L);
T = G*K*S;
KS = K*S;


hold on
grid
bode(S(1, 1))
bode(1/W_p1)
legend("S", "W_p1");
hold off
% step(feedback(L1, I), feedback(L2, I))

figure;
hold on
bode(KS(2, 1));
bode(1/W_u2);
grid;
legend("KS", "W_u2");
hold off

figure;
hold on
sigma(T)
hold off
