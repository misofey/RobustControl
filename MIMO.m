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
P = minreal(augw(G, W_p, W_u, []));
P = [zeros(2) W_u;
    W_p W_p*G;
    -eye(2) -G];
[K,CL,gamma,INFO] = hinfsyn(P, 2, 2);
% [K,CL,gamma,INFO] = mixsyn(G, W_p, W_u, []);

% sigma(CL, ss(gamma))
%%
K.OutputName = ["\beta", "\tau_r"];
L=minreal(G*K);
S = feedback(eye(2),L);
S.InputName = ["d_{\omega_r}", "d_z"];
S.OutputName = ["\omega_r", "z"];
T = G*K*S;
KS = K*S;

K_poles = pole(K);
G_poles = pole(G);

if(max(real(K_poles))<0 && max(real(G_poles))<0)
    disp("The system is internally stable")
else
    disp("The system is not internally stable")
    max_K_pole = max(real(K_poles))
    max_G_pole = max(real(G_poles))
end

%%
opts = bodeoptions;
opts.PhaseVisible = "off";

K_MS = K;
save("K_MS", "K_MS");
plot_sensitivity = false;
plot_controller_sensitivity = false;
plot_time_simulations = true;

if plot_sensitivity
    figure;
    hold on
    grid
    bodeplot(S(1, 1), opts)
    bodeplot(1/W_p1, opts)
    legend("S", "1/W_{p11}");
    title("")
    hold off
    
    figure;
    hold on
    grid
    bodeplot(S(2, 1), opts)
    bodeplot(1/W_p(2, 2), opts)
    legend("S", "1/W_{p22}");
    title("")
    hold off
    
    figure;
    hold on
    grid
    bodeplot(S(1, 2), opts)
    bodeplot(1/W_p(1, 1), opts)
    legend("S", "1/W_{p11}");
    title("")
    hold off
    
    figure;
    hold on
    grid
    bodeplot(S(2, 2), opts)
    bodeplot(1/W_p(2, 2), opts)
    legend("S", "1/W_{p22}");
    title("")
    hold off
end

if plot_controller_sensitivity
    figure;
    hold on
    grid
    bodeplot(KS(1, 1), opts)
    bodeplot(1/W_u(1, 1), opts)
    legend("KS", "1/W_{u11}");
    title("")
    hold off
    
    figure;
    hold on
    grid
    bodeplot(KS(2, 1), opts)
    bodeplot(1/W_u(2, 2), opts)
    legend("KS", "1/W_{u22}");
    title("")
    hold off
    
    figure;
    hold on
    grid
    bodeplot(KS(1, 2), opts)
    bodeplot(1/W_u(1, 1), opts)
    legend("KS", "1/W_{u11}");
    title("")
    hold off
    
    figure;
    hold on
    grid
    bodeplot(KS(2, 2), opts)
    bodeplot(1/W_u(2, 2), opts)
    legend("KS", "1/W_{u22}");
    title("")
    hold off
end

if plot_time_simulations
    stepWind(FWT, -K)
end