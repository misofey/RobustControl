load('Assignment_Data_SC42145_2025.mat');
load("K_MS");
mimo = FWT(1:2, 1:2);
G = tf(mimo);
Gd = FWT(1:2, 3);
% G = G;

G.InputName = ["beta", "tau"];
G.OutputName = ["omega", "z"];
Gd.InputName = ["d"];
Gd.OutputName = ["d_omega", "d_z"];

%% Performance weights

s = tf("s");
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


% W_p.u = ["omega_dist", "z_dist"];
% W_p.y = "z1";
W_p.u = ["e_omega", "e_z"];
W_p.y = ["z1"];
W_u.u = ["beta", "tau"];
W_u.y = "z2";

%% interconnection
% beta_controller = tunableTF('K_beta', 2, 2) + tunablePID("K_beta2", "PD");
beta_omega_LL = tunableTF('beta_omega_LL', 2, 2);
beta_omega_PID = tunablePID("beta_omega_PID", "PD");


% tau_controller = tunableTF('K_tau', 4, 4);
tau_omega_PID = tunablePID("tau_omega_PID", "PD");
tau_omega_LL = tunableTF("tau_omega_LL", 2, 2);

beta_z_pid = tunablePID("beta_z_pid", "PD");


% tau_z_PID = tunablePID("gain_tauz", "P");
% tau_z = tunableTF('K_beta', 4, 4);
% tau_z.u = 'e_z';
% tau_z.y = 'tau';

K = [beta_omega_LL + beta_omega_PID, 0; tau_omega_PID + tau_omega_LL, 0];
K.InputName = ["e_omega", "e_z"];
K.OutputName = ["beta", "tau"];


%%
dist_omega = sumblk('omega_dist = omega + d_omega');
dist_z = sumblk('z_dist = z + d_z');

err_omega = sumblk('e_omega = -omega_dist');
err_z = sumblk('e_z =- z_dist');

mimo_cl = connect(G, K, ...
    err_omega, err_z, dist_omega, dist_z, ...
    W_u, W_p,...
    {"d_omega", "d_z"}, ...
    {"z1", "z2"});

%%
opt = hinfstructOptions('Display', 'final', "RandomStart", 10);
N = hinfstruct(mimo_cl, opt);

%% recover sensitivity
% K = minreal([tf(beta_controller), 0; tf(tau_controller), 0]);
% S = connect(G, Gd, ...
%     beta_controller, tau_controller, ...
%     err_omega, err_z, dist_omega, dist_z, ...
%     W_u, W_p,...
%     {"d_omega", "d_z"}, ...
%     {"omega_dist", "z_dist"} ...
%     );

K_tuned = tf([N.Blocks.beta_omega_LL + N.Blocks.beta_omega_PID, 0; N.Blocks.tau_omega_PID + N.Blocks.tau_omega_LL, 0]);
K_tuned.InputName = ["e_omega", "e_z"];
K_tuned.OutputName = ["beta", "tau"];
    

L=minreal(G*K_tuned);
S = minreal(feedback(eye(2), L));

% S = connect(G, K, ...
%     err_omega, err_z, dist_omega, dist_z,...
%     {"d_omega", "d_z"}, ...
%     {"omega_dist", "z_dist"});
S.InputName = ["d_{\omega}", "d_z"];
S.OutputName = ["\omega{dist}", "z_{dist}"];
T = G*K_tuned*S;
KS = K_tuned*S;
hinfnorm(W_p *S)
opts = bodeoptions;
opts.PhaseVisible = "off";

L_MS = minreal(G*K_MS);
S_MS = minreal(feedback(eye(2), L_MS));
S_MS.InputName = ["d_{\omega}", "d_z"];
S_MS.OutputName = ["\omega{dist}", "z_{dist}"];
KS_MS = K_MS * S_MS;
T_MS = eye(2)-S_MS;

plot_sensitivity = false;
plot_controller_sensitivity = false;
plot_complementary_sensitivity = false;
plot_time_simulations = false;

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

if plot_complementary_sensitivity
    figure;
    hold on
    grid
    bodeplot(T(1, 1), opts)
    bodeplot(T_MS(1, 1), opts)
    legend("T_{FS}", "T_{MS}");
    title("")
    hold off
    
    figure;
    hold on
    grid
    bodeplot(T(2, 1), opts)
    bodeplot(T_MS(2, 1), opts)
    legend("T_{FS}", "T_{MS}");
    title("")
    hold off
    
    figure;
    hold on
    grid
    bodeplot(T(1, 2), opts)
    bodeplot(T_MS(1, 2), opts)
    legend("T_{FS}", "T_{MS}");
    title("")
    hold off
    
    figure;
    hold on
    grid
    bodeplot(T(2, 2), opts)
    bodeplot(T_MS(2, 2), opts)
    legend("T_{FS}", "T_{MS}");
    title("")
    hold off
end

if plot_controller_sensitivity
    figure;
    hold on
    grid
    bodeplot(KS(1, 1), opts)
    bodeplot(KS_MS(1, 1), opts)
    bodeplot(1/W_u(1, 1), opts)
    legend("S", "1/W_{u11}");
    title("KS_{fixed}", "KS_{mixed}", "W_{u11}")
    hold off
    
    figure;
    hold on
    grid
    bodeplot(KS(2, 1), opts)
    bodeplot(KS_MS(2, 1), opts)
    bodeplot(1/W_u(2, 2), opts)
    legend("S", "1/W_{u22}");
    title("KS_{fixed}", "KS_{mixed}", "W_{u22}")
    hold off
    
    figure;
    hold on
    grid
    bodeplot(KS(1, 2), opts)
    bodeplot(KS_MS(1, 2), opts)
    bodeplot(1/W_u(1, 1), opts)
    legend("S", "1/W_{u11}");
    title("KS_{fixed}", "KS_{mixed}", "W_{u11}")
    hold off
    
    figure;
    hold on
    grid
    bodeplot(KS(2, 2), opts)
    bodeplot(KS_MS(2, 2), opts)
    bodeplot(1/W_u(2, 2), opts)
    legend("S", "1/W_{u22}");
    title("KS_{fixed}", "KS_{mixed}", "W_{u22}")
    hold off
end

if plot_time_simulations
    stepWind(FWT, - K_tuned);
end