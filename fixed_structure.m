load('Assignment_Data_SC42145_2025.mat');
load("K_MS");
mimo = FWT(1:2, 1:2);
G = tf(mimo);
Gd = FWT(1:2, 3);
G = [1 0; 0 500] * G;

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


W_p.u = ["omega_dist", "z_dist"];
W_p.y = "z1";
W_u.u = ["beta", "tau"];
W_u.y = "z2";

%% interconnection

tf_order = 8;

beta_controller = tunableTF('K_beta', 4, 4) + tunablePID("K_beta2", "PD");
% beta_controller = tunablePID("K_beta", "PD");
% beta_controller.u = 'e_omega';
% beta_controller.y = 'beta';

tau_controller = tunableTF('K_tau', 4, 4);
% tau_controller = tunablePID("K_tau", "PD") + tunableTF("tftau", 2, 2);
% tau_controller.u = 'e_omega';
% tau_controller.y = 'tau';

beta_z = tunablePID("gain_betaz", "PID");
% beta_z.u = 'e_omega';
% beta_z.y = 'beta';

tau_z = tunablePID("gain_tauz", "P");
% tau_z = tunableTF('K_beta', 4, 4);
% tau_z.u = 'e_z';
% tau_z.y = 'tau';

K = [beta_controller, beta_z; tau_controller, 0];
K.InputName = ["e_omega", "e_z"];
K.OutputName = ["beta", "tau"];


%%
% blk = tunablePID('tunableTF', 'pid');
% Kp = realp('Kp' ,1);


err_omega = sumblk('e_omega = w_omega - omega_dist');
err_z = sumblk('e_z = w_z - z_dist');
dist_omega = sumblk('omega_dist = omega + d_omega');
dist_z = sumblk('z_dist = z + d_z');

% mimo_complete = connect(G, Gd,...
%     beta_controller, tau_controller, ...
%     err_omega, err_z, dist_omega, dist_z, ...
%     W_u, W_p,...
%     {"w_omega", "w_z", "d_omega", "d_z"}, ...
%     {"z1", "z2", "omega_dist", "z_dist", "beta", "tau"});
mimo_cl = connect(G, ...
    K, ...
    err_omega, err_z, dist_omega, dist_z, ...
    W_u, W_p,...
    {"w_omega", "w_z"}, ...
    {"z1", "z2"});

%%
opt = hinfstructOptions('Display', 'final', "RandomStart", 5);
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

    
L=minreal(G*K);
S = minreal(feedback(eye(2),L));
S.InputName = ["d_omega", "d_z"];
S.OutputName = ["omega_dist", "z_dist"];
T = G*K*S;
KS = K*S;
hinfnorm(S)
opts = bodeoptions;
opts.PhaseVisible = "off";

plot_sensitivity = true;
plot_controller_sensitivity = false;
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

if plot_controller_sensitivity
    figure;
    hold on
    grid
    bodeplot(KS(1, 1), opts)
    bodeplot(1/W_u(1, 1), opts)
    legend("S", "1/W_{u11}");
    title("")
    hold off
    
    figure;
    hold on
    grid
    bodeplot(KS(2, 1), opts)
    bodeplot(1/W_u(2, 2), opts)
    legend("S", "1/W_{u22}");
    title("")
    hold off
    
    figure;
    hold on
    grid
    bodeplot(KS(1, 2), opts)
    bodeplot(1/W_u(1, 1), opts)
    legend("S", "1/W_{u11}");
    title("")
    hold off
    
    figure;
    hold on
    grid
    bodeplot(KS(2, 2), opts)
    bodeplot(1/W_u(2, 2), opts)
    legend("S", "1/W_{u22}");
    title("")
    hold off
end

if plot_time_simulations
    stepWind(FWT, -K)
end