load('Assignment_Data_SC42145_2025.mat');
load("K_MS");
mimo = FWT(1:2, 1:2);
G = tf(mimo);
Gd = FWT(1:2, 3);

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

%% Uncertainity weights

W_i1 = (1/16/pi*s + 0.2) / (1/64/pi*s + 1);
W_i2 = W_i1;

W_o1 = (0.05*s + 0.25) / (0.1*s + 1);
W_o2 = W_o1;

W_i = [W_i1 0; 0 W_i2];
W_o = [W_o1 0; 0 W_o2];

plot_weights = false;

if plot_weights
    figure;
    bodemag(W_i1);
    title('');

    figure;
    bodemag(W_o1);
    title('');
end

delta_i = [ultidyn("delta_i1", 1) 0; 0 ultidyn("delta_i2", 1)];
delta_o = [ultidyn("delta_o1", 1) 0; 0 ultidyn("delta_o2", 1)];

G_p = (eye(2)+W_o*delta_o) * G * (eye(2)+W_i*delta_i);

plot_singular_values = true;
[sv, wout] = sigma(G_p);

if plot_singular_values
    figure;
    sigma(G_p);
    title('');

    figure;
    loglog(wout, sv(1, :)./sv(2, :));
    title('')
    xlabel("Frequency(rad/s)");
    ylabel('Condition number')
end

%% generalized plant

P = [zeros(2, 6) W_i;
    W_o * G zeros(2, 4) W_o*G;
    zeros(2, 6) W_u;
    W_p*G W_p W_p W_p*G;
    -G -eye(2) -eye(2) -G];

N = minreal(lft(P, K_MS));

M = N(1:4, 1:4);

%% NS, NP, RS, RP


% NS, NP
% negative controller when not using ge
L_cert = minreal(G * K_MS);

a = minreal(eye(2) + L_cert);
unstable_open_loop_zeros = sum(real(pole(L_cert))>=0)


omega = logspace(-4, 4, 300);
% NP
blk = [2 4];
[bounds_NP, muinfo_NP] = mussv(frd(N(5:8, 5:6), omega), blk, 'as');
[mag_NP, ~, ~] = bode(bounds_NP(1, 1));
maximum_ssv_NP = max(mag_NP(:))
if max(mag_NP)<1
    nominal_performance = true
else
    nominal_performance = false
end



% RS

blk = [4, 0];
[bounds_RS, muinfo_RS] = mussv(frd(M, omega), blk, 's');

[mag_RS, ~, ~] = bode(bounds_RS(1, 1));
if max(mag_RS)<1
    robust_stability = true
else
    robust_stability = false
end

% RP
blk = [[4 0]; [2 4]];
[bounds_RP] = mussv(frd(N, omega), blk, 's');

[mag_RP, ~, ~] = bode(bounds_RP(1, 1));
if max(mag_RP) < 1
    robust_performance = true;
else
    robust_performance = false;
end

% NP
blk_NP = [2 4];
mu_NP_max = peak_mu(N(5:8, 5:6), blk_NP, omega, 'as');

% RS
blk_RS = [4 0];
mu_RS_max = peak_mu(M, blk_RS, omega, 's');

% RP
blk_RP = [[4 0]; [2 4]];
mu_RP_max = peak_mu(N, blk_RP, omega, 's');

fprintf('\n--- μ PEAK VALUES (Before D-K) ---\n');
fprintf('Nominal Performance μ_max = %.4f\n', mu_NP_max);
fprintf('Robust Stability   μ_max = %.4f\n', mu_RS_max);
fprintf('Robust Performance μ_max = %.4f\n', mu_RP_max);

plot_mu_analysis = true;

if plot_mu_analysis
    figure;
    nyquist(minreal((a(1,1)*a(2,2) - a(2,1)*a(1,2))));
    title('');

    figure;
    semilogx(omega, mag2db(mag_NP(:)))
    hold on
    semilogx(omega, mag2db(mag_RS(:)))
    semilogx(omega, mag2db(mag_RP(:)))
    title('');
    legend("NP", "RS", "RP");
    grid()
    hold off
end


%% D-K iterations

NDelta = lft(blkdiag(delta_i, delta_o), P);

[K_DK, CLPerf] = musyn(NDelta, 2, 2);

 %% NS, NP, RS, RP
N_DK = minreal(lft(P, K_DK));
M_DK = N_DK(1:4, 1:4);
% NS, NP
% negative controller when not using ge
L_DK = minreal(G * K_DK);

a = minreal(eye(2) + L_DK);
unstable_open_loop_zeros = sum(real(pole(L_DK))>=0)


omega = logspace(-4, 4, 300);
% NP
blk = [2 4];
[bounds_NP, muinfo_NP] = mussv(frd(N_DK(5:8, 5:6), omega), blk, 'as');
[mag_NP, ~, ~] = bode(bounds_NP(1, 1));
maximum_ssv_NP = max(mag_NP(:))
if max(mag_NP)<1
    nominal_performance = true
else
    nominal_performance = false
end



% RS

blk = [4, 0];
[bounds_RS, muinfo_RS] = mussv(frd(M, omega), blk, 's');

[mag_RS, ~, ~] = bode(bounds_RS(1, 1));
if max(mag_RS)<1
    robust_stability = true
else
    robust_stability = false
end

% RP
blk = [[4 0]; [2 4]];
[bounds_RP] = mussv(frd(N_DK, omega), blk, 's');

[mag_RP, ~, ~] = bode(bounds_RP(1, 1));
if max(mag_RP) < 1
    robust_performance = true;
else
    robust_performance = false;
end

% NP
mu_NP_DK_max = peak_mu(N_DK(5:8, 5:6), blk_NP, omega, 'as');

% RS
mu_RS_DK_max = peak_mu(M_DK, blk_RS, omega, 's');

% RP
mu_RP_DK_max = peak_mu(N_DK, blk_RP, omega, 's');

fprintf('\n--- μ PEAK VALUES (After D-K) ---\n');
fprintf('Nominal Performance μ_max = %.4f\n', mu_NP_DK_max);
fprintf('Robust Stability   μ_max = %.4f\n', mu_RS_DK_max);
fprintf('Robust Performance μ_max = %.4f\n', mu_RP_DK_max);


plot_mu_analysis = true;

if plot_mu_analysis
    figure;
    nyquist(minreal((a(1,1)*a(2,2) - a(2,1)*a(1,2))));
    title('');

    figure;
    semilogx(omega, mag2db(mag_NP(:)))
    hold on
    semilogx(omega, mag2db(mag_RS(:)))
    semilogx(omega, mag2db(mag_RP(:)))
    title('');
    legend("NP", "RS", "RP");
    grid()
    hold off
end

function mu_max = peak_mu(sys, blk, omega, mu_type)
    % Compute peak structured singular value over frequency
    bounds = mussv(frd(sys, omega), blk, mu_type);
    mag = squeeze(abs(freqresp(bounds(1,1), omega)));
    mu_max = max(mag);
end

%% TIME-DOMAIN SIMULATIONS

figure('Name', 'Nominal Response Comparison', 'Position', [100 100 1200 800]);

% Mixed-Sensitivity Controller
subplot(2,2,1);
%f1 = figure;
stepWind(FWT, -K_MS);
title('Mixed-Sensitivity: Step Wind Response');
%exportgraphics(f1, "/home/mahargardr/master_courses/robust_control/part_2/images/ms_stepwind_response.png", 'Resolution',300);

% D-K Controller
%f2 = figure;
subplot(2,2,2);
stepWind(FWT, -K_DK);
title('D-K Iteration: Step Wind Response');
%exportgraphics(f2, "/home/mahargardr/master_courses/robust_control/part_2/images/dk_stepwind_response.png", 'Resolution',300);

%omega_r
subplot(2,2,3);
G_sim = G; G_sim.u = 'u'; G_sim.y = 'y_p';
Gd_sim = Gd; Gd_sim.u = 'V'; Gd_sim.y = 'y_d';
Sum_sim = sumblk('y = y_p + y_d', 2);

K_MS_sim = -K_MS; K_MS_sim.u = 'y'; K_MS_sim.y = 'u';
CL_MS = connect(G_sim, K_MS_sim, Gd_sim, Sum_sim, 'V', 'y');

K_DK_sim = -K_DK; K_DK_sim.u = 'y'; K_DK_sim.y = 'u';
CL_DK = connect(G_sim, K_DK_sim, Gd_sim, Sum_sim, 'V', 'y');

t = 0:0.01:100;
[y_MS, ~] = step(CL_MS, t);
[y_DK, ~] = step(CL_DK, t);

plot(t, y_MS(:,1), 'b-', 'LineWidth', 1.5); hold on;
plot(t, y_DK(:,1), 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('\omega_r');
title('Generator Speed Comparison');
legend('MS', 'D-K'); grid on;

subplot(2,2,4);
plot(t, y_MS(:,2), 'b-', 'LineWidth', 1.5); hold on;
plot(t, y_DK(:,2), 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('z (m)');
title('Platform Displacement Comparison');
legend('MS', 'D-K'); grid on;

sgtitle('Nominal Step Wind Disturbance Response');

%% ROBUST SIMULATIONS

num_samples = 20;

figure('Name', 'Uncertain Response', 'Position', [100 100 1400 600]);

% MS Robustness
subplot(1,2,1); hold on;
for i = 1:num_samples
    G_s = usample(G_p);
    G_s.u = 'u'; G_s.y = 'y_p';
    try
        CL_s = connect(G_s, K_MS_sim, Gd_sim, Sum_sim, 'V', 'y');
        [y_s, ~] = step(CL_s, t);
        plot(t, y_s(:,1), 'b-', 'Color', [0 0 1 0.2]);
    catch; end
end
plot(t, y_MS(:,1), 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\omega_r (rad/s)');
title('MS Controller Under Uncertainty'); grid on;

% D-K Robustness
subplot(1,2,2); hold on;
for i = 1:num_samples
    G_s = usample(G_p);
    G_s.u = 'u'; G_s.y = 'y_p';
    try
        CL_s = connect(G_s, K_DK_sim, Gd_sim, Sum_sim, 'V', 'y');
        [y_s, ~] = step(CL_s, t);
        plot(t, y_s(:,1), 'r-', 'Color', [1 0 0 0.2]);
    catch; end
end
plot(t, y_DK(:,1), 'r-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\omega_r (rad/s)');
title('D-K Controller Under Uncertainty'); grid on;

sgtitle('Step Response Under Uncertainty (20 samples)');

%% 3. FREQUENCY-DOMAIN COMPARISON

L_MS = minreal(G * K_MS);
L_DK = minreal(G * K_DK);
S_MS = minreal(feedback(eye(2), L_MS));
S_DK = minreal(feedback(eye(2), L_DK));
T_MS = eye(2) - S_MS;
T_DK = eye(2) - S_DK;

opts = bodeoptions; opts.PhaseVisible = 'off'; opts.Grid = 'on';

figure('Name', 'Sensitivity Comparison');
bodemag(S_MS(1,1), 'b-', S_DK(1,1), 'r--', 1/W_p1, 'k:', opts);
legend('S_{MS}', 'S_{DK}', '1/W_p', 'Location', 'best');
title('Sensitivity S(1,1): d_\omega \rightarrow \omega');

figure('Name', 'Complementary Sensitivity Comparison');
bodemag(T_MS(1,1), 'b-', T_DK(1,1), 'r--', opts);
legend('T_{MS}', 'T_{DK}', 'Location', 'best');
title('Complementary Sensitivity T(1,1)');
