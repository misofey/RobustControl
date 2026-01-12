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

W_i1 = (1/16/pi*s + 0.) / (1/64/pi*s + 1);
W_i2 = W_i1;

W_o1 = (0.05*s + 0.25) / (0.01*s + 1);
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
if plot_singular_values
    figure;
    sigma(G_p);
    title('');
end

%% generalized plant

P = [zeros(2, 6) W_i;
    W_o * G zeros(2, 4) W_o*G;
    zeros(2, 6) W_u;
    W_p*G W_p zeros(2) W_p*G;
    -G -eye(2) eye(2) -G];

N = lft(P, K_MS);

M = N(1:2, 1:2);