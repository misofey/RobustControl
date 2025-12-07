load("/home/mahargardr/master_courses/robust_control/Assignment_Data_SC42145_2025.mat")

% make ss model
sys = ss(A,B,C,D);

% indices
u_control = 1;
u_disturb  = 3;
y = 1;

G_pitch = tf(sys(1,1));
G_pitch = G_pitch * -1;
G_pitch = minreal(G_pitch);
cl_pitch = feedback(G_pitch,1);

G_wind = tf(sys(1,3));

% Bode plot both
figure
h = bodeplot(G_pitch);
setoptions(h,'FreqUnits','Hz')  
grid on
figure

% Pole zero plot
pzplot(G_pitch);
%% 

% Simulation
load("/home/mahargardr/master_courses/robust_control/siso_controller.mat")

FWT = ss(A, B, C, D);

FWT.InputName = {'beta', 'tau_e', 'V'};
FWT.OutputName = {'omega_r', 'z'};

disp('Controller:');
siso_controller

figure;
stepWind(FWT, siso_controller);