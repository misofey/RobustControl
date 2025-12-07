function [] = stepWind(FWT,K)
% STEPWIND Applies a wind step disturbance to the closed loop consisting of
% the FWT and a controller K of choice (may be SISO or MIMO).

if size(K,1) == 1 % adapt if SISO
    K = blkdiag(K,0);
end
G = FWT(:,1:2);

% Contstruct CL for simulation
G_sim = G;
G_sim.u = 'u';
G_sim.y = 'y_plant';

K_sim = K;
K_sim.u = 'y_meas';
K_sim.y = 'u';

Gd_sim = FWT(:,3);
Gd_sim.u = 'V';
Gd_sim.y = 'y_dist';

Sum_sim = sumblk('y_meas = y_plant + y_dist', 2);

P_sim = connect(G_sim,K_sim,Gd_sim,Sum_sim,{'V'},{'y_meas','u'});

P_sim.InputName = {'V (m/s)'};
P_sim.OutputName = [G.OutputName; G.InputName];

step(P_sim)

end