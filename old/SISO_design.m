
% States: omega, z1_dot, z1, z2_dot, z2
% Inputs: beta, tau(generator torque), wind_speed
mimo = ss(A, B, C, D);

siso = tf(mimo);
siso = minreal(-siso(1,1));
zero(mimo)
pole(mimo)
%step(siso)
%bode(siso);
%pzplot(siso)
%margin(siso)
controlSystemDesigner('bode', siso)
pzmap(mimo)
grid
