load('Assignment_Data_SC42145.mat');
s = tf('s');


G = ss(A, B, C, D);
G.u = {'beta', 'tau', 'V'};
G.y = {'omega', 'z'};

pid_1 = tunablePID("pid_1", 'pid');
pid_2 = tunablePID("pid_2", 'pid');
pid_3 = tunablePID("pid_3", 'pid');

K1 = pid_1;
K1.u = 'omega';
K1.y = 'tau';

sum_z = sumblk('z_sum = z + tau');
K2 = pid_2;
K2.u = 'z_sum';
K2.y = 'beta_ctrl';

K3 = pid_3;
K3.u = 'V';
K3.y = 'V_ctrl';
sum_beta = sumblk('beta = V_ctrl + beta_ctrl');

