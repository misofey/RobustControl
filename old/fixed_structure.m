matfile("Assignment_Data_SC42145_2025.mat");
s = tf('s');

Asiso = A;
Bsiso = B(:, 1);
Csiso = C(1, :);
Dsiso = D(1, 1);

siso_sys = ss(Asiso, Bsiso, Csiso, Dsiso);
Wp_simple = 0.95*(s+0.02*2*pi)/(0.0001*2*pi+s);
Wp_siso = Wp_simple;
Wp_siso.u = 'y_act';
Wp_siso.y = 'z1';

G_siso = -1* siso_sys;
G_siso.u = 'u';
G_siso.y = 'y_plant';
    
C_struct = 0.3410 + 1.7475/s + (1.0698*s)/(1.4856*s+1);
C_siso = C_struct;
C_siso.u = 'e';
C_siso.y = 'u';

Sum1 = sumblk('e=w_ref-y_act');
Sum2 = sumblk('y_act=y_plant+y_dist');

og_system = connect(G_siso, Wp_siso, C_siso, Sum1, Sum2,{'w_ref', 'y_dist'},{'y_plant' , 'z1'});
Gsiso = tf(siso_sys);
step(og_system(1, 1))
%% Fixed Structure SISO controller
Kp = realp ('Kp' ,1);
Ki = realp ('Ki' ,1);
Kd = realp ('Kd' ,1);
Tf = realp ('Tf' ,1);

Wp_simple = 0.95*(s+0.02*2*pi)/(0.0001*2*pi+s);
C_struct = Kp + Ki/s + (Kd*s)/(Tf*s+1);

% SISO hinfstruct

Wp_siso = Wp_simple;
Wp_siso.u = 'y_act';
Wp_siso.y = 'z1';

G_siso = -1* Gsiso;
G_siso.u = 'u';
G_siso.y = 'y_plant';

C_siso = C_struct;
C_siso.u = 'e';
C_siso.y = 'u';

Sum1 = sumblk('e=w_ref-y_act');
Sum2 = sumblk('y_act=y_plant+y_dist');
Siso_Con1 = connect (G_siso, Wp_siso, C_siso, Sum1, Sum2,{'w_ref', 'y_dist'},{'y_plant' , 'z1'});

opt = hinfstructOptions('Display', 'final', 'RandomStart',5);
N_siso = hinfstruct(Siso_Con1, opt);

% Extract controller gains :
Kp_opt = N_siso.Blocks.Kp.Value;
Ki_opt = N_siso.Blocks.Ki.Value;
Kd_opt = N_siso.Blocks.Kd.Value;
Tf_opt = N_siso.Blocks.Tf.Value;

Kfb_opt = Kp_opt + Ki_opt/s + (Kd_opt*s)/(Tf_opt*s+1);
Kfb_opt.u = 'e';
Kfb_opt.y = 'u';
final1 = connect (G_siso, Wp_siso, Kfb_opt, Sum1, Sum2,{'w_ref', 'y_dist'},{'y_plant' , 'z1'});
%% Custom structure controller

Kp = realp('Kp' ,1);
Ki = realp ('Ki' ,1);
Kd = realp ('Kd' ,1);
Tf = realp ('Tf' ,1);
p1 = realp("p1", 1);
p2 = realp("p2", 1);
z1 = realp("z1", 1);
z2 = realp("z2", 1);


filt = (s-z1)*(s-z2) / (s-p1)/(s-p2);
% C_struct = filt * (Kd*s)/(Tf*s+1);
C_struct = filt*Kp;

% SISO hinfstruct

Wp_siso = Wp_simple;
Wp_siso.u = 'y_act';
Wp_siso.y = 'z1';

G_siso = -1* Gsiso;
G_siso.u = 'u';
G_siso.y = 'y_plant';

C_siso = C_struct;
C_siso.u = 'e';
C_siso.y = 'u';

Sum1 = sumblk('e=w_ref-y_act');
Sum2 = sumblk('y_act=y_plant+y_dist');
Siso_Con2 = connect (G_siso, Wp_siso, C_siso, Sum1, Sum2,{'w_ref', 'y_dist'},{'y_plant' , 'z1'});

opt = hinfstructOptions('Display', 'final', 'RandomStart',5);
N_siso = hinfstruct(Siso_Con2, opt);

% Extract controller gains :
Kp_opt = N_siso.Blocks.Kp.Value;
% Ki_opt = N_siso.Blocks.Ki.Value;
% Kd_opt = N_siso.Blocks.Kd.Value;
% Tf_opt = N_siso.Blocks.Tf.Value;
p1_opt = N_siso.Blocks.p1.Value;
p2_opt = N_siso.Blocks.p2.Value;
z1_opt = N_siso.Blocks.z1.Value;
z2_opt = N_siso.Blocks.z2.Value;


filt_opt = (s-z1_opt)*(s-z2_opt) / (s-p1_opt)/(s-p2_opt);
% custom_opt = filt_opt * (Kd_opt*s)/(Tf_opt*s+1);
custom_opt = filt_opt * Kd;
custom_opt.u = 'e';
custom_opt.y = 'u';
final2 = connect (G_siso, Wp_siso, custom_opt, Sum1, Sum2,{'w_ref', 'y_dist'},{'y_plant' , 'z1'});

%% Tunable state space controller

% blk = tunableTF('tunableTF', 2, 3);
blk = tunablePID('tunableTF', 'pid');
Kp = realp('Kp' ,1);


C_struct = blk;

% SISO hinfstruct

Wp_siso = Wp_simple;
Wp_siso.u = 'y_act';
Wp_siso.y = 'z1';

G_siso = -1* Gsiso;
G_siso.u = 'u';
G_siso.y = 'y_plant';

C_siso = C_struct;
C_siso.u = 'e';
C_siso.y = 'u';

Sum1 = sumblk('e=w_ref-y_act');
Sum2 = sumblk('y_act=y_plant+y_dist');
Siso_Con2 = connect (G_siso, Wp_siso, C_siso, Sum1, Sum2,{'w_ref', 'y_dist'},{'y_plant' , 'z1'});

opt = hinfstructOptions('Display', 'final', 'RandomStart',5);
N_siso = hinfstruct(Siso_Con2, opt);

% Extract controller gains :
blk_opt = N_siso.Blocks.tunableTF;
% Kp_opt = N_siso.Blocks.Kp.Value;



custom_opt = blk_opt;
custom_opt.u = 'e';
custom_opt.y = 'u';
final3 = connect (G_siso, Wp_siso, custom_opt, Sum1, Sum2,{'w_ref', 'y_dist'},{'y_plant' , 'z1'});

%% SISO system performance
stepinfo(final1(1, 1))
stepinfo(final2(1, 1))
stepinfo(final3(1, 1))

step(final1, final2, final3)

hinfnorm(final1)
hinfnorm(final2)
hinfnorm(final3)


% controlSystemDesigner(G_siso, C_struct)
%% MIMO System

%% MIMO System performance

