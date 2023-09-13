P.L = 0.89785;
P.alpha1 = 30 * pi/180;
P.alpha2 = 14.64 * pi/180;
P.MD2RCM = 0.57922;
P.MOA = deg2rad(30);
P.a = 20 * 1e-3;
P.b = 9.85 * 1e-3;
P.c = 545.17 * 1e-3;
P.virtual_hold_factor = 0;

q_RCM = [P.L*cos(P.alpha1); 0; -P.L*sin(P.alpha1)];
q_H = q_RCM + [0; 0; P.MD2RCM];
q_ROT2 = q_H + [0; 0; -P.c];
q_ROT3 = q_ROT2 + [0; 0; -P.b];
q_T = q_ROT3 + [0; 0; -P.a * P.virtual_hold_factor];
q_end = q_ROT3 + [0; 0; -P.a];

s = zeros(6, 7);
w = zeros(3, 7);
q = zeros(3, 7);
g0 = zeros(4,4,7);

w(:, 1) = [cos(P.alpha1); 0; -sin(P.alpha1)];
q(:, 1) = q_RCM;
g0(:,:,1) = trvec2tform(q(:, 1)'); 
s(:, 1) = [w(:, 1); -cross(w(:, 1), q(:, 1))];

w(:, 2) = [-cos(P.alpha1 + P.alpha2); 0; sin(P.alpha1 + P.alpha2)];
q(:, 2) = q_RCM;
g0(:,:,2) = trvec2tform(q(:, 2)');
s(:, 2) = [w(:, 2); -cross(w(:, 2), q(:, 2))];


w(:, 3) = [0; 1; 0];
q(:, 3) = q_RCM;
g0(:,:,3) = trvec2tform(q(:, 3)');
s(:, 3) = [w(:, 3); -cross(w(:, 3), q(:, 3))];


w(:, 4) = [0; 0; 0];
q(:, 4) = q_H;
g0(:,:,4) = trvec2tform(q(:, 4)');
s(:, 4) = [w(:, 4); [0; 0; -1]];

w(:, 5) = [0; 0; -1];
q(:, 5) = q_ROT2;
g0(:,:,5) = trvec2tform(q(:, 5)');
s(:, 5) = [w(:, 5); -cross(w(:, 5), q(:, 5))];


w(:, 6) = [0; 1; 0];
q(:, 6) = q_ROT2;
g0(:,:,6) = trvec2tform(q(:, 6)');
s(:, 6) = [w(:, 6); -cross(w(:, 6), q(:, 6))];


w(:, 7) = [1; 0; 0];
q(:, 7) = q_ROT3;
R = [0 0 -1;
     0 1 0;
     1 0 0];
M = [[R, q_T]; [0 0 0 1]];
g0(:,:,7) = M;
s(:, 7) = [w(:, 7); -cross(w(:, 7), q(:, 7))];
