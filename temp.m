clear
clc

config;
theta = [0.1;0.2;0.3;0.4;0.5;0.6;0.7];
T = EXP(theta, s) * M;


gst = T;
theta1 = theta(1);
T1 = expm(VecTose3(s(:, 1))*theta1);
g1 = TransInv(T1)*gst*TransInv(M);

q = q_RCM;
p = TransInv(g1)*[q;1];
p = p(1:3);

% R2R3P4eR5R6R7 交换为 R2R3R5P4R6R7 求解467
v4 = s(4:6, 4);
w6 = s(1:3, 6);
w7 = s(1:3, 7);
p6 = q_ROT2;
p7 = q_ROT3;

[thetalist, success] = PKsub_PRR(p, q, v4,w6,w7,p6,p7);
theta4 = thetalist(1,end);
theta6 = thetalist(1,end);
theta7 = thetalist(1,end);

% R2R3R5 求解235







function T = EXP(theta, s)
T = eye(4);
for k = 1:7
    T = T*expm(VecTose3(s(:, k))*theta(k));
end


end