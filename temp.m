clear
clc

config;

tic
for cout = 1:1e4
%       6   7   8   9   10  11  12
theta = [0.1;0.2;0.3;0.4;0.5;0.6;0.7];
T = EXP(theta, s) * M;
[T1_, T2_, T3_, T4_, T5_, T6_, T7_] = EXP2(theta, s);
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
theta6 = thetalist(2,end);
theta7 = thetalist(3,end);

%% R2R3R5 求解235
% 求解23
T4 = expm(VecTose3(s(:, 4))*theta4);
T6 = expm(VecTose3(s(:, 6))*theta6);
T7 = expm(VecTose3(s(:, 7))*theta7);
g2 = g1*TransInv(T4*T6*T7);
q_r = q_RCM + [0;0;1];
p = q_r;
q_coor = g2*[p;1];
q = q_coor(1:3);

w2 = s(1:3, 2);
w3 = s(1:3, 3);

[result1, result2, success] = PKsub_Q2(p,q,w2,w3,q_RCM);

theta2 = result1(1);
theta3 = result1(2);

% 求5
T2 = expm(VecTose3(s(:, 2))*theta2);
T3 = expm(VecTose3(s(:, 3))*theta3);
g3 = TransInv(T2*T3)*g2;

p = q_RCM + [0;1;0];
q_coor = g3 * [p;1];
q = q_coor(1:3);
w5 = s(1:3, 5);
[result, success] = PKsub_Q1(p,q,w5,q_RCM);
theta5 = result;

thetalist = [theta1, theta2, theta3, theta4, theta5, theta6, theta7];
end
toc



function T = EXP(theta, s)
T = eye(4);
for k = 1:7
    T = T*expm(VecTose3(s(:, k))*theta(k));
end
end
function [T1, T2, T3, T4, T5, T6, T7] = EXP2(theta, s)
T1=expm(VecTose3(s(:, 1))*theta(1));
T2=expm(VecTose3(s(:, 2))*theta(2));
T3=expm(VecTose3(s(:, 3))*theta(3));
T5=expm(VecTose3(s(:, 5))*theta(5));
T4=expm(VecTose3(s(:, 4))*theta(4));
T6=expm(VecTose3(s(:, 6))*theta(6));
T7=expm(VecTose3(s(:, 7))*theta(7));
end

