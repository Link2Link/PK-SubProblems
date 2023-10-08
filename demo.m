clear
clc
rng(123);

p = [1 0 0]';
w1 = rand(3,1);
w1 = w1/norm(w1);

w2 = rand(3,1);
w2 = w2/norm(w2);

w3 = rand(3,1);
w3 = w3/norm(w3);

r1 = rand(3,1);
r2 = rand(3,1);
r3 = rand(3,1);

s1 = [w1 ; -cross(w1, r1)];
s2 = [w2 ; -cross(w2, r2)];
s3 = [w3 ; -cross(w3, r3)];

theta = [0.1;0.2;0.3];

EXP = @(theta) expm(VecTose3(s1)*theta(1)) * expm(VecTose3(s2)*theta(2)) * expm(VecTose3(s3)*theta(3));
T = EXP(theta);
q_coor = T*[p; 1];
q = q_coor(1:3);
[theta1, theta2, theta3] = PKsub_Q5(p,q,w1,w2,w3,r1,r2,r3);
[theta1;theta2;theta3]
