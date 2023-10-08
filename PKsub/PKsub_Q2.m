function [result1, result2, success] = PKsub_Q2(p,q,w1,w2,r)
% Paden-Kahan子问题 subproblem 2
% w1,w2为参考轴，r为两轴交点,从p转到q所需角度theta
% theta theta_为两组解
%   Inputs:
%       p: 3x1 vector
%       q: 3x1 vector
%       w1: 3x1 vector with norm(w) = 1
%       w2: 3x1 vector with norm(w) = 1
%       r:  3x1 vector
%   Outputs:
%       result1: 2x1 angle set1 (in radians)
%       result2: 2x1 angle set2 (in radians)
%       success: 1x1 boolean
%
% Demo
% w1 = [1 0 0]';
% w2 = [0 1 0]';
% theta1 = pi/4;
% theta2 = pi/6;
% p = [1 0 0]';
% q = expm(VecToso3(w1)*theta1)*expm(VecToso3(w2)*theta2)*p;
% r = [0 0 0]';
% 
% [theta, theta_, success] = PKsub_RR(p,q,w1,w2,r);

tol = 1e-8;

p = p(1:3);
q = q(1:3);
r = r(1:3);

p2 = p-r;
p1 = q-r;
k1 = -w1;
k2 = w2;

p1_nrm = p1/norm(p1);
p2_nrm = p2/norm(p2);

[theta1, t1_is_LS] = PKsub_Q4(k2, p1_nrm, k1, dot(k2,p2_nrm));
[theta2, t2_is_LS] = PKsub_Q4(k1, p2_nrm, k2, dot(k1,p1_nrm));

theta1 = [theta1(1) theta1(end)];
theta2 = [theta2(end) theta2(1)];

is_LS = abs(norm(p1) - norm(p2)) > tol || t1_is_LS || t2_is_LS;
if is_LS
    success = 0;
else
    success = 1;
end

result1 = [theta1(1); theta2(1)];
result2 = [theta1(2); theta2(2)];


end