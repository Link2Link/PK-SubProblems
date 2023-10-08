function [theta, success] = PKsub_Q1(p,q,w,r)
% Paden-Kahan子问题 subproblem 1
% w为参考轴，r为轴上一点,从p转到q所需角度theta
%   Inputs:
%       p: 3x1 vector
%       q: 3x1 vector
%       w: 3x1 vector with norm(w) = 1
%       r: 3x1 vector
%   Outputs:
%       theta: 1x1 angle (in radians)
%       success: 1x1 boolean
% Demo:
% p = [1 0 0]';       
% q = [0 1 0]';       
% r = [0 0 1]';       % 轴上一点
% w = [0 0 1]';       % 轴
% 
% [theta, success] = PKsub_Q1(p, q, r, w);
%

tol = 1e-8;

p = p(1:3);
q = q(1:3);
r = r(1:3);

p1 = p-r;
p2 = q-r;
k = w;
KxP = cross(k, p1);
A = [KxP -cross(k,KxP)];

x = A'*p2;
theta = atan2(x(1),x(2));
if abs(norm(p1) - norm(p2)) > tol || abs(dot(k,p1) - dot(k,p2)) > tol
    success = 0;
else
    success = 1;
end

end