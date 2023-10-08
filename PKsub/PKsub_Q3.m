function [result1,result2,success] = PKsub_Q3(p,q,w,r,d)
% Paden-Kahan子问题 subproblem 3
% w参考轴，r为轴上一点,从p转到距q点d的位置，所需角度theta
% 
% theta theta_为两组解

% p = [1 0 0]';
% q = [1 1 1]';
% r = [0 0 0]';
% w = [0 0 1]';
% d = 2;
% 
% [theta,theta_,success] = PKsub3(p,q,w,r,d);
% p_ = expm(VecToso3(w)*theta_)*p;
% norm(p_ - q)

p1 = p-r;
p2 = q-r;
k = w;


[theta, is_LS] = PKsub_Q4(p2, p1, k, 1/2 * (dot(p1,p1)+dot(p2,p2)-d^2));
result1 = theta(1);
result2 = theta(2);
if is_LS
    success = 0;
else
    success = 1;
end


end