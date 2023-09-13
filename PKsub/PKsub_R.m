function [theta, success] = PKsub_R(p,q,r,w)
% Paden-Kahan子问题 subproblem 1
% w为参考轴，r为轴上一点,从p转到q所需角度theta
% p = [1 0 0]';       
% q = [0 1 0]';       
% r = [0 0 1]';       % 轴上一点
% w = [0 0 1]';       % 轴
% 
% [theta, success] = PKsub_R(p, q, r, w);

success = 1;
epis = 1e-8;

p = p(1:3);
q = q(1:3);
r = r(1:3);

u = p - r;
v = q - r;
u_ = u - w*w'*u;
v_ = v - w*w'*v;

err_1 = w'*u - w'*v;
err_2 = norm(u_) - norm(v_);
err = abs(err_1) + abs(err_2);
if err > epis
    success = 0;    % 无解
    theta = 0;
    return
end

theta = atan2(w'*cross(u_, v_), u_'*v_);
end