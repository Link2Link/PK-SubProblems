function [thetalist, success] = PKsub_PRR(p, q, v1,w2,w3,p2,p3)
% PK子问题PRR，将点p经过PRR得到q,要求三个轴两两垂直
% 给定平移轴v1，旋转轴w2，旋转轴w3, w2轴上点p2, w3轴上点p3
% Demo:
% theta = [0.1;0.2;0.3];
% 
% v1 = [0;0;1];
% w2 = [1;0;0];
% w3 = [0;1;0];
% 
% p2 = [0;0;0.1];
% p3 = [0;0;-0.1];
% 
% s1 = [zeros(3,1); v1];
% s2 = [w2; -cross(w2, p2)];
% s3 = [w3; -cross(w3, p3)];
% 
% p = [1;2;3];
% p_coor = [p;1];
% 
% EXP = @(theta) expm(VecTose3(s1)*theta(1))*expm(VecTose3(s2)*theta(2))*expm(VecTose3(s3)*theta(3));
% q_coor = EXP(theta) * p_coor;
% q = q_coor(1:3);
% 
% [thetalist, valid] = PKsub_PRR(p, q, v1,w2,w3,p2,p3);
% thetalist


success = [0,0,0,0];
flag3 = 0;
flag2_1 = 0;
flag2_2 = 0;
thetalist = zeros(3, 4);


%% 构造投影
r3 = p3 + w3*w3'*(p-p3);
u31 = p - r3;

r2 = p2 + w2*w2'*(q-p2);
u11 = q-r2;

u12 = u11 - v1*v1'*u11;

%% theta3
fc1 = w2*w2'*(r2-r3);
fr3_length_square = u31'*u31 - fc1'*fc1;
if fr3_length_square >= 0 
    fr3_length = sqrt(fr3_length_square);           %开方之前需要为正
    flag3 = 1;
else
    fr3_length = 0;
    flag3 = 0;        %无解
end

u32_1 = fc1 + fr3_length * v1;
u32_2 = fc1 - fr3_length * v1;

cos_theta3_1 = u31'*u32_1/norm(u31) / norm(u32_1);   
cos_theta3_2 = u31'*u32_2/norm(u31) / norm(u32_2);

theta3_1 = acos(cos_theta3_1);
theta3_2 = acos(cos_theta3_2);

%判断方向
if cross(u31,u32_1)'*w3 < 0
    theta3_1 = -theta3_1;
end
if cross(u31,u32_2)'*w3 < 0
    theta3_2 = -theta3_2;
end

s3 = [w3; -cross(w3, p3)];
p_coor = [p;1];
c1_coor_1 = expm(VecTose3(s3)*theta3_1)*p_coor;
c1_coor_2 = expm(VecTose3(s3)*theta3_2)*p_coor;
c1_1 = c1_coor_1(1:3);                          % c1第一个解
c1_2 = c1_coor_2(1:3);                          % c1第二个解

%% theta2 for case theta3_1
c1 = c1_1;              % 对c1第一个解求后续
u21 = c1 - r2;
nc2_length_square = u21'*u21 - u12'*u12;         % 这里需要判断 开方之前需要为正
if nc2_length_square >= 0 
    nc2_length = sqrt(nc2_length_square);
    flag2_1 = 1;
else
    nc2_length = 0;
    flag2_1 = 0;
end

u22_1 = u12 + nc2_length*v1;
u22_2 = u12 - nc2_length*v1;
cos_theta2_1 = u21'*u22_1/norm(u21) / norm(u22_1);     
cos_theta2_2 = u21'*u22_2/norm(u21) / norm(u22_2);      
theta2_1 = acos(cos_theta2_1);             
theta2_2 = acos(cos_theta2_2);

%判断方向
if cross(u21,u22_1)'*w2 < 0
    theta2_1 = -theta2_1;
end
if cross(u21,u22_2)'*w2 < 0
    theta2_2 = -theta2_2;
end

s2 = [w2; -cross(w2, p2)];
c2_coor_1 = expm(VecTose3(s2)*theta2_1)*c1_coor_1;
c2_coor_2 = expm(VecTose3(s2)*theta2_2)*c1_coor_1;
c2_1 = c2_coor_1(1:3);                          % c2第一个解
c2_2 = c2_coor_2(1:3);                          % c2第二个解

theta1_1 = (q - c2_1)'*v1;                      % 对应theta1第一个解
theta1_2 = (q - c2_2)'*v1;                      % 对应theta2第二个解

thetalist(:, 1) = [theta1_1; theta2_1; theta3_1];
thetalist(:, 2) = [theta1_2; theta2_2; theta3_1];

%% theta2 for case theta3_2
c1 = c1_2;              % 对c1第一个解求后续
u21 = c1 - r2;
nc2_length_square = u21'*u21 - u12'*u12;         % 这里需要判断 开方之前需要为正
if nc2_length_square >= 0
    nc2_length = sqrt(nc2_length_square);
    flag2_2 = 1;
else
    nc2_length = 0;
    flag2_2 = 0;
end

u22_1 = u12 + nc2_length*v1;
u22_2 = u12 - nc2_length*v1;
cos_theta2_1 = u21'*u22_1/norm(u21) / norm(u22_1);      % 除之前需要判断除0
cos_theta2_2 = u21'*u22_2/norm(u21) / norm(u22_2);      % 除之前需要判断除0
theta2_1 = acos(cos_theta2_1);              % 反三角函数需要确保[-1 1]
theta2_2 = acos(cos_theta2_2);

%判断方向
if cross(u21,u22_1)'*w2 < 0
    theta2_1 = -theta2_1;
end
if cross(u21,u22_2)'*w2 < 0
    theta2_2 = -theta2_2;
end

s2 = [w2; -cross(w2, p2)];
c2_coor_1 = expm(VecTose3(s2)*theta2_1)*c1_coor_2;
c2_coor_2 = expm(VecTose3(s2)*theta2_2)*c1_coor_2;
c2_1 = c2_coor_1(1:3);                          % c2第一个解
c2_2 = c2_coor_2(1:3);                          % c2第二个解

theta1_1 = (q - c2_1)'*v1;                      % 对应theta1第一个解
theta1_2 = (q - c2_2)'*v1;                      % 对应theta2第二个解
thetalist(:, 3) = [theta1_1; theta2_1; theta3_2];
thetalist(:, 4) = [theta1_2; theta2_2; theta3_2];

%% 判断各个解是否有效
if flag3 && flag2_1 && flag2_2
    success = [1,1,1,1];
elseif flag3 && flag2_1 && ~flag2_2
    success = [1,1,0,0];
elseif flag3 && ~flag2_1 && flag2_2
    success = [0,0,1,1];
else
    success = [0,0,0,0];
end



end

function so3mat = VecToso3(omg)
so3mat = [0, -omg(3), omg(2); omg(3), 0, -omg(1); -omg(2), omg(1), 0];
end

function se3mat = VecTose3(V)
se3mat = [VecToso3(V(1: 3)), V(4: 6); 0, 0, 0, 0];
end


