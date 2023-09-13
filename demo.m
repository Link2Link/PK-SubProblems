clear
clc
config;


%% FK 
thetalist = zeros(7,1);
T = FKinSpace(M, s, thetalist);




theta1 = thetalist(1);
Tst = T;





