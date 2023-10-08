runTimes = 1e3;
tol = 1e-8;

%% Test 1 : accurate
for k = 1:runTimes
    theta = (rand(1,1) - 0.5)*2;

    w1 = [0;0;1];
    r = [0;0;rand()];
    s1 = [w1; -cross(w1, r)];
    
    p = [1;2;3];
    p_coor = [p;1];
    
    EXP = @(theta) expm(VecTose3(s1)*theta(1));
    q_coor = EXP(theta) * p_coor;
    q = q_coor(1:3);
    
    [thetalist, valid] = PKsub_Q1(p,q,w1,r);
    
    q_coor1 = EXP(thetalist) * p_coor;
    q1 = q_coor1(1:3);
    
    assert(norm(q - q1) < tol, 'result not consistence')

end

