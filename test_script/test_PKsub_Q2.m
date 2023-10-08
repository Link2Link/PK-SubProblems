% test PKsub_PRR
runTimes = 1e3;
tol = 1e-8;

%% Test 1 : accurate
for k = 1:runTimes
    theta = (rand(2,1) - 0.5)*2;
    w1 = [1;0;0];
    w2 = [0;1;0];
    r = [0;0;rand()];
    s1 = [w1; -cross(w1, r)];
    s2 = [w2; -cross(w2, r)];
    
    p = [1;2;3];
    p_coor = [p;1];
    
    EXP = @(theta) expm(VecTose3(s1)*theta(1)) * expm(VecTose3(s2)*theta(2));
    q_coor = EXP(theta) * p_coor;
    q = q_coor(1:3);
    
    [thetalist1, thetalist2, valid] = PKsub_Q2(p,q,w1,w2,r);
    
    q_coor1 = EXP(thetalist1) * p_coor;
    q1 = q_coor1(1:3);
    assert(norm(q - q1) < tol, 'result not consistence')

    q_coor2 = EXP(thetalist2) * p_coor;
    q2 = q_coor2(1:3);
    assert(norm(q - q2) < tol, 'result not consistence')

end

%% Test 2 : random axis
for k = 1:runTimes
    theta = (rand(2,1) - 0.5)*2;

    w1 = rand(3,1);
    w1 = w1/norm(w1);

    w2 = rand(3,1);
    w2 = w2/norm(w2);
    
    r = rand(3,1);
    s1 = [w1; -cross(w1, r)];
    s2 = [w2; -cross(w2, r)];
    
    p = [1;2;3];
    p_coor = [p;1];
    
    EXP = @(theta) expm(VecTose3(s1)*theta(1)) * expm(VecTose3(s2)*theta(2));
    q_coor = EXP(theta) * p_coor;
    q = q_coor(1:3);
    
    [thetalist1, thetalist2, valid] = PKsub_Q2(p,q,w1,w2,r);
    
    q_coor1 = EXP(thetalist1) * p_coor;
    q1 = q_coor1(1:3);
    assert(norm(q - q1) < tol, 'result not consistence')

    q_coor2 = EXP(thetalist2) * p_coor;
    q2 = q_coor2(1:3);
    assert(norm(q - q2) < tol, 'result not consistence')

end

