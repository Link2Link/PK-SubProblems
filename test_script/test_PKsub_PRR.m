% test PKsub_PRR
runTimes = 1e3;
tol = 1e-8;

%% Test 1 : accurate
for k = 1:runTimes
    theta = (rand(3,1) - 0.5)*2;

    v1 = [0;0;1];
    w2 = [1;0;0];
    w3 = [0;1;0];
    
    p2 = [0;0;0.1];
    p3 = [0;0;-0.1];
    
    s1 = [zeros(3,1); v1];
    s2 = [w2; -cross(w2, p2)];
    s3 = [w3; -cross(w3, p3)];
    
    p = [1;2;3];
    p_coor = [p;1];
    
    EXP = @(theta) expm(VecTose3(s1)*theta(1))*expm(VecTose3(s2)*theta(2))*expm(VecTose3(s3)*theta(3));
    q_coor = EXP(theta) * p_coor;
    q = q_coor(1:3);
    
    [thetalist, valid] = PKsub_PRR(p, q, v1,w2,w3,p2,p3);
    
    q_coor1 = EXP(thetalist(:, 1)) * p_coor;
    q1 = q_coor1(1:3);
    
    q_coor2 = EXP(thetalist(:, 2)) * p_coor;
    q2 = q_coor2(1:3);
    
    q_coor3 = EXP(thetalist(:, 3)) * p_coor;
    q3 = q_coor3(1:3);
    
    q_coor4 = EXP(thetalist(:, 4)) * p_coor;
    q4 = q_coor4(1:3);

    assert(norm(q - q1) < tol, 'result not consistence')
    assert(norm(q - q2) < tol, 'result not consistence')
    assert(norm(q - q3) < tol, 'result not consistence')
    assert(norm(q - q4) < tol, 'result not consistence')
end


%% Test 2 : random param 
for k = 1:runTimes
    theta = (rand(3,1) - 0.5)*2;

    v1 = [0;0;1];
    w2 = [1;0;0];
    w3 = [0;1;0];
    
    p2 = [0;0;rand()];
    p3 = [0;0;rand()];
    
    s1 = [zeros(3,1); v1];
    s2 = [w2; -cross(w2, p2)];
    s3 = [w3; -cross(w3, p3)];
    
    p = (rand(3,1) - 0.5)*2;
    p_coor = [p;1];
    
    EXP = @(theta) expm(VecTose3(s1)*theta(1))*expm(VecTose3(s2)*theta(2))*expm(VecTose3(s3)*theta(3));
    q_coor = EXP(theta) * p_coor;
    q = q_coor(1:3);
    
    [thetalist, valid] = PKsub_PRR(p, q, v1,w2,w3,p2,p3);
    
    q_coor1 = EXP(thetalist(:, 1)) * p_coor;
    q1 = q_coor1(1:3);
    
    q_coor2 = EXP(thetalist(:, 2)) * p_coor;
    q2 = q_coor2(1:3);
    
    q_coor3 = EXP(thetalist(:, 3)) * p_coor;
    q3 = q_coor3(1:3);
    
    q_coor4 = EXP(thetalist(:, 4)) * p_coor;
    q4 = q_coor4(1:3);

    if valid(1)
        assert(norm(q - q1) < tol, 'result not consistence')
    end
    if valid(2)
        assert(norm(q - q2) < tol, 'result not consistence')
    end
    if valid(3)
        assert(norm(q - q3) < tol, 'result not consistence')
    end
    if valid(4)
        assert(norm(q - q4) < tol, 'result not consistence')
    end
    
    assert(max(valid) == 1, 'no result');
    
end

%% Test 3 : random axis
for k = 1:runTimes
    theta = (rand(3,1) - 0.5)*2;

    v1 = rand(3,1);
    v1 = v1/norm(v1);
    
    p2 = v1*rand();
    p3 = v1*rand();
    w2 = rand(3,1);
    w2 = cross(v1, w2);
    w2 = w2/norm(w2);
    w3 = cross(v1, w2);
    w3 = w3/norm(w3);
    
    
    s1 = [zeros(3,1); v1];
    s2 = [w2; -cross(w2, p2)];
    s3 = [w3; -cross(w3, p3)];
    
    p = (rand(3,1) - 0.5)*2;
    p_coor = [p;1];
    
    EXP = @(theta) expm(VecTose3(s1)*theta(1))*expm(VecTose3(s2)*theta(2))*expm(VecTose3(s3)*theta(3));
    q_coor = EXP(theta) * p_coor;
    q = q_coor(1:3);
    
    [thetalist, valid] = PKsub_PRR(p, q, v1,w2,w3,p2,p3);
    
    q_coor1 = EXP(thetalist(:, 1)) * p_coor;
    q1 = q_coor1(1:3);
    
    q_coor2 = EXP(thetalist(:, 2)) * p_coor;
    q2 = q_coor2(1:3);
    
    q_coor3 = EXP(thetalist(:, 3)) * p_coor;
    q3 = q_coor3(1:3);
    
    q_coor4 = EXP(thetalist(:, 4)) * p_coor;
    q4 = q_coor4(1:3);

    if valid(1)
        assert(norm(q - q1) < tol, 'result not consistence')
    end
    if valid(2)
        assert(norm(q - q2) < tol, 'result not consistence')
    end
    if valid(3)
        assert(norm(q - q3) < tol, 'result not consistence')
    end
    if valid(4)
        assert(norm(q - q4) < tol, 'result not consistence')
    end
    
    assert(max(valid) == 1, 'no result');
    
end
