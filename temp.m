p = [1 0 0]';       
q = [0 1 0]';       
r = [0 0 1]';       % 轴上一点
w = [0 0 1]';       % 轴

[theta, success] = PKsub_R(p, q, r, w);