clear
clc

rt = [];

result = runtests('test_PKsub_Q1');
rt = [rt; table(result)];

result = runtests('test_PKsub_Q2');
rt = [rt; table(result)];

result = runtests('test_PKsub_PRR');
rt = [rt; table(result)];
