%% Recover known state space 

% initialize
numberOfPoles = 5.;
lengthOfIR = 200;

poles = (1:numberOfPoles)'/(numberOfPoles+1);
nullvec = randn(3, numberOfPoles);
[M, N, v] = spmtib(poles, nullvec);
A = full(M \ N);
B = full(M \ v);
C = randn(3, numberOfPoles);
ir = mtibir(poles, nullvec, C, lengthOfIR);

% recover
order = numberOfPoles;
[rw, ru, rC] = msvdreduce(ir, order);
[rM, rN, rv] = spmtib(rw, ru);
rA = full(rM \ rN);
rir = mtibir(rw, ru, rC, lengthOfIR);

[sort(diag(A)), sort(diag(rA))]
norm(reshape(rir, 1, []) - reshape(ir, 1, [])) / norm(reshape(ir, 1, []))
