%% an IR

lengthOfIR = 1000;
order = 20;

ir = zeros(3, 3, lengthOfIR);
for i=1:500;
    ir(1, 1, i) = i^(-1.1);
    ir(1, 2, i) = i^(-1.4);
    ir(1, 3, i) = i^(-1.7);
    ir(2, 1, i) = i^(-1.2);
    ir(2, 2, i) = i^(-1.5);
    ir(2, 3, i) = i^(-1.8);
    ir(3, 1, i) = i^(-1.3);
    ir(3, 2, i) = i^(-1.6);
    ir(3, 3, i) = i^(-1.9);
end

[rw, ru, rC] = msvdreduce(ir, order);
rir = mtibir(rw, ru, rC, lengthOfIR);

norm(reshape(rir, 1, []) - reshape(ir, 1, [])) / norm(reshape(ir, 1, []))
