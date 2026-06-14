%% an IR

lengthOfIR = 1500;
order = 10;

ir = zeros(3, 3, lengthOfIR);
for i=1:1500;
    ir(1, 1, i) = 0.1^i;
    ir(1, 2, i) = 0.2^i;
    ir(1, 3, i) = 0.3^i;
    ir(2, 1, i) = 0.4^i;
    ir(2, 2, i) = 0.5^i;
    ir(2, 3, i) = 0.6^i;
    ir(3, 1, i) = 0.7^i;
    ir(3, 2, i) = 0.8^i;
    ir(3, 3, i) = 0.9^i;
end


[rw, ru, rC] = msvdreduce(ir, order);
rir = mtibir(rw, ru, rC, lengthOfIR);

norm(reshape(rir, 1, []) - reshape(ir, 1, [])) / norm(reshape(ir, 1, []))