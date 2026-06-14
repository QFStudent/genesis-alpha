%% an IR

lengthOfIR = 1500;
order = 10;

ir = zeros(1, 1, lengthOfIR);
for i=1:lengthOfIR;
    ir(1, 1, i) = i^(-1.1);
end

[rw, ru, rC] = msvdreduce(ir, order);
rir = mtibir(rw, ru, rC, lengthOfIR);

norm(reshape(rir, 1, []) - reshape(ir, 1, [])) / norm(reshape(ir, 1, []))