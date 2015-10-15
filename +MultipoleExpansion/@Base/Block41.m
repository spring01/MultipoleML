function res = Block41(~, vecZYX)
raz = vecZYX(1);
ray = vecZYX(2);
rax = vecZYX(3);
rb = -[raz ray rax];
cMat = eye(3);
cz = cMat(1, :);
cy = cMat(2, :);
cx = cMat(3, :);

t40_1zyx = 5/8 * (63*raz^4*rb - 42*raz^2*rb + 28*raz^3*cz + 3*rb - 12*raz*cz);
t41c_1zyx = 1/4 * sqrt(10) * (63*rax*raz^3*rb - 21*rax*raz*rb + 21*rax*raz^2*cz + 7*raz^3*cx - 3*rax*cz - 3*raz*cx);
t41s_1zyx = 1/4 * sqrt(10) * (63*ray*raz^3*rb - 21*ray*raz*rb + 21*ray*raz^2*cz + 7*raz^3*cy - 3*ray*cz - 3*raz*cy);
t42c_1zyx = 1/4 * sqrt(5) * (63*rax^2*raz^2*rb - 63*ray^2*raz^2*rb - 7*rax^2*rb + 7*ray^2*rb + 14*rax^2*raz*cz + ...
    14*rax*raz^2*cx - 14*ray^2*raz*cz - 14*ray*raz^2*cy - 2*rax*cx + 2*ray*cy);
t42s_1zyx = 1/2 * sqrt(5) * (63*rax*ray*raz^2*rb - 7*rax*ray*rb + 14*rax*ray*raz*cz + 7*rax*raz^2*cy + 7*ray*raz^2*cx - ...
    rax*cy - ray*cx);
t43c_1zyx = 1/4 * sqrt(70) * (9*rax^3*raz*rb - 27*rax*ray^2*raz*rb + rax^3*cz + 3*rax^2*raz*cx - 3*rax*ray^2*cz - ...
    6*rax*ray*raz*cy - 3*ray^2*raz*cx);
t43s_1zyx = 1/4 * sqrt(70) * (27*rax^2*ray*raz*rb - 9*ray^3*raz*rb + 3*rax^2*ray*cz + 3*rax^2*raz*cy + 6*rax*ray*raz*cx - ...
    ray^3*cz - 3*ray^2*raz*cy);
t44c_1zyx = 1/8 * sqrt(35) * (9*rax^4*rb - 54*rax^2*ray^2*rb + 9*ray^4*rb + 4*rax^3*cx - 12*rax^2*ray*cy - 12*rax*ray^2*cx + ...
    4*ray^3*cy);
t44s_1zyx = 1/2 * sqrt(35) * (9*rax^3*ray*rb - 9*rax*ray^3*rb + rax^3*cy + 3*rax^2*ray*cx - 3*rax*ray^2*cy - ray^3*cx);

res = [t40_1zyx; t41s_1zyx; t41c_1zyx; t42s_1zyx; t42c_1zyx; t43s_1zyx; t43c_1zyx; t44s_1zyx; t44c_1zyx];
end