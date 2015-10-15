function res = Block21(~, vecZYX)
raz = vecZYX(1);
ray = vecZYX(2);
rax = vecZYX(3);
rb = -[raz ray rax];
cMat = eye(3);
cz = cMat(1, :);
cy = cMat(2, :);
cx = cMat(3, :);

t20_1zyx = 1/2 * (15*raz^2*rb + 6*raz*cz - 3*rb);
t21c_1zyx = sqrt(3) * (rax*cz + cx*raz + 5*rax*raz*rb);
t21s_1zyx = sqrt(3) * (ray*cz + cy*raz + 5*ray*raz*rb);
t22c_1zyx = 1/2 * sqrt(3) * (5*(rax^2 - ray^2)*rb + 2*rax*cx - 2*ray*cy);
t22s_1zyx = sqrt(3) * (5*rax*ray*rb + rax*cy + ray*cx);

res = [t20_1zyx; t21s_1zyx; t21c_1zyx; t22s_1zyx; t22c_1zyx];
end