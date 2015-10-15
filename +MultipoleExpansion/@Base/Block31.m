function res = Block31(~, vecZYX)
raz = vecZYX(1);
ray = vecZYX(2);
rax = vecZYX(3);
rb = -[raz ray rax];
cMat = eye(3);
cz = cMat(1, :);
cy = cMat(2, :);
cx = cMat(3, :);

t30_1zyx = 1/2 * (35*raz^3*rb + 15*raz^2*cz - 15*raz*rb - 3*cz);
t31c_1zyx = 1/4 * sqrt(6) * (35*rax*raz^2*rb + 5*raz^2*cx + 10*rax*raz*cz - 5*rax*rb - cx);
t31s_1zyx = 1/4 * sqrt(6) * (35*ray*raz^2*rb + 5*raz^2*cy + 10*ray*raz*cz - 5*ray*rb - cy);
t32c_1zyx = 1/2 * sqrt(15) * ( (rax^2 - ray^2)*(7*raz*rb + cz) + 2*raz*(rax*cx - ray*cy) );
t32s_1zyx = sqrt(15) * ( rax*ray*(7*raz*rb + cz) + raz*(rax*cy + ray*cx) );
t33c_1zyx = 1/4 * sqrt(10) * (7*rax^3*rb + 3*(rax^2 - ray^2)*cx - 21*rax*ray^2*rb - 6*rax*ray*cy);
t33s_1zyx = 1/4 * sqrt(10) * (-7*ray^3*rb + 3*(rax^2 - ray^2)*cy + 21*rax^2*ray*rb + 6*rax*ray*cx);

res = [t30_1zyx; t31s_1zyx; t31c_1zyx; t32s_1zyx; t32c_1zyx; t33s_1zyx; t33c_1zyx];
end