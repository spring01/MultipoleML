function res = Block30(~, vecZYX)
raz = vecZYX(1);
ray = vecZYX(2);
rax = vecZYX(3);

t30_00 = 1/2 * (5*raz^3 - 3*raz);
t31c_00 = 1/4 * sqrt(6) * rax * (5*raz^2 - 1);
t31s_00 = 1/4 * sqrt(6) * ray * (5*raz^2 - 1);
t32c_00 = 1/2 * sqrt(15) * raz * (rax^2 - ray^2);
t32s_00 = sqrt(15) * rax*ray*raz;
t33c_00 = 1/4 * sqrt(10) * rax * (rax^2 - 3*ray^2);
t33s_00 = 1/4 * sqrt(10) * ray * (3*rax^2 - ray^2);

res = [t30_00; t31s_00; t31c_00; t32s_00; t32c_00; t33s_00; t33c_00];
end