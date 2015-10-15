function res = Block20(~, vecZYX)
raz = vecZYX(1);
ray = vecZYX(2);
rax = vecZYX(3);

t20_00 = 1/2 * (3*raz^2 - 1);
t21c_00 = sqrt(3) * rax * raz;
t21s_00 = sqrt(3) * ray * raz;
t22c_00 = 1/2 * sqrt(3) * (rax^2 - ray^2);
t22s_00 = sqrt(3) * rax * ray;

res = [t20_00; t21s_00; t21c_00; t22s_00; t22c_00];
end