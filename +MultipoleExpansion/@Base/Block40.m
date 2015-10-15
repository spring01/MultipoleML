function res = Block40(~, vecZYX)
raz = vecZYX(1);
ray = vecZYX(2);
rax = vecZYX(3);

t40_00 = 1/8 * (35*raz^4 - 30*raz^2 + 3);
t41c_00 = 1/4 * sqrt(10) * (7*rax*raz^3 - 3*rax*raz);
t41s_00 = 1/4 * sqrt(10) * (7*ray*raz^3 - 3*ray*raz);
t42c_00 = 1/4 * sqrt(5) * (7*raz^2 - 1) * (rax^2 - ray^2);
t42s_00 = 1/2 * sqrt(5) * (7*raz^2 - 1) * rax* ray;
t43c_00 = 1/4 * sqrt(70) * rax * raz * (rax^2 - 3*ray^2);
t43s_00 = 1/4 * sqrt(70) * ray * raz * (3*rax^2 - ray^2);
t44c_00 = 1/8 * sqrt(35) * (rax^4 - 6*rax^2*ray^2 + ray^4);
t44s_00 = 1/2 * sqrt(35) * rax * ray * (rax^2 - ray^2);

res = [t40_00; t41s_00; t41c_00; t42s_00; t42c_00; t43s_00; t43c_00; t44s_00; t44c_00];
end