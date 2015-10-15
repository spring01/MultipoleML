function res = Block50(~, vecZYX)
raz = vecZYX(1);
ray = vecZYX(2);
rax = vecZYX(3);

t50_00 = 1/8 * (63*raz^5 - 70*raz^3 + 15*raz);
t51c_00 = 1/8 * sqrt(15) * (21*rax*raz^4 - 14*rax*raz^2 + rax);
t51s_00 = 1/8 * sqrt(15) * (21*ray*raz^4 - 14*ray*raz^2 + ray);
t52c_00 = 1/4 * sqrt(105) * (3*rax^2*raz^3 - 3*ray^2*raz^3 - rax^2*raz + ray^2*raz);
t52s_00 = 1/2 * sqrt(105) * (3*rax*ray*raz^3 - rax*ray*raz);
t53c_00 = 1/16 * sqrt(70) * (9*rax^3*raz^2 - 27*rax*ray^2*raz^2 - rax^3 + 3*rax*ray^2);
t53s_00 = 1/16 * sqrt(70) * (27*rax^2*ray*raz^2 - 9*ray^3*raz^2 - 3*rax^2*ray + ray^3);
t54c_00 = 3/8 * sqrt(35) * (rax^4*raz - 6*rax^2*ray^2*raz + ray^4*raz);
t54s_00 = 3/2 * sqrt(35) * (rax^3*ray*raz - rax*ray^3*raz);
t55c_00 = 3/16 * sqrt(14) * (rax^5 - 10*rax^3*ray^2 + 5*rax*ray^4);
t55s_00 = 3/16 * sqrt(14) * (5*rax^4*ray - 10*rax^2*ray^3 + ray^5);

res = [t50_00; t51s_00; t51c_00; t52s_00; t52c_00; t53s_00; t53c_00; t54s_00; t54c_00; t55s_00; t55c_00];
end