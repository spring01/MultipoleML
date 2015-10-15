function res = Block22(~, vecZYX)
raz = vecZYX(1);
ray = vecZYX(2);
rax = vecZYX(3);
rbz = -vecZYX(1);
rby = -vecZYX(2);
rbx = -vecZYX(3);
cxx = 1; cxy = 0; cxz = 0;
cyx = 0; cyy = 1; cyz = 0;
czx = 0; czy = 0; czz = 1;

t20_20 = 3/4 * (35*raz^2*rbz^2 - 5*raz^2 - 5*rbz^2 + 20*raz*rbz*czz + 2*czz^2 + 1);
t20_21c = 1/2 * sqrt(3) * (35*raz^2*rbx*rbz - 5*rbx*rbz + 10*raz*rbx*czz + 10*raz*rbz*czx + 2*czx*czz);
t20_21s = 1/2 * sqrt(3) * (35*raz^2*rby*rbz - 5*rby*rbz + 10*raz*rby*czz + 10*raz*rbz*czy + 2*czy*czz);
t20_22c = 1/4 * sqrt(3) * (35*raz^2*rbx^2 - 35*raz^2*rby^2 - 5*rbx^2 + 5*rby^2 + 20*raz*rbx*czx - 20*raz*rby*czy + 2*czx^2 - 2*czy^2);
t20_22s = 1/2 * sqrt(3) * (35*raz^2*rbx*rby - 5*rbx*rby + 10*raz*rbx*czy + 10*raz*rby*czx + 2*czx*czy);
t21c_21c = (35*rax*raz*rbx*rbz + 5*rax*rbx*czz + 5*rax*rbz*czx + 5*raz*rbx*cxz + 5*raz*rbz*cxx + cxx*czz + cxz*czx);
t21c_21s = (35*rax*raz*rby*rbz + 5*rax*rby*czz + 5*rax*rbz*czy + 5*raz*rby*cxz + 5*raz*rbz*cxy + cxy*czz + cxz*czy);
t21c_22c = 1/2 * (35*rax*raz*rbx^2 - 35*rax*raz*rby^2 + 10*rax*rbx*czx - 10*rax*rby*czy + 10*raz*rbx*cxx - 10*raz*rby*cxy + ...
    2*cxx*czx - 2*cxy*czy);
t21c_22s = (35*rax*raz*rbx*rby + 5*rax*rbx*czy + 5*rax*rby*czx + 5*raz*rbx*cxy + 5*raz*rby*cxx + cxx*czy + cxy*czx);
t21s_21s = (35*ray*raz*rby*rbz + 5*ray*rby*czz + 5*ray*rbz*czy + 5*raz*rby*cyz + 5*raz*rbz*cyy + cyy*czz + cyz*czy);
t21s_22c = 1/2 * (35*ray*raz*rbx^2 - 35*ray*raz*rby^2 + 10*ray*rbx*czx - 10*ray*rby*czy + 10*raz*rbx*cyx - 10*raz*rby*cyy + ...
    2*cyx*czx - 2*cyy*czy);
t21s_22s = (35*ray*raz*rbx*rby + 5*ray*rbx*czy + 5*ray*rby*czx + 5*raz*rbx*cyy + 5*raz*rby*cyx + cyx*czy + cyy*czx);
t22c_22c = 1/4 * (35*rax^2*rbx^2 - 35*rax^2*rby^2 - 35*ray^2*rbx^2 + 35*ray^2*rby^2 + 20*rax*rbx*cxx - 20*rax*rby*cxy - ...
    20*ray*rbx*cyx + 20*ray*rby*cyy + 2*cxx^2 - 2*cxy^2 - 2*cyx^2 + 2*cyy^2);
t22c_22s = 1/2 * (35*rax^2*rbx*rby - 35*ray^2*rbx*rby + 10*rax*rbx*cxy + 10*rax*rby*cxx - 10*ray*rbx*cyy - 10*ray*rby*cyx + ...
    2*cxx*cxy - 2*cyx*cyy);
t22s_22s = (35*rax*ray*rbx*rby + 5*rax*rbx*cyy + 5*rax*rby*cyx + 5*ray*rbx*cxy + 5*ray*rby*cxx + cxx*cyy + cxy*cyx);

res = [ ...
    t20_20  t20_21s  t20_21c  t20_22s  t20_22c; ...
    t20_21s t21s_21s t21c_21s t21s_22s t21s_22c; ...
    t20_21c t21c_21s t21c_21c t21c_22s t21c_22c; ...
    t20_22s t21s_22s t21c_22s t22s_22s t22c_22s; ...
    t20_22c t21s_22c t21c_22c t22c_22s t22c_22c];
end