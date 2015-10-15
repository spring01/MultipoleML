function res = Block11(~, vecZYX)
colVec = [vecZYX(1); vecZYX(2); vecZYX(3)];

res = 3 .* colVec * (-colVec)' + eye(3);
end