% Convert a Psi4 style occupied molecular orbital matrix to 
% a Gaussian style density matrix
function density = Psi4OccOrb2GaussianDensity(obj, psi4OccOrb)
numShell = length(obj.shellNfuncs);
shellNumFunc_ = obj.shellNfuncs;
shell2StartFunc = cumsum([1 shellNumFunc_]);
shell2StartFunc = shell2StartFunc(1:end-1);

if (obj.psi4_sphericalAM) % spherical
    for i = 1:numShell
        if (shellNumFunc_(i) == 3) % 3p; need to change [z x y] -> [x y z]
            psi4OccOrb((1:3) + shell2StartFunc(i) - 1, :) = ...
                psi4OccOrb([2, 3, 1] + shell2StartFunc(i) - 1, :);
        end
    end
else % cartesian
    for i = 1:numShell
        if (shellNumFunc_(i) == 6)
            % 6d; [xx xy xz yy yz zz] -> [xx yy zz xy xz yz]
            psi4OccOrb((1:6) + shell2StartFunc(i) - 1, :) = ...
                psi4OccOrb([1, 4, 6, 2, 3, 5] + shell2StartFunc(i) - 1, :);
            psi4OccOrb((4:6) + shell2StartFunc(i) - 1, :) = ...
                psi4OccOrb((4:6) + shell2StartFunc(i) - 1, :) ./ sqrt(3);
        end
        if (shellNumFunc_(i) == 10)
            % 10f
            % [xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz] ->
            % [xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz]
            psi4OccOrb((1:10) + shell2StartFunc(i) - 1, :) = ...
                psi4OccOrb([1 7 10 4 2 3 6 9 8 5] + shell2StartFunc(i) - 1, :);
            psi4OccOrb((1:10) + shell2StartFunc(i) - 1, :) = ...
                psi4OccOrb((1:10) + shell2StartFunc(i) - 1, :) ./ ...
                repmat([ones(1, 3), sqrt(5) * ones(1, 6), sqrt(15)]', ...
                       1, size(psi4OccOrb, 2));
        end
        if (shellNumFunc_(i) == 15) % 15g; reverse order
            psi4OccOrb((1:15) + shell2StartFunc(i) - 1, :) = ...
                psi4OccOrb((15:-1:1) + shell2StartFunc(i) - 1, :);
            psi4OccOrb((1:15) + shell2StartFunc(i) - 1, :) = ...
                psi4OccOrb((1:15) + shell2StartFunc(i) - 1, :) ...
                ./ repmat([1 sqrt(7) sqrt(35/3) sqrt(7) 1, ...
                           sqrt(7) sqrt(35) sqrt(35) sqrt(7) sqrt(35/3), ...
                           sqrt(35) sqrt(35/3) sqrt(7) sqrt(7) 1]', ...
                          1, size(psi4OccOrb, 2));
        end
    end
end
density = 2 .* (psi4OccOrb * psi4OccOrb');
end