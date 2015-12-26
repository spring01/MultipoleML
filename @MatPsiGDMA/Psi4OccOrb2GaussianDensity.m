% Convert a Psi4 style occupied molecular orbital matrix to 
% a Gaussian style density matrix
function density = Psi4OccOrb2GaussianDensity(obj, psi4_occOrb)
nshell = length(obj.shellNfuncs);
shellNfuncs_ = obj.shellNfuncs;
shell2startFunc = cumsum([1 shellNfuncs_]);
shell2startFunc = shell2startFunc(1:end-1);

if(obj.psi4_sphericalAM) % spherical
    for i = 1:nshell
        if(shellNfuncs_(i) == 3) % 3p; need to change [z x y] -> [x y z]
            psi4_occOrb((1:3)+shell2startFunc(i)-1, :) = ...
                psi4_occOrb([2 3 1]+shell2startFunc(i)-1, :);
        end
    end
else % cartesian
    for i = 1:nshell
        if(shellNfuncs_(i) == 6) % 6d; need to change [xx xy xz yy yz zz] -> [xx yy zz xy xz yz]
            psi4_occOrb((1:6)+shell2startFunc(i)-1, :) = ...
                psi4_occOrb([1 4 6 2 3 5]+shell2startFunc(i)-1, :);
            psi4_occOrb((4:6)+shell2startFunc(i)-1, :) = ...
                psi4_occOrb((4:6)+shell2startFunc(i)-1, :) ./ sqrt(3);
        end
        if(shellNfuncs_(i) == 10) % 10f; need to change
            % [xxx xxy xxz xyy xyz xzz yyyyyz yzz zzz]
            % -> [xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz]
            psi4_occOrb((1:10)+shell2startFunc(i)-1, :) = ...
                psi4_occOrb([1 7 10 4 2 3 6 9 8 5]+shell2startFunc(i)-1, :);
            psi4_occOrb((1:10)+shell2startFunc(i)-1, :) = ...
                psi4_occOrb((1:10)+shell2startFunc(i)-1, :) ...
                ./ repmat([ones(1,3) sqrt(5)*ones(1,6) sqrt(15)]',1,size(psi4_occOrb,2));
        end
        if(shellNfuncs_(i) == 15) % 15g; reverse order
            psi4_occOrb((1:15)+shell2startFunc(i)-1, :) = ...
                psi4_occOrb((15:-1:1)+shell2startFunc(i)-1, :);
            psi4_occOrb((1:15)+shell2startFunc(i)-1, :) = ...
                psi4_occOrb((1:15)+shell2startFunc(i)-1, :) ...
                ./ repmat( ...
                [1 sqrt(7) sqrt(35/3) sqrt(7) 1, ...
                sqrt(7) sqrt(35) sqrt(35) sqrt(7) sqrt(35/3), ...
                sqrt(35) sqrt(35/3) sqrt(7) sqrt(7) 1]',1,size(psi4_occOrb,2));
        end
    end
end
density = 2.*psi4_occOrb*psi4_occOrb';
end