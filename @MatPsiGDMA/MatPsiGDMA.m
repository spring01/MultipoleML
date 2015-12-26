classdef MatPsiGDMA < handle
    
    % Note for dimension check:
    % We "implicitly" use "limit", "shellNfuncs", and "shellNprims" to
    % determine the dimension of other vectors and matrices; this may result
    % in that if these 3 properties had some wrong dimension, the program
    % would rather complain about other properties, but not them.
    
    % The multipoles are stored in the order Q00, Q10, Q11c, Q11s, Q20, ...
    
    properties (SetAccess = private)
        
        %! input section; these field names are assumed by the mex driver
        % sites info
        nucleiCharges    % nuclei charges of sites; 0 for non-atom sites
        xyzSites         % cartesian coordinates of sites
        
        % shell info
        shellNfuncs     % how many basis functions each site has
        shell2atom      % on which atom each site locates
        shellNprims     % how many primitive gaussians each site has
        
        % gaussian primitive info
        primExps        % exponentials of primitive gaussians
        primCoefs       % contraction coefficients of primtive gaussians
        
        % density matrix
        density
        %! end input section
        
        
        %! output
        % GDMA output
        multipoles;
        pairs;
        
    end
    
    properties
        
        % maximum order of multipoles of sites
        limit;
        % algorithm threshold, if primExps(i) > bigexp then use old algorithm
        bigexp;
        
    end
    
    properties (Access = private)
        
        psi4_sphericalAM; % used in psi4->gaussian density format converting
        coresIncluded;
        
    end
    
    methods
        
        function obj = MatPsiGDMA(matpsi2)
            obj.shellNprims = matpsi2.BasisSet_ShellNumPrimitives();
            obj.shell2atom = matpsi2.BasisSet_ShellToCenter();
            obj.shellNfuncs = matpsi2.BasisSet_ShellNumFunctions();
            obj.nucleiCharges = matpsi2.Molecule_AtomicNumbers();
            obj.xyzSites = matpsi2.Molecule_Geometry()'; % xyzSites in Bohr
            obj.primExps = matpsi2.BasisSet_PrimExp();
            obj.primCoefs = matpsi2.BasisSet_PrimCoeffUnnorm(); % un-normalized
            
            obj.psi4_sphericalAM = matpsi2.BasisSet_IsSpherical();
            
            % default limit is 4 for non-H and 1 for H
            obj.SetLimitHeavyHydrogen(4, 1);
            
            % default algorithm threshold is 0 (always use old algorithm)
            obj.bigexp = 0.0;
            
            obj.coresIncluded = true;
        end
        
        function SetLimitHeavyHydrogen(obj, heavyLimit, hydrogenLimit)
            heavyOrNot = (obj.nucleiCharges ~= 1);
            obj.limit = zeros(size(heavyOrNot));
            obj.limit(heavyOrNot == 1) = heavyLimit;
            obj.limit(heavyOrNot == 0) = hydrogenLimit;
        end
        
        function multipoles_ = RunGDMA(obj, occOrb)
            
            % GDMA wants a Gaussian style density matrix
            obj.density = obj.Psi4OccOrb2GaussianDensity(occOrb);
            
            %!!! GDMA DRIVER MEX !!!
            [multipoles_, pos, notMoved, moved] = MatPsiGDMA.matgdma_mex(obj);
            
            % build obj.multipoles
            multipoles_ = multipoles_(2:end, 1:length(obj.limit));
            obj.multipoles = multipoles_;
            
            % build obj.pairs
            numPrims = sum(obj.shellNprims);
            for i = 1:numPrims
                for j = 1:numPrims
                    index = (i - 1) * numPrims + j;
                    obj.pairs{i, j}.xyz = reshape(pos(:, index), 1, []);
                    obj.pairs{i, j}.notMoved = notMoved(:, index);
                    obj.pairs{i, j}.moved = ...
                        reshape(moved(:, index), size(multipoles_, 1), []);
                end
            end
        end
        
        function RemoveCores(obj)
            if (obj.coresIncluded == false)
                disp('Warning: cores are not included');
                return;
            end
            obj.multipoles(1, :) = obj.multipoles(1, :) - obj.nucleiCharges;
            obj.coresIncluded = false;
        end
        
        function AddCoresBack(obj)
            if (obj.coresIncluded == true)
                disp('Warning: cores are already included');
                return;
            end
            obj.multipoles(1, :) = obj.multipoles(1, :) + obj.nucleiCharges;
            obj.coresIncluded = true;
        end
        
        function res = NthOrderMthSite(obj, nthOrder, mthSite)
            % Qn0 Qn1s Qn1c Qn2s Qn2c ...
            res = obj.multipoles((nthOrder^2 + 1):((nthOrder + 1)^2), mthSite);
        end
        
        % Convert a Psi4 style occupied molecular orbital matrix to 
        % a Gaussian style density matrix
        density = Psi4OccOrb2GaussianDensity(obj, psi4_occOrb);
        
    end
    
    methods (Static, Access = private)
        
        % mex declared as a static mathod (to wrap it under @folder)
        [multipoles_, pairsPos, mpCoeff, mpMovedCoeff] = matgdma_mex(struct);
        
    end
    
end
