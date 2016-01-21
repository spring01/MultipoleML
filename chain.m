molCart = ...
[6                 -5.33234853   -0.15509601    0.00000000
 1                 -4.97569411   -1.16390601    0.00000000
 1                 -4.97567569    0.34930218   -0.87365150
 1                 -6.40234853   -0.15508283    0.00000000
 6                 -4.81900631    0.57086026    1.25740497
 1                 -5.17566095    1.57967019    1.25740510
 1                 -5.17567901    0.06646161    2.13105626
 6                 -3.27900631    0.57084109    1.25740459
 1                 -2.92235253   -0.43796724    1.25936086
 1                 -2.92233341    1.07354414    0.38277663
 6                 -2.76566407    1.29923515    2.51339895
 1                 -3.12270220    0.79679017    3.38802619
 1                 -3.12195314    2.30817189    2.51121874
 6                 -1.22566421    1.29865953    2.51372134
 1                 -0.86862693    1.80108686    1.63908361
 1                 -0.86937490    0.28972291    2.51592171
 6                 -0.71232195    2.02707851    3.76970124
 1                 -1.06936295    1.52465291    4.64433845
 1                 -1.06860814    3.03601621    3.76749925
 6                  0.82767791    2.02649849    3.77002618
 1                  1.18471807    2.52892990    2.89539197
 1                  1.18396434    1.01756085    3.77222125
 6                  1.34102017    2.75490868    5.02601118
 1                  0.98398060    2.25247595    5.90064487
 1                  2.41102007    2.75450721    5.02623606
 1                  0.98473253    3.76384590    5.02381713];


mol = Molecule(molCart);

matpsi = MatPsi2(molCart, '6-31g');

matpsi.SCF_RunSCF();


gdma = MatPsiGDMA(matpsi);
gdma.SetLimitHeavyHydrogen(5, 5);
orb = matpsi.SCF_OrbitalAlpha();
occOrb = orb(:, 1:matpsi.Molecule_NumElectrons()/2);
gdma.RunGDMA(occOrb);


fock = matpsi.SCF_FockAlpha();
tei = matpsi.Integrals_AllTEIs();
dens = matpsi.SCF_DensityAlpha();

nbf = matpsi.BasisSet_NumFunctions();
func2center = matpsi.BasisSet_FuncToCenter();
func2shell = matpsi.BasisSet_FuncToShell();
shellNumPrims = matpsi.BasisSet_ShellNumPrimitives();


cutoffAtomNum = 11;
iFunc = 4;
jFunc = 1;
iShell = func2shell(iFunc);
jShell = func2shell(jFunc);
iPrimOffset = sum(shellNumPrims(1:iShell-1));
jPrimOffset = sum(shellNumPrims(1:jShell-1));

cutoffFunc = find(func2center == cutoffAtomNum);
cutoffFunc = cutoffFunc(1);
cutoffShell = func2shell(cutoffFunc);
cutoffPrim = sum(shellNumPrims(1:cutoffShell-1));

localIndices = 1:cutoffFunc-1;
densRemote = dens;
for i = localIndices
    for j = localIndices
        densRemote(i, j) = 0;
    end
end

% remoteIndices = cutoffFunc:nbf;
% densRemote = zeros(size(dens));
% for i = remoteIndices
%     for j = remoteIndices
%         densRemote(i, j) = dens(i, j);
%     end
% end

coulombRemoteRef = reshape(tei(iFunc,jFunc,:,:), 1, []) * reshape(densRemote, [], 1);

%%
iNnumPrim = shellNumPrims(func2shell(iFunc));
jNnumPrim = shellNumPrims(func2shell(jFunc));
pairFunc = cell(iNnumPrim, jNnumPrim);
for i = 1:iNnumPrim
    for j = 1:jNnumPrim
        pair = gdma.pairs{i+iPrimOffset, j+jPrimOffset};
        if norm(pair.xyz) ~= 0
            pairFunc{i, j} = MultipoleExpansion.MaxOrder3();
            pairFunc{i, j}.InitializeFromMultipoles(pair.notMoved(1:16), pair.xyz);
            pairFunc{i, j}.AddCoeffsBy([-pairFunc{i, j}.coeffs(1:3); 0; -pairFunc{i, j}.coeffs(5:end)]);
        end
    end
end


pairRemote = cell(size(gdma.pairs));
for k = 1:size(gdma.pairs, 1)
    for l = 1:k
        pairRemote{k, l} = MultipoleExpansion.MaxOrder3();
        if k > cutoffPrim || l > cutoffPrim
            pairRemote{k, l}.InitializeFromMultipoles(gdma.pairs{k, l}.notMoved(1:16), ...
                                                      gdma.pairs{k, l}.xyz);
        end
    end
end

coulombRemote = 0;
% for i = 1:iNnumPrim
%     for j = 1:jNnumPrim
%         if ~isempty(pairFunc{i, j})
%             for k = 1:size(gdma.pairs, 1)
%                 for l = 1:k
%                     if k > cutoffPrim || l > cutoffPrim
%                         coulombRemote = coulombRemote + pairFunc{i, j}.InteractionWith(pairRemote{k, l});
%                     end
%                 end
%             end
%             disp(j);
%         end
%     end
% end


pairRemoteMoved = cell(1, size(gdma.xyzSites, 2));
for p = 1:size(gdma.xyzSites, 2)
    pairRemoteMoved{p} = MultipoleExpansion.MaxOrder3();
    pairRemoteMoved{p}.InitializeFromMultipoles(zeros(16, 1), gdma.xyzSites(:, p));
end


for k = 1:size(gdma.pairs, 1)
    for l = 1:k
        pair = gdma.pairs{k, l};
        if k > cutoffPrim && l > cutoffPrim
            indMoveTo = find(sum(pair.moved.^2) ~= 0);
            for ind = indMoveTo
                pairRemoteMoved{ind}.AddCoeffsBy(pair.moved(1:16, ind));
            end
        end
    end
end

coulombRemoteMoved = 0;
for i = 1:iNnumPrim
    for j = 1:jNnumPrim
        if ~isempty(pairFunc{i, j})
            for p = 1:length(pairRemoteMoved)
                if norm(pairRemoteMoved{p}.coeffs) ~= 0
                    coulombRemoteMoved = coulombRemoteMoved + pairFunc{i, j}.InteractionWith(pairRemoteMoved{p});
                end
            end
            disp(j);
        end
    end
end

coulombRemoteRef
coulombRemote / 4 / dens(iFunc, jFunc) / ((iFunc == jFunc) + 2 * (iFunc ~= jFunc))
coulombRemoteMoved / 4 / dens(iFunc, jFunc) / ((iFunc == jFunc) + 2 * (iFunc ~= jFunc))



