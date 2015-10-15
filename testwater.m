cart = [
    8 0         0   -0.0617
    1 0   -0.7116    0.4893
    1 0    0.7116    0.4893
    8 0         0    5.9383
    1 0   -0.7116   6.4893
    1 0    0.7116   6.4893];


mp = MatPsi2(cart, '6-31+g*');
mp.SCF_RunSCF();


gdma = MatPsiGDMA(mp);
gdma.SetLimitHeavyHydrogen(5, 5);
orb = mp.SCF_OrbitalAlpha();
occOrb = orb(:, 1:mp.Molecule_NumElectrons()/2);
gdma.RunGDMA(occOrb);
gdma.RemoveCore();
for i = 1:length(gdma.limit)
    sites{i} = MultipoleExpansion.Create(gdma, i);
end

sum2 = 0;
for i = 1:3
    for j = 4:6
        sum2 = sum2 + sites{i}.InteractionWith(sites{j});
    end
end

basisSetInfo.basisSetAO = '6-31gs';
basisSetInfo.basisSetAMBO = 'sto-3g-cartesian';
basisSetInfo.path = '/home/haichen/working/MultipoleML';

[quambo, matpsi2AO] = QUAMBO.MatPsi2Interface(cart, basisSetInfo);
ao2quambo = quambo.AOtoQUAMBO;
