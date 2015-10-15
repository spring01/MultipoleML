function multipoleExpansion = Create(matpsiGDMA, iSite)
switch(matpsiGDMA.limit(iSite))
    case(0)
        multipoleExpansion = MultipoleExpansion.MaxOrder0();
    case(1)
        multipoleExpansion = MultipoleExpansion.MaxOrder1();
    case(2)
        multipoleExpansion = MultipoleExpansion.MaxOrder2();
    case(3)
        multipoleExpansion = MultipoleExpansion.MaxOrder3();
    case(4)
        multipoleExpansion = MultipoleExpansion.MaxOrder4();
    case(5)
        multipoleExpansion = MultipoleExpansion.MaxOrder5();
    otherwise
        throw(MException('MultipoleExpansion:Create', 'Order not implemented yet.'));
end
multipoleExpansion.InitializeFromGDMA(matpsiGDMA, iSite);
end