function dy =  ode_growth_linkingT3A(t,y,par,flist)
    
       
    kmet = par(2);
    uR_aa =  par(3); 
    NR = par(4);
    NP = par(5);
    A = par(6); %accuracy 
    b = par(7);
    c = par(8);
    
    Mf = flist{1};
    Vf = flist{2};
    ktrans = flist{3};
    T3convert = flist{4}; 
    kmeteff = flist{5};
    q = flist{6};  
    t3a = flist{7};
    degrade = flist{8};
%     degrade = @(x) 0.01;
    
    dy = NaN*ones(3,1);
    %1 - 8: R - P - aa 
    if T3convert(y) > 0
    dy(1) = 1/NR*uR_aa*y(1)*ktrans(y)*t3a(A,b,c) - degrade(A)*y(1);
    dy(2) = 1/NP*(1-uR_aa)*y(1)*ktrans(y)*t3a(A,b,c) - degrade(A)*y(2);
    dy(3) = y(2)*kmeteff(y)-y(1)*ktrans(y)*t3a(A,b,c); %-0.03*Mf(y);    
    end
    
end