function [krw,kro,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability2D(Sw,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o)

    krw = krwe.*((Sw-Swc)/(1-Swc-Sor)).^nw; 
    kro = kroe.*((1-Sw-Sor)/(1-Swc-Sor)).^no;
    
    lambdaw = krw./visco_w;
    lambdao = kro./visco_o; 
    
    dlambdaw = nw*(krwe./visco_w).*((Sw-Swc)./(1-Swc-Sor)).^(nw-1);
    dlambdao = -no*(kroe./visco_o).*((1-Sw-Sor)./(1-Swc-Sor)).^(no-1);
return