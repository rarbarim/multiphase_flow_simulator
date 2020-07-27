function [Sw_ex,Sw_im] = changeviscosity2D(NX,NY,q,phi,UX,UY,dx,dy,dt,nt,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o);

%% Explicit Method
Sw = ones(NY,NX)*Swc;

[krw,kro,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability2D(Sw,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o);
[fw,fo,dfwds] = fractionalflow2D(lambdaw,lambdao,dlambdaw,dlambdao);
A = zeros(NX*NY,NX*NY);

for Ndt = 1 : nt
    Sw(1,1) = 1 - Sor;
    
    for i = 1 : NX
        for j = 1 : NY
            if i == 1 && j > 1
                Sw(j,i) = Sw(j,i) + dt*q(j,i)/phi - dt/(phi*dy)*(fw(j,i).*UY(j+1,i)-fw(j-1,i).*UY(j,i)) - dt/(phi*dx)*(fw(j,i).*UX(j,i+1)-0.*UX(j,i));
            elseif j == 1 && i > 1
                Sw(j,i) = Sw(j,i) + dt*q(j,i)/phi - dt/(phi*dy)*(fw(j,i).*UY(j+1,i)-0.*UY(j,i)) - dt/(phi*dx)*(fw(j,i).*UX(j,i+1)-fw(j,i-1).*UX(j,i));
            elseif i > 1 && j > 1
                Sw(j,i) = Sw(j,i) + dt*q(j,i)/phi - dt/(phi*dy)*(fw(j,i).*UY(j+1,i)-fw(j-1,i).*UY(j,i)) - dt/(phi*dx)*(fw(j,i).*UX(j,i+1)-fw(j,i-1).*UX(j,i));
            end
        end
    end
    
    [krw,kro,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability2D(Sw,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o);
    [fw,fo,dfwds] = fractionalflow2D(lambdaw,lambdao,dlambdaw,dlambdao);
    
    Sw_ex(:,:,Ndt) = Sw;
end

%% Implicit Method

Sw = ones(NX*NY,1)*Swc;
Sw(1) = 1-Sor;
A = zeros(NX*NY,NX*NY);
R = zeros(NX*NY,1);

for Ndt = 1:round(nt)
    % Newton Loop
    converged = 0;
    Swn = Sw;
    [krw,kro,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability2D(Swn,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o);
    [fw,fo,dfwds] = fractionalflow2D(lambdaw,lambdao,dlambdaw,dlambdao);
    R = residual2D(q,phi,dt,dx,dy,Sw,Swn,fw,UX,UY,NX,NY);
%     iter =1;
    while converged == 0
        
        for i = 1:NX
            for j = 1:NY
                [I,IN,IE,IS,IW] = findindex2D(j,i,NX);
                gridindex(j,i) = (j-1)*NX+i;                
                if i > 1 && j > 1
                    A(I,I)=  1 + dfwds(I)*UX(j,i+1)*dt/(phi*dx) + dfwds(I)*UY(j+1,i)*dt/(phi*dy);
                    A(I,IW)= -dfwds(IW)*UX(j,i)*dt/(phi*dx);
                    A(I,IN)= -dfwds(IN)*UY(j,i)*dt/(phi*dy);
                elseif i == 1 && j ~= 1
                    A(I,I)= 1 + dfwds(I)*UX(j,i+1)*dt/(phi*dx) + dfwds(I)*UY(j+1,i)*dt/(phi*dy);
                    A(I,IN)= -dfwds(IN)*UY(j,i)*dt/(phi*dy);
                elseif j == 1 && i ~= 1
                    A(I,I)= 1 + dfwds(I)*UX(j,i+1)*dt/(phi*dx) + dfwds(I)*UY(j+1,i)*dt/(phi*dy);
                    A(I,IW)= -dfwds(IW)*UX(j,i)*dt/(phi*dx);
                end
            end
        end
        A(1,1) = 1;
        
        ds = A\R;
        Swn = Swn+ds;
        
        [krw,kro,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability2D(Swn,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o);
        [fw,fon,dfwds] = fractionalflow2D(lambdaw,lambdao,dlambdaw,dlambdao);
        R = residual2D(q,phi,dt,dx,dy,Sw,Swn,fw,UX,UY,NX,NY);
        
        if norm(R,inf) < 1e-6;
            converged = 1;
            Sw = Swn;
        end
%         iter=iter+1;
    end
    
    Sw_im(:,Ndt) = Sw;
    
end

return