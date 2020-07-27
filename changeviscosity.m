function [Sw_ex,Sw_im] = changeviscosity(N,q,phi,U,dx,dt,nt,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o);

%% Explicit Method
Sw = ones(N,1)*Swc;
Sw(1) = 1-Sor;
Sw(N) = Swc;

[krw,kro,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability(Sw,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o);
[fw,fo] = fractionalflow(lambdaw,lambdao,dlambdaw,dlambdao);

Sw_im = zeros(N,nt);
Sw_im(1,:)= 1 - Sor;
Sw_im(end,:) = Swc;

for i = 1 : nt
    for j = 2 : N
    Sw(j) = Sw(j)+dt*q(j)/phi-dt/(phi*dx)*(fw(j).*U(j+1)-fw(j-1).*U(j));
    Sw_im(j,i) = Sw(j);
    end
[krw,kro,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability(Sw,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o);
[fw,fo] = fractionalflow(lambdaw,lambdao,dlambdaw,dlambdao);

end

%% Implicit Method
Sw = ones(N,1)*Swc;
Sw(1) = 1-Sor;
Sw(N) = Swc;

A = zeros(N,N);
R = zeros(N,1);

for Ndt = 1:round(nt)
    % Newton Loop
    converged = 0;
    Swn = Sw;
    [krw,kro,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability(Swn,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o);
    [fw,fo,dfwds] = fractionalflow(lambdaw,lambdao,dlambdaw,dlambdao);
    R = residual(q,phi,dt,dx,Sw,Swn,fw,U,N);
    while converged == 0
        
        A(1,1) = 1;
        for i = 2:N
            A(i,i)= 1 + dfwds(i)*U(i+1)*dt/(phi*dx);
            A(i,i-1)= -dfwds(i-1)*U(i)*dt/(phi*dx);
        end

        ds = A\R;
        Swn = Swn+ds;
                
        [krw,kro,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability(Swn,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o);
        [fw,fon,dfwds] = fractionalflow(lambdaw,lambdao,dlambdaw,dlambdao);
        R = residual(q,phi,dt,dx,Sw,Swn,fw,U,N);
        
        if norm(R,2) < 1e-6;
            converged = 1;
            Sw = Swn;
        end
    end

Sw_ex(:,Ndt) = Sw;  
% plot(x,Sw); 
end

return