% |---*---|---*---|---*---|
% Note : number of interface = number of cells + 1

close all 
clear all
format short

%% Default input
LX = 100; 
LY = 100;
NX = 20;
NY = 20;
t = 300;
rho = 1000; % water density in kg/m3
phi = 0.3;
k = 1e-12; % absolute permeability in m2
visco_w = 1e-3;
visco_o = 10e-3;
nw = 3.5;
no = 2;
Swc = 0.2;
Sor = 0.1;
krwe = 0.7;
kroe = 0.8;
CFL = 0.5;

Sw_dummy = linspace(Swc,1-Sor,NX*NY);
% UX = ones(NY,NX+1);
% UY = ones(NY+1,NX);
[UX,UY] = incompressible_flow_solver2D(NX,NY,LX,LY);
% UX = 2*UX;
% UY = 2*UY;
q = zeros(NY,NX);
% q(1) = 1;
% q(end) = 1;

%% Initialization

dx = LX/NX; %grid size
x = linspace(dx/2,LX-dx/2,NX); %location of grid center
xi = linspace(0,LX,NX+1); %location of interfaces

dy = LY/NY; %grid size
y = linspace(dy/2,LY-dy/2,NY); %location of grid center
yi = linspace(0,LY,NY+1); %location of interfaces

[xx yy] = meshgrid(x,y);

[krw,kro,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability2D(Sw_dummy,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o);
[fw,fo,dfwds] = fractionalflow2D(lambdaw,lambdao,dlambdaw,dlambdao);
speed  = max(max(UX)/phi*max(dfwds));
dt = CFL*dx/speed;
nt = round(t/dt);

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
    iter =1;
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
        iter=iter+1;
    end
    
    Sw_im(:,Ndt) = Sw;
    
end

%% Graph Explicit Method

figure
subplot(1,2,1)
view(0,0);
surf(xx,yy,Sw_ex(:,:,end))
xlabel('x location')
ylabel('y location')
zlabel('Sw')
zlim([0 1])
title('2D Saturation Profile in Explicit Method')

% Graph Implicit Method
subplot(1,2,2)
view(0,90);
surf(xx,yy,vec2mat(Sw_im(:,end),NX))
xlabel('x location')
ylabel('y location')
zlabel('Sw')
zlim([0 1])
title('2D Saturation Profile in Implicit Method')

%% Sensitivity Part
[Sw_sens1_ex,Sw_sens1_im] = changeviscosity2D(NX,NY,q,phi,UX,UY,dx,dy,dt,nt,Swc,Sor,no,nw,kroe,krwe,10*visco_o,visco_o);
[Sw_sens2_ex,Sw_sens2_im] = changeviscosity2D(NX,NY,q,phi,UX,UY,dx,dy,dt,nt,Swc,Sor,no,nw,kroe,krwe,visco_o,visco_o);

%% Graph  Sensitivity Explicit Method

figure
subplot(1,3,1)
view(0,0);
surf(xx,yy,Sw_sens1_ex(:,:,end))
xlabel('x location')
ylabel('y location')
zlabel('Sw')
zlim([0 1])
title('2D Explicit Method Sensitivity \mu_w > \mu_o')

subplot(1,3,2)
view(0,90);
surf(xx,yy,Sw_sens2_ex(:,:,end))
xlabel('x location')
ylabel('y location')
zlabel('Sw')
zlim([0 1])
title('2D Explicit Method Sensitivity \mu_w = \mu_o')

subplot(1,3,3)
view(0,90);
surf(xx,yy,Sw_ex(:,:,end))
xlabel('x location')
ylabel('y location')
zlabel('Sw')
zlim([0 1])
title('2D Explicit Method Sensitivity \mu_w < \mu_o')

figure
subplot(1,3,1)
view(0,0);
surf(xx,yy,vec2mat(Sw_sens1_im(:,end),NX))
xlabel('x location')
ylabel('y location')
zlabel('Sw')
zlim([0 1])
title('2D Implicit Method Sensitivity \mu_w > \mu_o')

subplot(1,3,2)
view(0,0);
surf(xx,yy,vec2mat(Sw_sens2_im(:,end),NX))
xlabel('x location')
ylabel('y location')
zlabel('Sw')
zlim([0 1])
title('2D Implicit Method Sensitivity \mu_w > \mu_o')

subplot(1,3,3)
view(0,0);
surf(xx,yy,vec2mat(Sw_im(:,end),NX))
xlabel('x location')
ylabel('y location')
zlabel('Sw')
zlim([0 1])
title('2D Implicit Method Sensitivity \mu_w > \mu_o')
