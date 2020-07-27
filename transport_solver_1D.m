% |---*---|---*---|---*---|
% Note : number of interface = number of cells + 1

close all
clear all
format short 

%% Default input
L = 100;
N = 100;
PL = 1;
PR = 0;
t = 1e3;
rho = 1000; % water density in kg/m3
phi = 0.3; 
k = 1e-12; % absolute permeability in m2
visco_w = 1e-3;
visco_o = 10e-3;
nw = 3.5;
no = 2;
% ut = 1;
Swc = 0.2;
Sor = 0.1;
krwe = 0.7;
kroe = 0.8;
CFL = 0.5;

[U] = incompressible_flow_solver1D(N,L,PL,PR);

%% Initialization
dx = L/N; %grid size
x = linspace(dx/2,L-dx/2,N); %location of grid center
xi = linspace(0,L,N+1); %location of interfaces

Sw_dummy = linspace(Swc,1-Sor,N);

[~,~,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability(Sw_dummy,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o);
[~,~,dfwds] = fractionalflow(lambdaw,lambdao,dlambdaw,dlambdao);
speed  = max(U)/phi*max(dfwds);
dt = CFL*dx/speed;
% dt = 1;
nt = round(t/dt);

q = zeros(N,1);

%% Explicit Method
Sw = ones(N,1)*Swc;
Sw(1) = 1-Sor;
Sw(N) = Swc;

[krw,kro,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability(Sw,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o);
[fw,fo] = fractionalflow(lambdaw,lambdao,dlambdaw,dlambdao);

Sw_im = zeros(N,nt);
Sw_im(1,:)= 1 - Sor;
Sw_im(end,:) = Swc;

figure(1)
subplot(1,2,1)
hold on
for i = 1 : nt
    for j = 2 : N
    Sw(j) = Sw(j)+dt*q(j)/phi-dt/(phi*dx)*(fw(j).*U(j+1)-fw(j-1).*U(j));
    Sw_im(j,i) = Sw(j);
    end
[krw,kro,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability(Sw,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o);
[fw,fo] = fractionalflow(lambdaw,lambdao,dlambdaw,dlambdao);

plot(x,Sw)
end
xlabel('x location')
ylabel('Sw')
title('Saturation Profile in Explicit Method')
ylim([0 1])
hold off

%% Implicit method

Sw = ones(N,1)*Swc;
Sw(1) = 1-Sor;
Sw(N) = Swc;

A = zeros(N,N);
R = zeros(N,1);

figure(1)
subplot(1,2,2)
hold on;
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
%         Sw = Sw+ds;
                
        [krw,kro,lambdaw,lambdao,dlambdaw,dlambdao] = relativepermeability(Swn,Swc,Sor,no,nw,kroe,krwe,visco_w,visco_o);
        [fw,fon,dfwds] = fractionalflow(lambdaw,lambdao,dlambdaw,dlambdao);
        R = residual(q,phi,dt,dx,Sw,Swn,fw,U,N);
        
        if norm(R,2) < 1e-6;
            converged = 1;
            Sw = Swn;
        end
    end
Sw_ex(:,Ndt) = Sw;  
plot(x,Sw); 
end
xlabel('x location')
ylabel('Sw')
title('Saturation Profile in Implicit Method')
ylim([0 1])
hold off

%% Comparison Explicit vs Implicit

figure
hold on
plot(x,Sw_im(:,end),'linewidth',2)
plot(x,Sw_ex(:,end),'--','linewidth',3)
xlabel('x location')
ylabel('Sw')
title('Saturation Profile in Implicit and Explicit Method')
legend('Implicit','Explicit')
ylim([0 1])

%% Comparison Explicit vs Implicit vs Analytical
[Sw_an,x_an] = analytical(visco_w,visco_o,nw,no,max(U),Swc,Sor,krwe,kroe,phi,t,N);

figure
hold on
plot(x,Sw_im(:,end),'rs','linewidth',1.5)
plot(x,Sw_ex(:,end),'k:','linewidth',2.5)
plot(x_an,Sw_an,'b','linewidth',2)
xlabel('x location')
ylabel('Sw')
title('Numerical and Analytical Saturation Profile')
legend('Implicit','Explicit','Analytical')
ylim([0 1])
xlim([0 L])

%% Error Part
x_an = x_an(4:end);
Sw_an = Sw_an(4:end);
t = round(1/3*N);
Sw2 = interp1(x_an,Sw_an,x(2:t+1),'spline');
error_avg = mean(abs(Sw2'-Sw_ex(2:t+1,end)));
error_max = max(abs(Sw2'-Sw_ex(2:t+1,end)));

%% Sensitivity Part

[Sw_sens1_ex,Sw_sens1_im] = changeviscosity(N,q,phi,U,dx,dt,nt,Swc,Sor,no,nw,kroe,krwe,5*visco_o,visco_o);
[Sw_sens2_ex,Sw_sens2_im] = changeviscosity(N,q,phi,U,dx,dt,nt,Swc,Sor,no,nw,kroe,krwe,visco_o,visco_o);

figure
subplot(1,2,1)
hold on
plot(x,Sw_sens1_im(:,end),'b','linewidth',2)
plot(x,Sw_sens2_im(:,end),'k','linewidth',2)
plot(x,Sw_im(:,end),'r','linewidth',2)
xlabel('x location')
ylabel('Sw')
legend('\mu_w > \mu_o','\mu_w = \mu_o','\mu_w < \mu_o')
title('Sensitivity by Implicit Method')
ylim([0 1])
xlim([0 L])
hold off

subplot(1,2,2)
hold on
plot(x,Sw_sens1_ex(:,end),'b','linewidth',2)
plot(x,Sw_sens2_ex(:,end),'k','linewidth',2)
plot(x,Sw_ex(:,end),'r','linewidth',2)
xlabel('x location')
ylabel('Sw')
legend('\mu_w > \mu_o','\mu_w = \mu_o','\mu_w < \mu_o')
title('Sensitivity by Explicit Method')
ylim([0 1])
xlim([0 L])
hold off