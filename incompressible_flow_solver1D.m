function [U] = incompressible_flow_solver1D(N,L,PL,PR)
%% Default input

% N = 20; %number of grid cells
% L = 1; %length of the reservoir (m)
% PL = 1; %input boundary condition
% PR = 0; %input boundary condition
well = 0;
% nwell = 2;
% pi = [1000 1000];
% pw = [3 0];
% grid = [1 N];
homogenous = 1;

%% Initialization
% homogenous = 0;
dx = L/N; %grid size
x = linspace(dx/2,L-dx/2,N); %location of grid center
xi = linspace(0,L,N+1); %location of interfaces
% ceff = 0.1; % effective compressibility in 1/Pa
% k = 10^-13; % absolute permeability in m2
% por = 0.15;
% visw = 10^-3; % water viscosity
% densw = 1000; % water density in kg/m3

lambda = zeros(N,1);
lambdahar = zeros(N+1,1);
T = zeros(N+1,1);

switch homogenous
    case(1)
    lambda(1:N)= 1; % input lambda
    case(0)
    lambda(1:N) = rand(N,1);
end

[T,lambdahar] = computetransmisibility1D(lambda,dx,N);

%% Pressure Solver

A = zeros(N,N);
P = zeros(N,1);
q = zeros(N,1);

for i = 1:N
    if i > 1 % there is left neighbor
        % T(i) * (P(i)-P(i-1))
        A(i,i)= T(i);
        A(i,i-1)= -T(i);
    end
    if i < N % there is right neighbor
        A(i,i)= A(i,i)+T(i+1);
        A(i,i+1)= -T(i+1);
    end
end

%% Insert BC
switch well
    case(0) % case no well
    i = 1; % left boundary
    % T(1) *(P(1)-PL)
    A(i,i) = A(i,i)+ T(i);
    q(i) = q(i)+T(i)*PL;

    i = N; % left boundary
    % T(N+1) *(P(N)-PR)
    A(i,i) = A(i,i)+ T(i+1);
    q(i) = q(i)+T(i+1)*PR;
    
    case(1) % case with wells, constant PI and Pw
    [A,q] = wells1D(pi,pw,lambda,A,q,grid);
end

P = A\q;
P = P';
U = computevelocity1D(lambdahar,dx,P,N);

% %% Final Pressure and Velocity Plot 
% 
% figure
% hold on
% plot(x,P,'linewidth',1.2)
% xlabel('x location')
% ylabel('Pressure ')
% title('Incompressible Flow Pressure Solver')
% hold off
% 
% figure
% hold on
% plot(xi,U,'linewidth',1.2)
% xlabel('x location')
% ylabel('Velocity')
% title('Incompressible Flow Velocity Solver')
% hold off
% 
% %% Validation
% P_an = (PR-PL)/L.*x+PL;
% figure
% hold on
% plot(x,P,'linewidth',1.2)
% plot(x,P,'.','markersize',15)
% xlabel('x location')
% ylabel('Pressure ')
% title('Incompressible Flow Solver Validation')
% legend('Numerical','Analytical')
% hold off

end
