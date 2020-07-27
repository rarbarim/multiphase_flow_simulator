function [UX,UY] = incompressible_flow_solver2D(NX,NY,LX,LY);

%%
% NX = 5; %number of grid cells in x direction
% NY = 4; %number of grid cell in y directionn
% LX = 10; %length of the reservoir (m)
% LY = 10; %height of the reservoir (m) 
PL = 1; %input boundary condition
PR = 0; %input boundary condition
PT = 1;
PB = 0;

well = 1;
nwell = 2;
grid = [1 NX*NY];
pi = [1000 1000];
pw = [1 0];
homogenous = 1;

%% Initialization
dx = LX/NX; %grid size
x = linspace(dx/2,LX-dx/2,NX); %location of grid center
xi = linspace(0,LX,NX+1); %location of interfaces

dy = LY/NY; %grid size
y = linspace(dy/2,LY-dy/2,NY); %location of grid center
yi = linspace(0,LY,NY+1); %location of interfaces

ceff = 0.1; % effective compressibility in 1/Pa
k = 10^-13; % absolute permeability in m2
por = 0.1;
visw = 10^-3; % water viscosity
densw = 1000; % water density in kg/m3

lambda = zeros(NY,NX);
lambdaharx = zeros(NY,NX+1);
lambdahary = zeros(NY+1,NX);
TX = zeros(NY,NX+1);
TY = zeros(NY+1,NX);

switch homogenous
    case(1)
    lambda = ones(NY,NX); % input lambda
    case(0)
    lambda = rand(NY,NX);
end

[TX,TY,lambdaharx,lambdahary] = computetransmisibility2D(lambda,dx,NX,dy,NY);

%% Pressure Solver

A = zeros(NX*NY,NX*NY);
P = zeros(NX*NY,1);
q = zeros(NX*NY,1);
I = zeros(NY,NX);

for i = 1 : NX
    for j = 1:NY
            gridindex(j,i) = (j-1)*NX+i;
            [I,IN,IE,IS,IW] = findindex2D(j,i,NX);
        if i > 1 % there is left neighbor
            A(I,I)= A(I,I)+TX(j,i);
            A(I,IW)= -TX(j,i);
        end
        if i < NX % there is right neighbor
            A(I,I)= A(I,I)+TX(j,i+1);
            A(I,IE)= -TX(j,i+1);
        end
        if j > 1 % there is top neighbor
            A(I,I)= A(I,I)+TY(j,i);
            A(I,IN)= -TY(j,i);
        end
        if j < NY % there is bottom neighbor
            A(I,I)= A(I,I)+TY(j+1,i);
            A(I,IS)= -TY(j+1,i);
        end
    end
end

%% Insert BC
switch well
    case(0) % case no well
    for i = 1 : NX
        for j = 1 : NY
        [I,IN,IE,IS,IW] = findindex2D(j,i,NX);
        if i == 1; % left boundary
        A(I,I) = A(I,I)+ TX(j,i);
        q(I) = q(I)+TX(j,i)*PL;

        elseif i == NX; % right boundary
        A(I,I) = A(I,I)+ TX(j,i+1);
        q(I) = q(I)+TX(j,i+1)*PR;
        
        elseif j == 1
        A(I,I)= A(I,I)+TY(j,i);
        q(I) = q(I)+TY(j,i)*PT;
        
        elseif j == NY
        A(I,I)= A(I,I)+TY(j+1,i);
        q(I) = q(I)+TY(j+1,i)*PB;
        
        end
        end
    end
    
    case(1) % case with wells, constant PI and Pw
    [A,q] = wells2D(pi,pw,lambda,A,q,grid);
end

P = A\q;
P = P';
[UX,UY] = computevelocity2D(lambdaharx,lambdahary,dx,dy,P,NX,NY,i,j);

return
