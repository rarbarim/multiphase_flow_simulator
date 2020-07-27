function [sw2,x] = analytical(visco_w,visco_o,nw,no,ut,Swc,Sor,krwe,kroe,porosity,t,nsw);

% ----------------------------------------------------% Hadi Hajibeygi – All Rights Reserved – 29 April 2015
%  ANALITYCAL SOLUTION OF BUCKLEY LEVERETT EQUATION
% ----------------------------------------------------clear all; close all;
% parameters
% visco_w = 1.0e-3;
% visco_o = 10.0e-3;
% nw = 3.5;
% no = 2;
% ut = 1;
% Swc = 0.2;
% Sor = 0.1;
% krwe    = 0.7;
% kroe = 0.8;
% porosity= 0.3;

% saturation
% nsw = 101;  % number of gridpoints
dsw = (1-Sor-Swc)/(nsw-1); 
sw = Swc:dsw:(1-Sor);  % saturation range
% relperms + mobility + fracflow
krw = krwe*((sw-Swc)/(1-Swc-Sor)).^nw;
kro = kroe*((1-sw-Sor)/(1-Swc-Sor)).^no;
mob_w = krw/visco_w;
mob_o = kro/visco_o;
fw = mob_w./(mob_w+mob_o);
% fracflow derivative (numerical)
rr = 2:(nsw-1);
dfw1 = [0 (fw(rr+1)-fw(rr-1))./(sw(rr+1)-sw(rr-1)) 0];
vw    = ut/porosity*dfw1;
% jump velocity 
dfw2 = (fw - fw(1))./(sw-sw(1));
[shock_vel, shock_index] = max(dfw2);
shock_sat = sw(shock_index);
rr2 = shock_index:nsw;
% valid saturation range within BL solution
sw2 = [sw(1) sw(1) sw(rr2)];
vw2   = vw(rr2);
vw2 = [1.5*vw2(1) vw2(1) vw2];
% t  = 10;
x  = vw2*t;

%%
% plot of saturation as function of position

% figure(1); 
% subplot(2,1,2); plot(x,sw2,'.-','linewidth',2);
% xlabel('x [m]'); ylabel('S_w [-]');
% % plot of fractional flow derivative + jump velocity
% subplot(2,1,1); plot(sw,dfw1,sw,dfw2,[shock_sat shock_sat],[0 
% max(dfw1)],'r:',[Swc (1-Sor)],[shock_vel shock_vel],'r:');
% xlabel('Sw'); ylabel('derivative');
% legend('d{f_w}/d{S_w}','\Deltaf_w/\DeltaS_w');

% %%
% figure
% hold on
% plot(sw,krw)
% plot(sw,kro)
% legend('krw','kro')
% xlabel('Sw')
% ylabel('krw & kro')
% xlim([0 1])
% ylim([0 1])
% hold off

return