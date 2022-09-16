function OCDMphi
%%
%% Last modified by Cristel Chandre (September 14, 2022)
%% Comments? cristel.chandre@cnrs.fr 
%%
close all

F0 = 0.015;
Delta_alpha = 10;
r = 4;
mu = 32548.53;
beta = 3e-10;
ti_ps = 5;
tf_ps = 20;

nx = 512;
ny = 512;

h = 1e-6;
ti = ti_ps/2.418884254e-5;
tf = tf_ps/2.418884254e-5;

p_phi = linspace(-10,100,ny);
phi = linspace(-pi,pi,nx);
[Phi,P_phi] = meshgrid(phi,p_phi);

alpha2 = F0^2/4*Delta_alpha;
mu_eff = mu*r^2;

Phi = Phi(:);
P_phi = P_phi(:);
tspan = [ti ti+1 tf];
Y0 = [Phi P_phi];
half = ceil(numel(Y0)/2);
options = odeset('RelTol',h,'AbsTol',h);
[~,yf] = ode45(@(t,y) [y(half+1:end)/mu_eff-beta*t; -alpha2*sin(2*y(1:half))],...
    tspan,Y0,options);
Pf = reshape(yf(end,half+1:end),ny,nx);
pcolor(phi,p_phi,Pf)
shading flat
hold on
plot(-asin(mu_eff*beta/alpha2)/2,mu_eff*beta*ti,'ro','MarkerSize',8,'LineWidth',3)
%plot(pi/2+asin(mu_eff*beta/alpha2)/2,mu_eff*beta*ti,'rx','MarkerSize',8,'LineWidth',3)
%plot(-asin(mu_eff*beta/alpha2)/2+pi,mu_eff*beta*ti,'ro','MarkerSize',8,'LineWidth',3)
plot(pi/2+asin(mu_eff*beta/alpha2)/2-pi,mu_eff*beta*ti,'rx','MarkerSize',8,'LineWidth',3)
colorbar
set(gca,'box','on','FontSize',20,'LineWidth',2)
xlabel(['$\phi$(' num2str(ti_ps) 'ps)'],'interpreter','latex','FontSize',26)
ylabel(['$p_\phi$(' num2str(ti_ps) 'ps)'],'interpreter','latex','FontSize',26)
ylabel(colorbar,['$p_\phi$(' num2str(tf_ps) 'ps)'],'interpreter','latex','FontSize',26)
%
% Copyright (c) 2022 Cristel Chandre.
% All rights reserved.
%
% Redistribution and use in source and binary forms are permitted provided 
% that the above copyright notice and this paragraph are duplicated in all 
% such forms and that any documentation, advertising materials, and other 
% materials related to such distribution and use acknowledge that the
% software was developed by the CNRS. The name of the CNRS may not be used 
% to endorse or promote products derived from this software without 
% specific prior written permission.

% THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED 
% WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF 
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.