function OCDMphi
%%
%% Last modified by Cristel Chandre (September 28, 2022)
%% Comments? cristel.chandre@cnrs.fr 
%%
close all

F0 = 0.03;
mu = 32548.53;
re = 3.756;
Dalpha = 17;
beta = 3e-10;

eta = 4*mu*re^2*beta/(F0^2*Dalpha);
disp(eta)

ti_ps = 5;
tf_ps = 60;

ti = ti_ps/2.418884254e-5*sqrt(beta);
tf = tf_ps/2.418884254e-5*sqrt(beta);

nx = 512;
ny = 512;
nt = 512;

h = 1e-6;

%% Phase space diagram of the angular Hamiltonian model
p = linspace(-10,10,ny);
x = linspace(-10,10,nx);
[X,P] = meshgrid(x,p);

X = X(:);
P = P(:);
tspan = [ti, tf/2, tf];
Y0 = [X P];
half = ceil(numel(Y0)/2);
options = odeset('RelTol',h,'AbsTol',h);
[~,yf] = ode45(@(t,y) [y(half+1:end); -eta^(-1)*sin(2*y(1:half))-1],...
    tspan,Y0,options);
Pf = reshape(yf(end,half+1:end),ny,nx);
imagesc(x,p,Pf)
cmap1 = gray(256).*[0.68 0.85 0.9];
cmap2 = flipud(1-gray(256)).*[1 0.8 0.8];
colormap([flip(cmap1(2:end,:));cmap2])
caxis([-40 40])
shading flat
hold on
for k = -10:10
    plot(-asin(eta)/2+k*pi,0,'ro','MarkerSize',10,'LineWidth',2)
    plot(pi/2+asin(eta)/2+k*pi,0,'rx','MarkerSize',10,'LineWidth',2)
end
colorbar
set(gca,'box','on','FontSize',20,'LineWidth',2)
xlabel(['$x$(' num2str(ti_ps) ' ps)'],'interpreter','latex','FontSize',26)
ylabel(['$p$(' num2str(ti_ps) ' ps)'],'interpreter','latex','FontSize',26)
ylabel(colorbar,['$p$(' num2str(tf_ps) ' ps)'],'interpreter','latex','FontSize',26)
xlim([min(x) max(x)])
ylim([min(p) max(p)])

%% Sample trajectories
X = -asin(eta)/2 + [0 0 0];
P = [2.8 2.7 0];
colors = ['b' 'r' 'k'];
linewidth = [3 3 3];

Y0 = [X P];
tspan = linspace(ti,tf,nt);
half = ceil(numel(Y0)/2);
options = odeset('RelTol',h,'AbsTol',h);
[~,yf] = ode45(@(t,y) [y(half+1:end); -eta^(-1)*sin(2*y(1:half))-1],...
    tspan,Y0,options);
t_au = tspan/sqrt(beta);
half = length(X);
p_phi = yf(:,half+1:end)*(mu*re^2*sqrt(beta));
phi = yf(:,1:half);
figure
for it = 1:half
    subplot(2,1,1)
    plot(t_au*2.418884254e-5,p_phi(:,it)+mu*re^2*beta*t_au.',colors(it),'linewidth',3)
    xlim([ti_ps tf_ps])
    hold on
    set(gca,'box','on','FontSize',20,'LineWidth',2)
    xlabel('$t$ (ps)','interpreter','latex','FontSize',26)
    ylabel('$p_\phi$ (a.u.)','interpreter','latex','FontSize',26)
    subplot(2,1,2)
    plot(t_au*2.418884254e-5,phi(:,it),colors(it),'linewidth',linewidth(it))
    xlim([ti_ps tf_ps])
    ylim([-3*pi pi])
    hold on
    set(gca,'box','on','FontSize',20,'LineWidth',2)
    xlabel('$t$ (ps)','interpreter','latex','FontSize',26)
    ylabel('$\phi$','interpreter','latex','FontSize',26)
end


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