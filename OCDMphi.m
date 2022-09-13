function OCDMphi
close all

F0 = 0.015;
Delta_alpha = 10;
r = 4;
mu = 32548.53;
beta = 3e-10;
tf = 20 ;

n = 512;

h = 1e-6;
tf = tf/2.418884254e-5;

p_phi = linspace(-50,50,n);
phi = linspace(-pi,pi,n);
[Phi,P_phi] = meshgrid(phi,p_phi);

alpha2 = F0^2/4*Delta_alpha;
mu_eff = mu*r^2;

Phi = Phi(:);
P_phi = P_phi(:);
tspan = [0 tf/2 tf];
Y0 = [Phi P_phi];
half = ceil(numel(Y0)/2);
options = odeset('RelTol',h,'AbsTol',h);
[~,yf] = ode45(@(t,y) [y(half+1:end)/mu_eff-beta*t; alpha2*sin(2*y(1:half))],...
    tspan,Y0,options);
Pf = reshape(yf(end,half+1:end),n,n);

pcolor(Pf)
shading flat
colorbar
