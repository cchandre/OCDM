re = 3.756;
De = 0.0915;
gam = 1.0755;
mu = 32548.53;

p_phi = 100;

N = 512;

r = linspace(2, 10, N);
eps = De * (1 - exp(-gam * (r - re))).^2 - De;

Veff = p_phi^2./(2*mu*r.^2)+eps;

plot(r,Veff,'LineWidth',2)