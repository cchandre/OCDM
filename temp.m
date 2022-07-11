re = 3.756;
De = 0.0915;
gam = 1.0755;
mu = 32548.53;

p_phi = 100;

N = 512;

r = linspace(4, 5, N);
eps = De * (1 - exp(-gam * (r - re))).^2 - De;
d_eps = 2*De*gam*(1 - exp(-gam * (r - re))).*exp(-gam * (r - re));

Veff = p_phi^2./(2*mu*r.^2)+eps;

plot(r,Veff,'LineWidth',2)
plot(r,d_eps./r,'r','LineWidth',3)