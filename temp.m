close all
re = 3.756;
De = 0.0915;
gam = 1.0755;
mu = 32548.53;

a_s = [42.13 15.4 5.4 -5.25];
a_m = [-1599.0948228665286 1064.701691434201 -262.7958617988855...
    31.287242627165202 -1.8164825476900417 0.04141328363593082];
r_a = [5 10];
b_s = [25.29 2.87 -0.09 -0.42];
b_m = [68.28948313113857 -56.39145030911223 25.523842390023567...
    -5.33917063891859 0.5409284377558191 -0.02155141584655734];
r_b = [3 6];
al_Cl = 15.5421;
al_Cl2 = 2 * 15.5421;

F0 = 0.013;
beta = 3e-10;

N = 2048;
r = linspace(3.5, 5, N);
potential = De * (1 - exp(-gam * (r - re))).^2 - De;
para_s = a_s(1)+a_s(2)*(r-re)+a_s(3)*(r-re).^2+a_s(4)*(r-re).^3;
perp_s = b_s(1)+b_s(2)*(r-re)+b_s(3)*(r-re).^2+b_s(4)*(r-re).^3;
para_m = a_m(1)+a_m(2)*r+a_m(3)*r.^2+a_m(4)*r.^3+a_m(5)*r.^4+a_m(6)*r.^5;
perp_m = b_m(1)+b_m(2)*r+b_m(3)*r.^2+b_m(4)*r.^3+b_m(5)*r.^4+b_m(6)*r.^5;
para_l = (al_Cl2+4*al_Cl^2./r.^3)./(1-4*al_Cl.^2./r.^6);
perp_l = (al_Cl2-2*al_Cl^2./r.^3)./(1-al_Cl.^2./r.^6);
al_para = para_m;
al_para(r<=r_a(1)) = para_s(r<=r_a(1));
al_para(r>=r_a(2)) = para_l(r>=r_a(2));
al_perp = perp_m;
al_perp(r<=r_b(1)) = perp_s(r<=r_b(1));
al_perp(r>=r_b(2)) = perp_l(r>=r_b(2));
Dal = al_para-al_perp;

%figure, plot(r,potential,'k','LineWidth',3)
%figure, plot(r,al_perp,'r','LineWidth',3)
%hold on
%plot(r,al_para,'b','LineWidth',3)


eta = 4*mu*r.^2*beta./(F0^2*Dal);
figure, plot(r,eta,'-','linewidth',3)
set(gca,'box','on','FontSize',20,'LineWidth',2)
xlabel('$r$','interpreter','latex','FontSize',26)
ylabel('$\eta$','interpreter','latex','FontSize',26)


p_phi = [0, 100, 200, 300, 400];
colors = ['k', 'b', 'g', 'c', 'r'];
r = linspace(2,7,N);
potential = De * (1 - exp(-gam * (r - re))).^2 - De;
p_st = sqrt(mu*r.^3*2*De*gam.*exp(-gam*(r-re)).*(1-exp(-gam*(r-re))));
Veff_st = potential+p_st.^2./(2*mu*r.^2);
figure
for it = 1:length(p_phi)
    Veff = potential+p_phi(it)^2./(2*mu*r.^2);
    plot(r,Veff,colors(it),'LineWidth',3)
    hold on
end
set(gca,'box','on','FontSize',20,'LineWidth',2)
ylabel('$V_\mathrm{eff}$','interpreter','latex','FontSize',26)
xlabel('$r$','interpreter','latex','FontSize',26)
ylim([-0.1 0.2])
plot(r(r<=4.87),Veff_st(r<=4.87),'k','LineWidth',2)
plot(r(r>4.87),Veff_st(r>4.87),'k--','LineWidth',2)
legend([string(p_phi) '' ''])

r = linspace(re,7,512);
p_phi = sqrt(mu*r.^3*2*De*gam.*exp(-gam*(r-re)).*(1-exp(-gam*(r-re))));
figure, plot(p_phi,r,'LineWidth',3)
set(gca,'box','on','FontSize',20,'LineWidth',2)
xlabel('$p_\phi$','interpreter','latex','FontSize',26)
ylabel('$r^*$','interpreter','latex','FontSize',26)

eta = 0.03;
x = linspace(-15,15,512);
p = linspace(-30,30,512);
[X,P] = meshgrid(x,p);
H = P.^2/2+X+1/eta*cos(X).^2;
figure, pcolor(X,P,H)
colorbar
shading flat
hold on
contour(X,P,H,[10, 50, 100, 200],'-k','LineWidth',2)


eta = 0.19;
xu = pi/2+asin(eta)/2;
x = linspace(xu,2.39*xu,512);
p_p = sqrt(2)*(xu-eta^(-1)*cos(xu)^2-x+eta^(-1)*cos(x).^2).^0.5;
p_m = -sqrt(2)*(xu-eta^(-1)*cos(xu)^2-x+eta^(-1)*cos(x).^2).^0.5;
figure, plot(x,p_p,'b','linewidth',3)
hold on 
plot(x,p_m,'b','linewidth',3)

eta = 0.17;
xu = pi/2+asin(eta)/2;
x = linspace(xu,2.43*xu,512);
p_p = sqrt(2)*(xu-eta^(-1)*cos(xu)^2-x+eta^(-1)*cos(x).^2).^0.5;
p_m = -sqrt(2)*(xu-eta^(-1)*cos(xu)^2-x+eta^(-1)*cos(x).^2).^0.5;
plot(x,p_p,'r','linewidth',3)
plot(x,p_m,'r','linewidth',3)