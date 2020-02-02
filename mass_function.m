function [bias,n]=mass_function(p_k,m,z)
G=4.3009e-9;
om0=0.27;
rho_cr=3*(100)^2/(8*pi*G);
rho=om0*rho_cr;
r=(3*m/(4*pi*rho))^(1/3);
delta_c=1.6865*(1+z);
p=0.3;
q=0.707;
k=p_k(:,1);
w=3./(k*r).*(sin(k*r)./(k*r).^2-cos(k*r)./(k*r));
sig_sq=k.^2.*p_k(:,2).*w.^2;
sig_sq=trapz(k,sig_sq);
sig_sq=sig_sq/(2*pi^2);
int=k.^2.*p_k(:,2).*w.*sin(k.*r)./(k.*r);
dlnu=1-1/(2*pi^2*sig_sq)*trapz(k, int);
v=delta_c^2/(sig_sq);
f_v=0.3222*sqrt(2*q*v/(pi))*(1+(q*v)^(-p))*exp(-q*v/2);
n=rho*dlnu/m*f_v;
bias=1+(q*v-1)/delta_c+2*p/(delta_c*(1+(q*v)^p));

a=0.707; b=0.35; c=0.8;
%bias=1+(sqrt(a)*a*v+sqrt(a)*b*(a*v)^(1-c)-(a*v)^c/((a*v)^c+b*(1-c)*(1-c/2)))/(sqrt(a)*delta_c);