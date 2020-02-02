ode0=0.723;
om0=0.277;
zi=19;
gi= 0.99994101;
dgi=-0.00017790;
om=om0*(1+zi)^3;
Hi=0.702*100*sqrt(om+ode0);
Di=gi/(1+zi);
Didot=Hi*(Di+dgi/(1+zi));
r=37.5;
corr=0;
Qpar=0;
Qper=0;
beta=0;
sig_ms=0;
sig_mp=0;
dsig_ms=0;
dsig_mp=0;
sigpar=0;
sigper=0;
i=1;
p_k=p_k19;
M_p=5e13;
M_s=5e13;
dM_s=1e12;
dM_p=1e12;
a=(vrad(:,1)==r & vrad(:,2)>=-2200 & vrad(:,2)<=2200);
p0(:,1)=vrad(a,2);
q0(:,1)=p0(:,1);
q0(:,2)=zeros(length(q0(:,1)),1);
te=vrad(a,:);
vcenter=te(vrad(a,5)==max(vrad(a,5)),2);
plot(p0(:,1),log(vrad(a,5)/sum(vrad(a,5)*100)),'r');
hold on;
rho_m=2.762e7*Hi^2*om;

while M_p<=5e14
    while M_s<=5e14
        r_p=(3*M_p/(4*pi*rho_m))^1/3;
        r_s=(3*M_s/(4*pi*rho_m))^1/3;
        while i<length(p_k)
            k=p_k(i,1);
            W_p=3*(sin(k*r_p)/(k*r_p)^3-cos(k*r_p)/(k*r_p)^2);
            W_s=3*(sin(k*r_s)/(k*r_s)^3-cos(k*r_s)/(k*r_s)^2);
            sig_ms=sig_ms+(1/2*pi^2)*k^3*abs(W_s)^2*p_k(i,2)*(p_k(i+1,1)-p_k(i,1));
            sig_mp=sig_mp+(1/2*pi^2)*k^3*abs(W_p)^2*p_k(i,2)*(p_k(i+1,1)-p_k(i,1));
            dsig_ms=dsig_ms+2*(1/2*pi^2)*k^3*abs(W_s)*p_k(i,2)*(p_k(i+1,1)-p_k(i,1))*(-9*(sin(k*r_s)/(k*r_s)^2-cos(k*r_s)/(k*r_s))+3*sin(k*r_s)/(k*r_s)^2);
            dsig_mp=dsig_mp+2*(1/2*pi^2)*k^3*abs(W_p)*p_k(i,2)*(p_k(i+1,1)-p_k(i,1))*(-9*(sin(k*r_p)/(k*r_p)^2-cos(k*r_p)/(k*r_p))+3*sin(k*r_p)/(k*r_p)^2);
            corr=corr+4*pi*p_k(i,2)*(k^2/r)*sin(k*r)*(p_k(i+1,1)-p_k(i,1));
            Qpar=Qpar+k*p_k(i,2)*(p_k(i+1,1)-p_k(i,1))*((W_p^2+W_s^2)/2-3*W_s*W_p*sin(k*r)/(k*r)+6*W_s*W_p*(sin(k*r)/(k*r)^3-cos(k*r)/(k*r)^2));
            Qper=Qper+k*p_k(i,2)*(p_k(i+1,1)-p_k(i,1))*((W_p^2+W_s^2)/2-3*W_s*W_p*(sin(k*r)/(k*r)^3-cos(k*r)/(k*r)^2));
            beta=beta+k*p_k(i,2)*(p_k(i+1,1)-p_k(i,1))*k*W_s*W_p*(sin(k*r)/(k*r)^3-cos(k*r)/(k*r)^2);
            i=i+1;
        end
        Qpar=(Qpar*(Didot/Di)^2)/(3*pi^2);
        Qper=(Qper*(Didot/Di)^2)/(3*pi^2);
        nu_s=1.673/sig_ms;
        nu_p=1.673/sig_mp;
        b_s=(0.75*nu_s^2-1)/1.673+0.6/(1.673*(1+(0.75*nu_s^2)^0.3));
        b_p=(0.75*nu_p^2-1)/1.673+0.6/(1.673*(1+(0.75*nu_p^2)^0.3));
        ns=-1.673/(4*pi*r_s^2*nu_s*M_s)*2*sqrt(nu_s^2/(2*pi))*exp(-nu_s^2/2)*dsig_ms;
        np=-1.673/(4*pi*r_p^2*nu_p*M_p)*2*sqrt(nu_p^2/(2*pi))*exp(-nu_p^2/2)*dsig_mp;
        corr=corr/(8*pi^3);
        corrhh=b_s*b_p*corr;
        sigpar=sigpar+ns*np*(1+corrhh)*Qpar;
        sigper=sigper+ns*np*(1+corrhh)*Qper;
        beta=(Didot/Di)*beta/(2*pi^2);
        beta1=-2*beta/sqrt(sigpar);
        beta2=beta^2/(sigpar);
        x=q0(:,1)./sqrt((sigpar));
        h100=x;
        h200=x.^2-1;
        p0(:,2)=exp(-(p0(:,1)-vcenter).^2/(2*sigpar))/sqrt(2*pi*sigpar);
        q0(:,2)=q0(:,2)+p0(:,2).*(1+corrhh+h100*beta1+h200*beta2)./(1+corrhh);
        M_s=M_s+dM_s;
    end
    M_p=M_p+dM_p;
end

plot(q0(:,1),log(q0(:,2)/sum(q0(:,2)*100)));
        