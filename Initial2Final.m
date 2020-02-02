ode0=0.723;
om0=0.277;
zi=19;
gi=0.99994101;
dgi=-0.00017790;
z=0.0;
g=0.764605383282211;
dg=-0.389798698868592;
Hi=0.702*100*sqrt(om0*(1+zi)^3+ode0);
H=0.702*100*sqrt(om0*(1+z)^3+ode0);
D=g/(1+z);
Ddot=H*(D+dg/(1+z));
Di=gi/(1+zi);
Didot=Hi*(Di+dgi/(1+zi));
R=52.5;
V=-2000:2000;
dr=5;
prob=zeros(2,length(V));
prob(1,:)=V;
for i=1:length(V)
    rstar=abs(R-D*V(i)/Ddot);
    r=(2*floor(rstar/5)+1)*2.5;
    pr=0;
    while r<=802.5
        x=(pinit(:,1)==r);
        vpari=Didot/D*(R/r*(R-D*V(i)/Ddot)-r);
        vperi=Didot/D*R*sqrt(1-(R-D*V(i)/Ddot)^2/r^2);
        pr=pr+exp(vpari^2*pinit(x,2)+vpari*pinit(x,3)+pinit(x,4)+vperi^2*pinit(x,5)+vperi*pinit(x,6)+pinit(x,7))*r*dr;
        r=r+dr;
    end
    prob(2,i)=pr;
end
    