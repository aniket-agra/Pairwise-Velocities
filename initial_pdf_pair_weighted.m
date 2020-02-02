r=12.5;
m_min=5e13;
m_max=1e14;
m1=m_min;
m2=m_min;
dm=(m_max-m_min)/1e3;
v=-1950:100:1950;
corr=correlation(p_k19,r);
prob=zeros(length(v),1);
tic
for i=1:length(v)
        norm=0;
        temp=0;
        while m1<=m_max
            while m2<=m_max
            [x,y,beta]=dispersion_smooth(p_k0,r,m1,m2);
            b1=interp1(m,b,m1);
            n1=interp1(m,n,m1);
            b2=interp1(m,b,m2);
            n2=interp1(m,n,m2);
            corr_hh=b1*b2*corr;
            h1=v(i)/sqrt(x);
            h2=h1^2-1;
            beta1=2*beta/sqrt(x);
            beta2=beta^2/x;
            p0=exp(-v(i)^2/(2*x));
            q0=(1+corr_hh+h1*beta1+h2*beta2)/(1+corr_hh)*p0;
            temp=temp+n1*n2*(1+corr_hh)*q0;
            norm=n+n1*n2*(1+corr_hh);
            end
        end
        temp=temp/norm;
        prob(i)=prob(i)+temp*100;
end
prob=prob/sum(prob*100);
plot(v,prob);
toc       
        
        
        
        
