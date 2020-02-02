function [mean_smoothed,sig_r,sig_t]=smoothed_disp(r,p_k,correl,rsample)
m_min=1e14;
m_max=5e15;
dm=(m_max-m_min)/1e3;
%dm=1e10;
    m1=m_min;
    m2=m_max;
    mean_smoothed=0;
    sig_r=0;
    sig_t=0;
    N=0;
    corr=interp1(rsample,correl,r);
    while m1<=m_max
        [b1,~]=mass_function(p_k,m1,19);
        while m2<=m_max
            %m2
            [r1vir,r2vir,x,y,beta]=dispersion_smooth(p_k,r,m1,m2);
            [~,n1]=mass_function(p_k,m1,0);
            [~,n2]=mass_function(p_k,m2,0);
            [b2,~]=mass_function(p_k,m2,19);
%             r1vir 
%             r2vir 
%             n1  
%             n2
%             b1 
%             b2
            if r>r1vir+r2vir
            corr_hh=b1*b2*corr;
            else
                corr_hh=-1;
            end
            if corr_hh~=-1
            [mean,x]=pair_weighted_mean_dispersion(x,corr_hh,beta);
            else
                mean=0;
            end
%             mean
%             x
            mean_smoothed=mean_smoothed+n1*n2*(1+corr_hh)*mean*dm^2;
            sig_r=sig_r+n1*n2*(1+corr_hh)*x*dm^2;
            sig_t=sig_t+n1*n2*(1+corr_hh)*y*dm^2;
            N=N+n1*n2*(1+corr_hh)*dm^2;
            m2=m2+dm;
        end
        m1=m1+dm;
    end
%     corr_hh
%     mean_smoothed
%     sig_r
%     sig_t
%     N
    mean_smoothed=mean_smoothed/N;
    if N==0 && sig_r==0 && sig_t==0
        sig_r=0;
        sig_t=0;
    else
    sig_r=sig_r/N;
    sig_t=sig_t/N;
    end

