function [mean,disp]=pair_weighted_mean_dispersion(uw_disp,corr_hh,beta)
beta1=2*beta/sqrt(uw_disp);
beta2=beta^2/uw_disp;
mean=beta1*sqrt(uw_disp)/(1+corr_hh);
fact=(1+2*beta2/(1+corr_hh)-(beta1/(1+corr_hh))^2);
if fact<0
    fact
    uw_disp
end
disp=uw_disp*fact;
end    