tic
clear avg sigr sigt;
r=2:2:202;
rcorr=1e-5:300;
avg=zeros(length(r),1);
sigr=zeros(length(r),1);
sigt=zeros(length(r),1);
corr=zeros(length(r),1);
p_k19=p_k0;
p_k19(:,2)=p_k19(:,2)*(0.05/0.76)^2;
for i=1:length(rcorr)
    corr(i)=correlation(p_k19,rcorr(i));
end
for i=1:length(r)
[avg(i),sigr(i),sigt(i)]=smoothed_disp(r(i),p_k0,corr,rcorr);
% avg(i)
% sigr(i)
% sigt(i)
end

toc