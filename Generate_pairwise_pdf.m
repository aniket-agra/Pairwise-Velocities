tic
clear prob
prob=zeros(4040,3);
for k=1:4040
    prob(k,1)=5*floor(k/101)+2.5;
    prob(k,2)=-5100+100*mod(k,101);
end
for i=1:length(red)
    x1=red(i,6); y1=red(i,7); z1=red(i,8); vx1=red(i,9); vy1=red(i,10); vz1=red(i,11);
    for j=i+1:length(red)
        x2=red(j,6); y2=red(j,7); z2=red(j,8); vx2=red(j,9); vy2=red(j,10); vz2=red(j,11);
        temp=[x2-x1;x2-x1+2400;x2-x1-2400];
        temp_min=min(abs(temp));
        temp_arr=abs(temp(:,:))==temp_min;
        delx=temp(temp_arr);
        temp=[y2-y1;y2-y1+2400;y2-y1-2400];
        temp_min=min(abs(temp));
        temp_arr=abs(temp(:,:))==temp_min;
        dely=temp(temp_arr);
        temp=[z2-z1;z2-z1+2400;z2-z1-2400];
        temp_min=min(abs(temp));
        temp_arr=abs(temp(:,:))==temp_min;
        delz=temp(temp_arr);
        sep=sqrt(delx^2+dely^2+delz^2);
        vpar=((vx2-vx1)*delx+(vy2-vy1)*dely+(vz2-vz1)*delz)/sep;
        m=floor(sep/5);
        n=floor((vpar+5050)/100)+1;
        ind=m*101+n;
        if ind<=4040
            prob(ind,3)=prob(ind,3)+1;
        end
    end
end
toc