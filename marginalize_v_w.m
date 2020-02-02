clear temp2
r=22.5;
data=pvwr_dr_5;
x=data(:,1)==r & abs(data(:,2))<2000 & abs(data(:,3))<2000;
temp2=zeros(length(v),1);
temp=data(x,:);
for i=1:length(v)
    y=temp(:,2)==v(i);
    temp2(i)=sum(temp(y,4));
end
temp2=temp2/(100*sum(temp2));
figure
subplot(1,2,1)
plot(v,log10(temp2), 'LineWidth',2);
str=sprintf('R=%gMpc/h',r);
title(str);
hold on
y=data_1e14(:,1)==r & abs(data_1e14(:,2))<2000;
plot(data_1e14(y,2),log10(data_1e14(y,3)),'ro', 'MarkerSize', 6, 'LineWidth', 2);
%legend('Our model', 'N-body measurements', 'Location','SouthWest');
xlabel('v(km/s)');
ylabel('log_{10}p(v)')
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')
set(gca, 'FontSize',14, 'FontWeight', 'bold')

temp3=zeros(length(w),1);
temp=data(x,:);
for i=1:length(w)
    y=temp(:,3)==w(i);
    temp3(i)=sum(temp(y,4));
end
temp3=temp3/(100*sum(temp3));
subplot(1,2,2)
plot(w,log10(temp3),'LineWidth',2);
str=sprintf('R=%gMpc/h',r);
title(str);
hold on
y=data_1e14t(:,1)==r & abs(data_1e14t(:,2))<2000;
plot(data_1e14t(y,2),log10(data_1e14t(y,3)),'ro', 'MarkerSize', 6, 'LineWidth', 2);
%legend('Our model', 'N-body measurements', 'Location','SouthWest');
xlabel('w(km/s)');
ylabel('log_{10}p(w)')
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')
set(gca, 'FontSize',14, 'FontWeight', 'bold')