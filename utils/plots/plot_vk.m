function plot_vk(x,z,v,vk,fileName)

%% plot inputs
fig9 = figure(10);

subaxis(1,1,1,'Spacing',0.03,'Margin',0.04)

imagesc(x/1000,z/1000,vk,[min(v.Base(:)) max(v.Base(:))]);title('Inverted average velocity (km/s)');axis equal tight
colorbar
ylabel('z(km)');
xlabel('x(km)');
set(gca,'FontSize',40)
%set(gca,'xticklabel',{[]})

hold on 
%rectangle('Position',[leftSalt topSalt rightSalt-leftSalt bottomSalt-topSalt])

load('mycmap.colormap','mycmap','-mat');
colormap(mycmap);

fig9.PaperPosition = [0 0 20 7];
set(gca,'FontSize',40)

print(fig9, fileName, '-depsc');

end