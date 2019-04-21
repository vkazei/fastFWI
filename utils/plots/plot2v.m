function plot2v(x,z,v,fileName)

%% plot inputs
fig9 = figure(9);

subaxis(2,1,1,'Spacing',0.03,'Margin',0.04)

imagesc(x/1000,z/1000,v.Base,[min(v.Base(:)) max(v.Base(:))]);title('Target velocity (km/s)');axis equal tight
colorbar
ylabel('z(km)');
%xlabel('x(km)');
set(gca,'FontSize',40)
set(gca,'xticklabel',{[]})

hold on 
%rectangle('Position',[leftSalt topSalt rightSalt-leftSalt bottomSalt-topSalt])


subaxis(2,1,2,'Spacing',0.03,'Margin',0.04);
imagesc(x/1000,z/1000,v.Init,[min(v.Base(:)) max(v.Base(:))]);title('Initial velocity (km/s)');axis equal tight
colorbar
ylabel('z(km)');
%xlabel('x(km)');
%caption('a');
set(gca,'FontSize',40)
set(gca,'xticklabel',{[]})

% subaxis(2,2,3,'Spacing',0.03,'Margin',0.04);
% imagesc(x/1000,z/1000,v.Mon,[min(v.Base(:)) max(v.Base(:))]);title('Monitor velocity (km/s)');axis equal tight
% colorbar
% ylabel('z(km)');
% xlabel('x(km)');
% set(gca,'FontSize',20)

fig9.PaperPosition = [0 0 20 15];
set(gca,'FontSize',40)

print(fig9, fileName, '-depsc');

end