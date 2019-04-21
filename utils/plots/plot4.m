function plot4 (x,z,v,fileName);

%% plot inputs
fig9 = figure(9);

subaxis(2,2,1,'Spacing',0.03,'Margin',0.04)

imagesc(x/1000,z/1000,v.Base,[min(v.Base(:)) max(v.Base(:))]);title('Baseline velocity (km/s)');axis equal tight
colorbar
ylabel('z(km)');
xlabel('x(km)');
set(gca,'FontSize',20)

hold on 
%rectangle('Position',[leftSalt topSalt rightSalt-leftSalt bottomSalt-topSalt])


subaxis(2,2,2,'Spacing',0.03,'Margin',0.04);
imagesc(x/1000,z/1000,v.Init,[min(v.Base(:)) max(v.Base(:))]);title('Initial velocity (km/s)');axis equal tight
colorbar
ylabel('z(km)');
xlabel('x(km)');
%caption('a');
set(gca,'FontSize',20)

subaxis(2,2,3,'Spacing',0.03,'Margin',0.04);
imagesc(x/1000,z/1000,v.Mon,[min(v.Base(:)) max(v.Base(:))]);title('Monitor velocity (km/s)');axis equal tight
colorbar
ylabel('z(km)');
xlabel('x(km)');
set(gca,'FontSize',20)

subaxis(2,2,4,'Spacing',0.03,'Margin',0.04);
imagesc(x/1000,z/1000,v.Mon-v.Base,[-1.2*max(max(abs(v.Mon-v.Base))) 1.2*max(max(abs(v.Mon-v.Base)))]);
title('Difference (km/s)');axis equal tight
ylabel('z(km)');
xlabel('x(km)');
colorbar
hold on 
     %rectangle('Position',[leftSalt topSalt rightSalt-leftSalt bottomSalt-topSalt])
%saveas(fig9, 'latex/Fig/model', 'epsc');

fig9.PaperPosition = [0 0 24 7];
set(gca,'FontSize',20)
print(fig9, fileName, '-depsc');

end