function plot3v(x,z,v,fileName)

%% plot inputs
fig9 = figure(9);

subaxis(3,1,1,'Spacing',0.03,'Margin',0.04)

imagesc(x/1000,z/1000,v.Base,[min(v.Base(:)) max(v.Base(:))]);title('Baseline velocity (km/s)');axis equal tight
colorbar
ylabel('z(km)');
%xlabel('x(km)');
set(gca,'FontSize',40)
set(gca,'xticklabel',{[]})https://www.cs.mcgill.ca/~jpineau/ReproducibilityChecklist.pdf

hold on 
%rectangle('Position',[leftSalt topSalt rightSalt-leftSalt bottomSalt-topSalt])


subaxis(3,1,2,'Spacing',0.03,'Margin',0.04);
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

subaxis(3,1,3,'Spacing',0.03,'Margin',0.04);
imagesc(x/1000,z/1000,v.Mon-v.Base,[-1.2*max(max(abs(v.Mon-v.Base))) 1.2*max(max(abs(v.Mon-v.Base)))]);
title('Difference (km/s)');axis equal tight
ylabel('z(km)');
xlabel('x(km)');
colorbar
hold on 
     %rectangle('Position',[leftSalt topSalt rightSalt-leftSalt bottomSalt-topSalt])
%saveas(fig9, 'latex/Fig/model', 'epsc');
colormap(rdbuMap());

fig9.PaperPosition = [0 0 20 15];
set(gca,'FontSize',40)

print(fig9, fileName, '-depsc');

end