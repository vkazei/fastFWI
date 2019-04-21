dx=20;
v  = dlmread(['marm2_' num2str(dx) '.dat']);

load('vBase8Hz')
   % v0 = v+0.01;
     v0 = FWIArr(2).final;

fig252 = figure(1);

%load('DDWI_MS_Clean_vFinalFinalFinal');
%i=21;
%vk = FWIArr(i).final;

%subaxis(sizeReg,sizeReg,i,'Spacing',0.01,'Margin',0.01);
%strbeta = sprintf('\\beta = %.1e,\\epsilon = %.1e', optsArr(i).R.betta, optsArr(i).R.epsilon);
%imagesc(x,z,vk-v0,[-1.2*max(max(abs(v1-v))) 1.2*max(max(abs(v1-v)))]);
imagesc(x/1000,z/1000,v0,[min(v(:)) max(v(:))]);
title('Average velocity (km/s)');axis equal tight
ylabel('z(km)');
xlabel('x(km)');
colorbar
hold on 
     %rectangle('Position',[leftSalt topSalt rightSalt-leftSalt bottomSalt-topSalt])
%saveas(fig9, 'latex/Fig/model', 'epsc');

fig252.PaperPosition = [0 0 12 3.5];
set(gca,'FontSize',20)
print(fig252, 'latex/Fig/invertedBase', '-depsc');