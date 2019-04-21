%% preload models
%load('vMon50Hz.mat');
load('CDWI_All_Full.mat');
%load('velElast.mat');
%jafkl
tic;

close all
%clear all

cd .. 

startup


cd example
fFlag = 1;
% load previous result for intitial model
previousFlag = 0;



% read model, dx = 10, 20 or 40
dx = 20;

%% reading Otway model
% cd ../CO2CRCOtway2C_data_for_KAUST/Velocity_model_2D
% vp = read_binary_matrix(1220,984,['B_Sodas_Ln.vp']);
% vp1 = read_binary_matrix(1220,984,['M_Sodas_Ln.vp']);
% cd ../../example

% resizing

dv1 = v.Mon-v.Base;
v0 = v.Init;
v = v.Base;
v1 = v+dv1;




% initial model
%v0 = @(zz,xx)v(1)+1.0e-3*max(zz-150,0);

% set frequency, do not set larger than min(1e3*v(:))/(7.5*dx) or smaller
% than 0.5

    


% frequencies
fMin  = 50; %0.5<fmin(1e3*v(:))/(7.5*dx)*0.8 this is the minimum frequency used for inversion
fFactor = 1.2;
fMax = 50;


% receivers, xr = .1 - 10km, with 2*dx spacing, zr = 2*dx
xr = 450:15:1860;
zr = 2*dx*ones(1,length(xr));

% sources, xr = .1 - 10km, with 4*dx spacing, zr = 2*dx
xs = 420:15:2010;
zs = 2*dx*ones(1,length(xs));

minMaxOffset = 10000;
maxMaxOffset = 10000;

% grid
n  = size(v);
h  = dx*[1 1];
z  = [0:n(1)-1]*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

%v0 = v0(zz,xx);
%(v+v1)/2;
%v0 = dlmread(['marm_' num2str(dx) '.dat']);


% if previousFlag
%    load('vMon50Hz.mat');
%    v0 = imresize(vk, 5/dx);
% %     load('vAverageOtway10Hz.mat')
% %    v0 = v+0.01;
% %      v0 = FWIArr(1).final;
% %     v0 = imresize(v0,size(v));
% end;


% % parameters
% model.f = fMin;
% model.n = n;
% model.h = h;
% model.zr = zr;
% model.xr = xr;
% model.zs = zs;
% model.xs = xs;
% model.dx = dx;
% model.maxOffset = minMaxOffset; 

% model - squared slowness
%baseline
m = 1./v(:).^2;
%initial
m0 = 1./v0(:).^2;



%monitor
m1 = 1./v1(:).^2;

%% LBFGS parameters

% lbfgs "depth" - number of gradients
opts.m = 2;
% max number of iterations - number of search directions

opts.maxIts = 1000;

opts.maxTotalIts = 10000;

%same including the number of line search steps
%opts.maxTotalIts = 2;
opts.printEvery = 1;
%stopping tolerance for the gradient
opts.pgtol = 10^(-10);

%relative decrease in misfit the larger the less precise
opts.factr = 1e4;

%% regularization parameters;

% the power of the gradient in norm
opts.R.p = 1.1;

opts.R.alphaTV = 1;

opts.R.epsilon = 2*10^-8; % W12  

opts.R.alpha = 0;

% 1 - regularize model update, 0 - regularize model 
opts.R.dmFlag = 0;
% 1 - use norms of the update of the model (Sobolev space norms), 
% 0 - use norms of the model itself

opts.R.gFlag = 0;

OFlag = 0; % 1 - Apply variation based Oleg's flooding (Salt oriented)
VFlag = 1; % 1 - Apply multistage regularization (Salt oriented)

% regFlag chooses the type of regularization to be applied
% 0 - minimum support regularization
% 1 - MGS 
% 2 - Wp (TV)
regFlag = 2; 

switch regFlag
    case -1 %Sobolev (Tikhonov with derivative regularizaton)
        opts.R.alpha = 100;
        opts.R.gFlag = 1;
        opts.R.alphaTV = 0;
        opts.R.betta = 0;
        opts.R.epsilon = 2*10^-5;
    case 0 % minimum support
        opts.R.betta = 5*10^-6; % Min Supp
        opts.R.epsilon = 10^-5; % Min Supp
        opts.R.gFlag = 0;
        opts.R.alphaTV = 0.1;
    case 1 % MGS
        opts.R.gFlag = 1;
        opts.R.epsilon = 1*10^-6; % MGS
        opts.R.betta = 5*10^-5; % MGS, TV
        opts.R.alphaTV = 1;
    case 2 
        opts.R.betta = 0; % TV
        opts.R.gFlag = 1;
        opts.R.alphaTV = 1;
        opts.R.epsilon = 10^-18;
        opts.R.p = 1.1; % W12load
end
        


%% plot inputs
fig9 = figure(9);
%load('vMon50Hz.mat');



subaxis(2,2,1,'Spacing',0.03,'Margin',0.04)

imagesc(x/1000,z/1000,v1-v,[-1.2*max(max(abs(v1-v))) 1.2*max(max(abs(v1-v)))]);
title('Difference (km/s)');axis equal tight
colorbar
ylabel('z(km)');
xlabel('x(km)');
set(gca,'FontSize',20)
hold on 
%rectangle('Position',[leftSalt topSalt rightSalt-leftSalt bottomSalt-topSalt])


subaxis(2,2,2,'Spacing',0.03,'Margin',0.04);
imagesc(x/1000,z/1000,vk-v0,[-1.2*max(max(abs(v1-v))) 1.2*max(max(abs(v1-v)))]);
title('Inverted difference (km/s)');axis equal tight
colorbar
ylabel('z(km)');
xlabel('x(km)');
%caption('a');
set(gca,'FontSize',20)


subaxis(2,2,3,'Spacing',0.03,'Margin',0.04);
imagesc(x/1000,z/1000,vk-v0,[-1.2*max(max(abs(v1-v))) 1.2*max(max(abs(v1-v)))]);
title('Error in difference (km/s)');axis equal tight
colorbar
ylabel('z(km)');
xlabel('x(km)');
set(gca,'FontSize',20)

subaxis(2,2,4,'Spacing',0.03,'Margin',0.04);
imagesc(x/1000,z/1000,(vk-v0)./v);
title('Relative error in difference');axis equal tight
ylabel('z(km)');
xlabel('x(km)');
colorbar
hold on 
     %rectangle('Position',[leftSalt topSalt rightSalt-leftSalt bottomSalt-topSalt])
%saveas(fig9, 'latex/Fig/model', 'epsc');

fig9.PaperPosition = [0 0 20 20];
set(gca,'FontSize',20)
print(fig9, 'Fig/spec_model', '-depsc');

hold off
opts.histAll = 0;


%% plotting Spectra


% prepare frequency vectors

% scale frequency vectors with 1/lambda_min
vAnom = 3000;
lambda_min=(vAnom/fMax);

fSpace.z = lambda_min*(dx^-1)*((0:n(1)-1)/n(1)-0.5);
fSpace.x = lambda_min*(dx^-1)*((0:n(2)-1)/n(2)-0.5);


%%    


figSpec = figure;

sPlot = subaxis(1,1,1,'Spacing',0.1,'Margin',0.15);

imagesc(fSpace.x,fSpace.z,abs(fftshift(fft2(v1-v))));
title('Pertubation spectrum');axis equal tight
colorbar
ylabel('K_z (1/m)');
xlabel('K_x (1/m)');
set(gca,'FontSize',20)
a=1.5*caxis(sPlot);
caxis(a);

figSpec.PaperPosition = [0 0 6 6];
set(gca,'FontSize',20)
print(figSpec, 'Fig/pertSpectrum', '-depsc');

hold on 

imagesc(fSpace.x,fSpace.z,abs(fftshift(fft2(vk-v0))));
title('Inverted spectrum');axis equal tight
colorbar
ylabel('K_z (1/m)');
xlabel('K_x (1/m)');
set(gca,'FontSize',20)
caxis(a);

figSpec.PaperPosition = [0 0 6 6];
set(gca,'FontSize',20)
print(figSpec, 'Fig/invSpectrum', '-depsc');

%%
dv_Inv = vk - v0;

dv_True = v1 - v;

dv_True(1:end/2,:) = 0;

dv_Inv(1:end/2,:) = 0;



title('Spectral mask');axis equal tight

relLog = fftshift(log10(abs(fft2(dv_Inv-dv_True)./fft2(dv_Inv))));
maskK = relLog<-0.5;
imagesc(fSpace.x,fSpace.z,maskK);
colorbar
ylabel('K_z (1/m)');
xlabel('K_x (1/m)');
set(gca,'FontSize',20)
caxis([0 1]);

figSpec.PaperPosition = [0 0 6 6];
set(gca,'FontSize',20)
print(figSpec, 'Fig/mask', '-depsc');

%%
maskK = ifftshift(maskK);

%%
% total variation based spectrum extrapolation
% we use l-bfgs and try to keep what is under mask fixed
% we optimize the standard misfit for TV))

opts.x0=dv_Inv(:);
fh = @(m)misfitSS(dv_Inv(:),m,maskK,opts,model);
l = -v1(:);
u = v1(:);

[mk, farsh, info] = lbfgsb(fh, l, u, opts);

vTV = reshape(mk,model.n);
figure; imagesc(vTV);

%%
figure(figSpec); 
imagesc(fSpace.x,fSpace.z,abs(fftshift(fft2(vTV))));
title('Extrapolated spectrum');axis equal tight
colorbar
ylabel('K_z (1/m)');
xlabel('K_x (1/m)');
set(gca,'FontSize',20)
caxis(a);

figSpec.PaperPosition = [0 0 6 6];
set(gca,'FontSize',20)
print(figSpec, 'Fig/procSpectrum', '-depsc');



%%
% subaxis(2,2,3,'Spacing',0.03,'Margin',0.04);
% imagesc(fSpace.x,fSpace.z,log(abs(fftshift(fft2(vk-v0)))));
% 
% title('Inverted spectrum (log scale)');axis equal tight
% colorbar
% ylabel('K_z (1/m)');
% xlabel('K_x (1/m)');
% set(gca,'FontSize',20)     %rectangle('Position',[leftSalt topSalt rightSalt-leftSalt bottomSalt-topSalt])
%saveas(fig9, 'latex/Fig/model', 'epsc');

%% filtering and TV

% mask creation

%maskK = relLog<-3;
% figure;
% imagesc(maskK);
% figure;
% %maskK = medfilt2(maskK,[2 5]);
% imagesc(maskK);
% title Mask

load('CDWI_All_Full.mat','v');

figure();
vFiltered = real(ifft2(fft2(vk-v0).*maskK));
plot_dvk(x,z,v,v.Init+vFiltered,'Fig/DV_Filtered');


figure();
plot_dvk(x,z,v,v.Init+vTV,'Fig/DV_TV');

% figure();
% plot_dvk(x,z,v, v.Init-v.dv1, 'Fig/DV_True');


