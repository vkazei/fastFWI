%% run FWI for marmousi II model
restoredefaultpath

tic; close all; clearvars; cd ..; startup; cd examples
spparms('bandden',0.0);
figFolder = 'Fig';


%% MODEL and GRID
% read Marmousi II model (Baseline), dx = 10, 20 or 40
dx = 50;
v.Base = dlmread('models/marm2/marm2_10.dat');
v.Base = imresize(v.Base, 10/dx, 'bilinear');

% append  N_ext_up points from above for the absorbing layer
N_ext_up = 5;
v.Base = [repmat(v.Base(1,:),5,1);  v.Base];
  
% grid creation
model.n  = size(v.Base);
model.h  = dx*[1 1];
z  = [-N_ext_up:model.n(1)-1-N_ext_up]*model.h(1);
x  = [0:model.n(2)-1]*model.h(2);
[zz,xx] = ndgrid(z,x);

v_Fun_Init = @(zz,xx)v.Base(1)+0.9e-3*max(zz-450,0);
v.Init = v_Fun_Init(zz,xx);

plot2v(x,z,v,[figFolder,'/inputModel']);

% model - converted to squared slowness
%baseline
m.Base = 1./v.Base(:).^2;
%initial
m.Init = 1./v.Init(:).^2;



%% DATA AACQUISITION
% set frequency range, not larger than min(1e3*v(:))/(7.5*dx) or smaller than 0.5
fMin  = 1; % this is the minimum frequency used for inversion
fFactor = 1.2; % factor to the next frequncy
fMax = 5; % this is the max frequency used for inversion

% receivers, xr = .1 - 10km, with 2*dx spacing, zr = 2*dx
model.xr = 5*dx:dx:16800;
model.zr = 5*dx*ones(1,length(model.xr));

% sources, xr = .1 - 10km, with 4*dx spacing, zr = 2*dx
model.xs = 5*dx:200:16800;
model.zs = 5*dx*ones(1,length(model.xs));

% for each frequency offsets are inverted sequentially each next is *Factor of perevious 
minMaxOffset = 8000;
maxMaxOffset = 8000;

%v.Init = dlmread(['marm_' num2str(dx) '.dat']);
% dx replicates h(1)
model.dx = dx;



%% REGULARIZATION
% regFlag chooses the type of regularization to be applied 
% -2 - Tikhonov;  -1 - MS; 1 - MGS;  2 - (TV); 3 - W_p^1; 
% 0 - no regularization
flags.R = 3;
opts.R = loadRegPresets(flags.R);
% 1 - regularize update (usual for time-lapse), 
% 0 - regularize model itself (blocky model) 
opts.R.dmFlag = 0;

% this restricts the updates l and u are multiplied by this
opts.R.dvMask = zz>400;

% you can modify the parameters after presetting to improve regularization
% performance e.g. "opts.R.p = 1.5;"

%% LBFGSB PARAMETERS
% lbfgs "depth" - number of gradients
opts.m = 2;
% max number of iterations - number of search directions
opts.maxIts = 50;
%same including the number of line search steps
opts.maxTotalIts = 1000;
%output of functional through iterations
opts.printEvery = 1;
%stopping tolerance for the gradient
opts.pgtol = 10^(-9);
%relative decrease in misfit the larger the less precise
opts.factr = 1e11;

opts.histAll = 0;


%% MAIN LOOP OVER FREQUENCIES

freq = fMin;
while freq<=fMax
    %% acquire noisy data
    model.maxOffset = minMaxOffset;
    model.f = freq;
    
    DClean = F(m.Base,model);
    
    D0 = F(m.Init,model);
    
    % Gaussian noise
    mySNR = 100;
    
    noiseStandDev = 0*1/sqrt(mySNR);
    noiseStandDev = noiseStandDev * sqrt(mean(mean(abs(DClean).*abs(DClean))));
    D = DClean + sqrt(1/2)*(randn(size(DClean))*noiseStandDev+1i*randn(size(DClean))*noiseStandDev); 
      
    snrDB = snr(DClean, D-DClean)
    matSNR  = db2pow(snrDB)
        
    %% LOOP OVER OFFSETS
    while model.maxOffset <= maxMaxOffset
        %% OFFSET LIMITATION
         for i=1:size(model.xr,2)
            for j=1:size(model.xs,2)
                if abs(model.xr(i)-model.xs(j))>model.maxOffset
                    D(i,j)=0;
                end
            end
        end
     
            fwiResult = fwiFunc_Clean(i, m.Init, D, model, opts);
            figure(100)
            imagesc(fwiResult.final);
            drawnow
       
        
        model.maxOffset = fFactor*model.maxOffset;
        opts.R.alphaTV = opts.R.alphaTV/fFactor;
        %opts.R.alpha = opts.R.alpha;
        %opts.histAll = optsArr(1).histAll;
        %m.InitArr(i) = 1./FWIArr(fix((sizeReg^2+1)/2)).final(:).^2;
        
        m.Init(:) = 1./fwiResult.final(:).^2;
        
    end
    m.Init = fwiResult.final(:).^-2;
    freq = freq*fFactor
end

%%
%%%%%%%%%%%%%%%%%
% FINAL PLOT
%%%%%%%%%%%%%%%%%
fig251 = figure(251);

    vk = fwiResult.final;
    %subaxis(sizeReg,sizeReg,i,'Spacing',0.01,'Margin',0.01);
    %strbeta = sprintf('\\alpha = %.1e,\\beta = %.1e,\\epsilon = %.1e',optsArr(i).R.alpha, optsArr(i).R.betta, optsArr(i).R.epsilon);
    %strbeta = sprintf('\\alpha = %.1e',optsArr(i).R.alpha);
    %strbeta = sprintf('\\beta = %.1e,\\epsilon = %.1e', optsArr(i).R.alphaTV, optsArr(i).R.epsilon);
    imagesc(x,z,vk,[min(v.Base(:)) max(v.Base(:))]);%title(strbeta);
    
    
    axis off equal tight

fig.PaperUnits = 'inches';
fig251.PaperPosition = [0 0 24 8];
print(fig251, [figFolder,'/fwiFinal'], '-depsc');




toc;
