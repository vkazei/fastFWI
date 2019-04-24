%% Example of FWI for Marmousi II benchmark model.
%
% Press "Run" to launch the script
%
% Vladimir Kazei and Oleg Ovcharenko, 2019

%% initialization
tic
% cleaning
close all; clearvars; restoredefaultpath;
spparms('bandden',0.0)

% Link core folders
addpath(genpath('../engine/'));
set(groot,'DefaultFigureColormap',rdbuMap())

% Create folder to store output figures
figFolder = 'Fig/';
mkdir_inform(figFolder);

%% MODEL and GRID
% read Marmousi II model (Baseline)
dx = 50;
v.True = dlmread('models/BP/bp_tooth_long.dat')/1e3;
v.True = imresize(v.True, 40/dx, 'bilinear');git 

% append  N_ext_up points from above for the absorbing layer
N_ext_up = 10;
v.True = [repmat(v.True(1,:),5,1);  v.True];

% grid creation
model.n  = size(v.True);
model.h  = dx*[1 1];
z  = [-N_ext_up:model.n(1)-1-N_ext_up]*model.h(1);
x  = [0:model.n(2)-1]*model.h(2);
[zz,xx] = ndgrid(z,x);

model.x = x;
model.z = z;

% build linear initial model for FWI
v_Fun_Init = @(zz,xx)v.True(1)+0.82e-3*max(zz-450,0);
v.Init = v_Fun_Init(zz,xx);

imagescc(v.True,model,'BP 2004',[figFolder,'true'])
model.caxis = caxis;
imagescc(v.Init,model,'Initial',[figFolder,'init'])

% model - converted to squared slowness
%baseline
m.Base = 1./v.True(:).^2;
%initial
m.Init = 1./v.Init(:).^2;



%% DATA ACQUISITION
% set frequency range, not larger than min(1e3*v(:))/(7.5*dx) or smaller than 0.5
fMin  = 0.5; % this is the minimum frequency used for inversion
fFactor = 1.2; % factor to the next frequncy
fMax = 2; % this is the max frequency used for inversion

% receivers
model.xr = N_ext_up*dx:dx:max(x(:))-N_ext_up*dx;
model.zr = 0*N_ext_up*dx*ones(1,length(model.xr));

% sources
model.xs = N_ext_up*dx:400:max(x(:))-N_ext_up*dx;
model.zs = 0*N_ext_up*dx*ones(1,length(model.xs));

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
opts.R.alphaTV = 0.1;
% 1 - regularize update (usual for time-lapse),
% 0 - regularize model itself (blocky model)
opts.R.dmFlag = 0;

% this restricts the updates l and u are multiplied by this
opts.R.dvMask = zz>100;

% you can modify the parameters after presetting to improve regularization
% performance e.g. "opts.R.p = 1.5;"

%% LBFGSB PARAMETERS
% lbfgs "depth" - number of gradients
opts.m = 10;
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

%
opts.tracking = true;
 opts.histFolder = 'JHist/';
    disp('Full misfit evaluations history will be saved, which affects performance')


%% MAIN LOOP OVER FREQUENCIES
it = 0;
iFWI = 0;
for iOut=1:3
freq = fMin;
for iOut=1:3
    while freq<=fMax
        it=it+1;
        % acquire noisy data
        model.maxOffset = minMaxOffset;
        model.f = freq;
        
        DClean = F(m.Base,model);
        D0 = F(m.Init,model);
        
        % Gaussian noise
        mySNR = 100;
        
        noiseStandDev = 1/sqrt(mySNR);
        noiseStandDev = noiseStandDev * sqrt(mean(mean(abs(DClean).*abs(DClean))));
        D = DClean + sqrt(1/2)*(randn(size(DClean))*noiseStandDev+1i*randn(size(DClean))*noiseStandDev);
        
        snrDB = snr(DClean, D-DClean)
        matSNR  = db2pow(snrDB)
        
        % Mute offsets beyond the max limit
        while model.maxOffset <= maxMaxOffset
            for i=1:size(model.xr,2)
                for j=1:size(model.xs,2)
                    if abs(model.xr(i)-model.xs(j))>model.maxOffset
                        D(i,j)=0;
                    end
                end
            end
            
            %  FWI ================================================
            fwiResult = fwiFunc(iFWI, m.Init, D, model, opts);
            % =====================================================
            
            figure;
            imagescc(fwiResult.final,model,[num2str(model.f) ' Hz'],[figFolder num2str(it)])
            
            % Relax regularization and increase affordable offset
            model.maxOffset = fFactor*model.maxOffset;
            opts.R.alphaTV = opts.R.alphaTV/fFactor;
            m.Init(:) = 1./fwiResult.final(:).^2;
            
        end
        m.Init = fwiResult.final(:).^-2;
        freq = freq*fFactor
    end
end
end
%%
%%%%%%%%%%%%%%%%%
% FINAL PLOT
%%%%%%%%%%%%%%%%%
figure;
imagescc(fwiResult.final,model,'Final',[figFolder 'final'])

if opts.tracking
    disp('Run create_J_movie.m later to gather all iterations in a movie')
end

toc;
