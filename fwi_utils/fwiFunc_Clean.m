%% setup
%close all
%clear all



function  FWI = fwiFunc_Clean(i, m0, D, model, opts)


% % read model, dx = 20 or 50
% dx = 50;
% v  = dlmread(['marm_' num2str(dx) '.dat']);
% 
% %v0 = 0.999*v;
% 
% % initial model
% v0 = @(zz,xx)v(1)+.7e-3*max(zz-350,0);
% 
% v1 = @(zz,xx)2*0.25*0.25*v(1)*((1-sign(zz-1500)).*(1+sign(zz-1000)).*(1-sign(xx-6000)).*(1+sign(xx-5000)));
% 
% % set frequency, do not set larger than min(1e3*v(:))/(7.5*dx) or smaller
% % than 0.5
% 
% f  = 2;%min(1e3*v(:))/(7.5*dx)*0.8
% 
% 
% 
% % receivers, xr = .1 - 10km, with 2*dx spacing, zr = 2*dx
% xr = 100:2*dx:10000;
% zr = 2*dx*ones(1,length(xr));
% 
% % sources, xr = .1 - 10km, with 4*dx spacing, zr = 2*dx
% xs = 100:4*dx:10000;
% zs = 2*dx*ones(1,length(xs));
% 
% regularization parameter
% 
% 
% %% observed data
% % grid
% n  = size(v);
% h  = dx*[1 1];
% z  = [0:n(1)-1]*h(1);
% x  = [0:n(2)-1]*h(2);
% [zz,xx] = ndgrid(z,x);
% v0 = v0(zz,xx);
% 
% v = v + v1(zz,xx);
% 
% % parameters
% model.f = f;
% model.n = n;
% model.h = h;
% model.zr = zr;
% model.xr = xr;
% model.zs = zs;
% model.xs = xs;

% model
%m = m0;
%m0 = 1./v0(:).^2;

% data
%D = F(m,model);

%% inversion

%initial model

%m0 = vec(1./v0(zz,xx).^2);
%m0 = mk;

%m = vec(1./v1(zz,xx).^2);

% misfit
fh = @(m)misfit(m0,m,D,opts,model);

% Simple BB iteration
%[mk,hist] = BBiter(fh,m0,1e-3,100);

% %projection - constraints
% l = 0.7*min(m)*ones(size(m,1),1);
% u = 1.5*max(m)*ones(size(m,1),1);


l = 0*0.5*m0;
u = 10*1.5*m0;

%%box constraints
vMin = 1.1;
vMax = 5.5;
%mWater = (1/1.5)^2;
u(:) = min(vMin^-2,u(:));
l(:) = max(vMax^-2,l(:)); 
% % grid
% n=model.n;
% h  = model.dx*[1 1];
% z  = [0:n(1)-1]*h(1);
% x  = [0:n(2)-1]*h(2);
% [zz,xx] = ndgrid(z,x);
% for ilu = 1:size(m0,1)
%     if zz(ilu)<100
%         l(ilu)=mWater*(1-0.0001);
%         u(ilu)=mWater*(1+0.0001);
%     end
%     if zz(ilu)<1200
%         l(ilu)=m0(ilu)*(1-0.0001);
%         u(ilu)=m0(ilu)*(1+0.0001);
%     end
% end
l(:)=l(:)-m0;
u(:)=u(:)-m0;

l = l.*opts.R.dvMask(:);
u = u.*opts.R.dvMask(:);

opts.x0 = 0*m0;
[mk, farsh, info] = lbfgsb(fh, l, u, opts );

mk = mk + m0;


%hist = 

% Matlab constrained iterations
% options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
% [mk,hist] = fminunc(fh,m0,options);

vk = reshape(real(1./sqrt(mk)),model.n);



FWI.opts = opts;

FWI.info = info;

FWI.final = vk;

end