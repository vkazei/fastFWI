%% setup
%close all
%clear all



function  FWI = fwiFunc(iFWI, m0, D, model, opts)

if opts.tracking
    % Create folder to store history of misfit (f, g) evaluations
   
    
    mkdir_inform(opts.histFolder)
        
    if iFWI == 0
        J.evalnum = 0;
        J.m = m0;
        J.gFWI = zeros(model.n);
        J.fFWI = 0;
        save([opts.histFolder,'J'],'J')
    else
        load([opts.histFolder,'J'],'J')
    end     
end

% misfit
fh = @(m)misfit(m0,m,D,opts,model);

%%box constraints
l = 0.25*m0;
u = 4*m0;

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