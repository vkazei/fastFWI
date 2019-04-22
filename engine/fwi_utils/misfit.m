function [f,g] = misfit(m0,m,D,opts,model)
% Evaluate least-squares misfit
%
%   0.5||P^TA^{-1}(m)Q - D||_{F}^2 
%   + 0.5\alpha0||m||_2^2 + 0.5\alpha1||Lm||_2^2
%   + 0.5\betta0 ||m/(m+\epsilon0)||_2^2,
%   + + 0.5\betta1 |||\nabla m|^2/(|\nabla m|^2+\epsilon1)||_1,
%
% where P, Q encode the receiver and source locations and L is the first-order FD matrix
%
% use:
%   [f,g,H] = misfit(m,D,model)
%
% input:
%   m - squared-slownes [s^2/km^2]
%   D - single-frequency data matrix
%   regt.alpha1 - regularization parameter (TV) - Tikhonov 1st derivative
%   .betta0, .epsilon0 - Minimum support
%   .betta1, .epsilon1 - MGS regularization parameters
%   .alpha0 - m^2 (Tikhonov - 0 order) 
%   .
%   model.h - gridspacing in each direction d = [d1, d2];
%   model.n - number of gridpoints in each direction n = [n1, n2]
%   model.f - frequency [Hz].
%   model.{zr,xr} - {z,x} locations of receivers [m] (must coincide with gridpoints)
%   model.{zs,xs} - {z,x} locations of sources [m] (must coincide with gridpoints)
%
%
% output:
%   f - value of misfit
%   g - gradient (vector of size size(m))
%   H - GN Hessian (function handle)

alpha1 = opts.R.alpha;
betta1 = opts.R.betta;
epsilon1 = opts.R.epsilon;
alphaTV = opts.R.alphaTV;
dmFlag = opts.R.dmFlag;
gFlag = opts.R.gFlag;
p = opts.R.p;


%epsilonTV = opts.R.epsilonTV;


%% 
%m = m0 + m;
m = m(:);
m = m0 + m;
%%  get matrices
L = getL(model.h,model.n);
A = getA(model.f,m,model.h,model.n);
P = getP(model.h,model.n,model.zr,model.xr);
Q = getQ(model.h,model.n,model.zs,model.xs);
G = @(u)getG(model.f,m,u,model.h,model.n);

%% forward banded matrix based solver
%U = A\Q;

spparms('bandden',0.0);
%spparms('spumoni',2);
%spparms('umfpack',0);
tic;
U = A\Q;
toc;
spparms('default');


%% MGS regularization functional
%m2d = reshape(m, model.n(1), model.n(2));

if dmFlag == 1
    % regularize the update of the model
    dm = m-m0;
else
    % regularize the final model
    dm = m;
end

% minimum gradient support 

if gFlag == 1
    % minimum gradient support 
    LReg = L;
else
    %minimum support
    LReg = 1;
end


%TVvec = dm.*((LReg'*LReg)*dm);

TVvec2 = (LReg*dm).^2;

fprintf('Average epsilon to gradient %.3e; ',epsilon1/max(TVvec2));

MGSfunctional = norm((TVvec2./(TVvec2 + epsilon1)).^0.5)^2;

%p=opts.R.p;

TVfunctional = norm((TVvec2 + epsilon1).^(p/4))^2;

%% compute f with regularization terms
% full data
DOBS = P'*U;

% crop from data
DOBS = DOBS.*sign(abs(D));


f = .5*norm(DOBS - D,'fro')^2 + .5*alpha1*norm(L*dm)^2 + .5*betta1*MGSfunctional + ...
    + .5*alphaTV*TVfunctional;

regPart = 1-.5*norm(DOBS - D,'fro')^2/f;
fprintf('Regularization ratio = %.3e \n',regPart);

%% adjoint solve

%tic;
spparms('bandden',0.0);
V = A'\(P*(D - DOBS));
spparms('default');
%toc;

%% compute g
g = alpha1*(L'*L)*dm ...
    + epsilon1*betta1*(LReg'*((LReg*dm)./(TVvec2 + epsilon1).^2)) ...
    + p*alphaTV*LReg'*((LReg*dm).*((TVvec2 + epsilon1).^(p/2-1)));

gAbs = zeros(size(g));
omega = 1e-3*2*pi*f;

for k = 1:size(U,2)
    %g = g + real(G(U(:,k))'*V(:,k));
    g = g + omega^2 * real(conj(U(:,k)).*V(:,k));
    gAbs = gAbs + real(abs(U(:,k)).*abs(U(:,k)));
end

%% project g

gMax2 = max(gAbs);

thresh = 0;%gMax2/5;

g(gAbs(:)<thresh) = NaN;

g2d = reshape(g, model.n);

%% fix the borders
% horizontally
edgeExp = 7; % 10
g2d(:,1:edgeExp) = NaN;
g2d(:,end-edgeExp+1:end)=NaN;


% vertically
edgeExp = 3;
g2d(1:edgeExp, :) = NaN;
g2d(end-edgeExp+1:end, :)=NaN;


g2d = fillmissing(g2d,'nearest',2);
g2d = fillmissing(g2d,'nearest');

%imagesc(g2d);
%drawnow

g = g2d(:);


% 
% figure(111);
% imagesc(gAbs2d);
% drawnow;
% 
% figure(112);
% imagesc(g2d);
% drawnow;
end