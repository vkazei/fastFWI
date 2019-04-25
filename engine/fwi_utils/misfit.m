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
%   [f,g] = misfit(m0, m, D, opts, model)
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
P = getP(model.zr,model.xr, model.z, model.x);
Q = getQ(model.zs,model.xs, model.z, model.x);

%% forward banded matrix based solver
spparms('bandden',0.0);
U = A\Q;
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

%fprintf('Average epsilon to gradient %.3e; ',epsilon1/max(TVvec2));

MGSfunctional = norm((TVvec2./(TVvec2 + epsilon1)).^0.5)^2;

%p=opts.R.p;

TVfunctional = norm((TVvec2 + epsilon1).^(p/4))^2;

%% compute f with regularization terms
% full data
D_m = P'*U;

% crop from data
D_m = D_m.*sign(abs(D));



f_reg = .5*alpha1*norm(L*dm)^2 + .5*betta1*MGSfunctional + ...
    + .5*alphaTV*TVfunctional;


fFWI = .5*norm(D_m - D,'fro')^2;

f = f_reg + fFWI;

%% adjoint solve

%tic;
spparms('bandden',0.0);
V = A'\(P*(D - D_m));
spparms('default');
%toc;

%% compute gReg -- gradient of the regularizer
% Tikhonov regularizer
gReg.Tikh = alpha1*(L'*L)*dm;

% Minimum gradient support regularizer
gReg.MGS = epsilon1*betta1*(LReg'*((LReg*dm)./(TVvec2 + epsilon1).^2));

% Sobolev space regularizer
gReg.SS = p*alphaTV*LReg'*((LReg*dm).*((TVvec2 + epsilon1).^(p/2-1)));

gReg.Full = gReg.Tikh + gReg.MGS + gReg.SS;



%% stacking gradients from different sources
omega = 1e-3*2*pi*f;
gFWI = omega^2 * real(sum(conj(U).*V,2));

g = gFWI + gReg.Full;

regPart=norm(gReg.Full)/norm(gFWI);
fprintf('Grad norm regularization ratio (norm(gReg.Full)/norm(gFWI)) = %.3e \n',regPart);

%% basic tracking gradients and models for Pavel
if opts.tracking

load([opts.histFolder,'J'],'J')

J.m = m;
J.gFWI = reshape(gFWI, model.n);
J.evalnum = J.evalnum + 1;
J.fFWI = fFWI;

save([opts.histFolder,'J'],'J')
save([opts.histFolder,'/J_',num2str(J.evalnum)],'J')

end
%% fix the boundaries
g2d = reshape(g, model.n);

% horizontally
edgeExp = 3; % 10
g2d(:,1:edgeExp) = NaN;
g2d(:,end-edgeExp+1:end)=NaN;

% vertically
edgeExp = 3;
g2d(1:edgeExp, :) = NaN;
g2d(end-edgeExp+1:end, :)=NaN;

% fill in
g2d = fillmissing(g2d,'nearest',2);
g2d = fillmissing(g2d,'nearest');

g2d = opts.R.dvMask.*g2d;

%figure(111); imagesc(g2d);

g = g2d(:);


end