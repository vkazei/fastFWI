function L = getL(h,n)
% FD gradient operator
%
% use:
%   L = getL(h,n)
%
% input:
%   h,n   - gridspacing and number of gridpoints
%
% output
%   L     - sparse matrix

L = getD_j(h,n,1);
for j=2:length(n)
    
    D_j = getD_j(h,n,j);
    L = [L; D_j];
    
end

end