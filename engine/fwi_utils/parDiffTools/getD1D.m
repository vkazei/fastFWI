% function provides 1D differential operator
% from length n and spatial grid step h
function D1D = getD1D(h,n) 
    D1D = spdiags(ones(n,1)*[-1 1]/h, 0:1, n-1, n);        
end