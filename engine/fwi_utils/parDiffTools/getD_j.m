% this function provides partial derivative operator matrix
% w.r.t. j_th direction taking as inputs 
% vector of grid spacing h and field size n
function D_j = getD_j(h,n,j)
    D_j = getD1D(h(j),n(j));
    if j == 1
        D_j = kron(speye(prod(n)/n(j)),D_j);
    % this should work only for 3-D and was not tested
    elseif (j > 1) && (j < length(n)) 
        D_j = kron(kron(speye(prod(n(j+1:end))),D_j),speye(prod(n(1:j-1))));
    elseif j == length(n)
        D_j  = kron(D_j,speye(prod(n)/n(j)));
    end
end