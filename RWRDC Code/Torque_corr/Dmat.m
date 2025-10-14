function D = Dmat(N)
% Compute the Chebyshev differentiation matrix of order N

x = chgslb(N+1);
c = [2; ones(N-1, 1); 2] .* (-1).^(0:N)';
X = repmat(x, 1, N+1);
dX = X - X';
D = (c * (1./c)')./(dX + eye(N+1));
D = D - diag(sum(D, 2));
end