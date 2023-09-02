function y = laplacian(x,m,b,l)

% laplace function to fit cross-correlograms
% Abhi

    y = A.*exp(-l.*abs((x-m))./b).*(l/2*b);

end