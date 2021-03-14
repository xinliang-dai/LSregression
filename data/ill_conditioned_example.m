function [A,b,xref] = ill_conditioned_example(n)
    % create a symmetric ill-conditioned matrix H
    A    = hilb(n);
%     A    = A'*A;
    xref = [1:n]';
    b    = A*xref;
end