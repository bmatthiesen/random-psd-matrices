function A = randpdm(dim, trace, num, type, method)
% RANDPDM randomly generates positive definite matrices with constant or
%        bounded trace according to a uniform distribution.
%
%   A = RANDPDM(dim, trace, num, type, method)
%   
%   Inputs:
%       dim    : Matrix Dimension.
%       trace  : Trace of the matrices. Either a scalar for a constant trace
%                or a vector with [lb ub] for a bounded trace with
%                lb < tr(A) <= ub.
%       num    : Number of matrices to generate.
%       type   : Optional. Either 'real' or 'complex'. Defaults to 'real'.
%       method : Optional. Used to select generation algorithm. Either
%                'rejection' or 'betadistr'. Defaults to 'rejection'.
%                Note that 'betadistr' requires the Statistics Toolbox.
%
%   Output:
%       3-dimensional array containing the generated matrices. First two
%       dimensions are the row and column of each matrix, repspectively.
%       The third dimension indexes the individual matrices. See EXAMPLES
%       below.
%
%   Examples:
%       A = randpdm(3, 1, 10) returns ten real 3x3 positive definite matrices
%       with unit trace drawn according to a uniform distribution using the 
%       rejection method. The i-th matrix in the generated set is accessed
%       with A(:,:,i). 
%
%       A = randpdm(2, [1 2], 100, 'complex','betadistr') returns 100 
%       complex 2-dimensional positive definite matrices with trace in the 
%       interval (1,2] drawn according to a uniform distribution using the 
%       method of transforming a beta distribution.
%
%   See also gallery

%   References:
%   [1] M. Mittelbach, B. Matthiesen, E. A. Jorswieck, "Sampling Uniformly 
%       From the Set of Positive Definite Matrices With Trace Constraint," 
%       accepted for publication in IEEE Transactions on Signal Processing,
%       Jan. 2012.
%
%   Bho Matthiesen <bho.matthiesen@gmail.com>
%   Copyright 2010-2012 Bho Matthiesen, Martin Mittelbach, Eduard A. Jorswieck
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU Lesser General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version.
%
%   This program is distributed in the hope that it will be useful, but WITHOUT
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
%   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
%   License for more details.
%
%   You should have received a copy of the GNU Lesser General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.


% setup
if nargin > 4
    switch lower(method)
        case 'rejection'
            rejection = true;
        case 'betadistr'
            rejection = false;
        otherwise
            error('Unknown method ''%s''. Use either ''rejection'' or ''betadistr''.', method)
    end
else
    rejection = true;
end

if nargin > 3
    switch lower(type)
        case 'real'
            complex = false;
        case 'complex'
            complex = true;
        otherwise
            error('Unknown type ''%s''. Use either ''real'' or ''complex''.', type);
    end
else
    complex = false;
    type = 'real';
end

if complex
    [h, g] = zf_complex(dim);
else
    [h, g] = zf_real(dim);
end
phi_n = length(h);

switch length(trace)
    case 1
        tau = trace * ones(1,num);
    case 2
        if trace(1) == trace(2)
            tau = trace * ones(1,num);
        else
            tau = rndtrace(dim, trace(1), trace(2), num, type);
        end
    otherwise
        error('Trace has to be a skalar or [ lb ub ]!')
end

% generate
A = zeros(dim, dim, num);
if rejection
    for ii = 1:num
        phi = zeros(1, phi_n);
        
        for l = 1:phi_n
        if h(l) == 0
            sigma = sqrt(1/g(l));
            
            while 1
            Xn = randn(1) .* sigma + pi/2;
            Un = rand(1);
            
            if 0 <= Xn && Xn <= pi
                tmp = sin(Xn)^g(l);
            else
                tmp = 0;
            end
            
            if Un <= tmp * exp((Xn-pi/2)^2 / (2/g(l)))
                break
            end
            end
        else
            mu = atan(sqrt(g(l)/h(l)));
            sigma = 1/(sqrt(h(l)) + sqrt(g(l)));
            sigmasq = sigma^2;
            
            a = sqrt(1+g(l)/h(l));
            b = sqrt(1+h(l)/g(l));
            
            while 1
            Xn = randn(1) .* sigma + mu;
            Un = rand(1);
            
            if 0 <= Xn && Xn <= pi/2
                tmp = (a * cos(Xn))^h(l) * (b * sin(Xn))^g(l);
            else
                tmp = 0;
            end
            
            if Un <= tmp * exp((Xn-mu)^2 / (2*sigmasq))
                break
            end
            end
        end
        
        phi(l) = Xn;
        end

        A(:,:,ii) = tau(ii) .* makeA(dim, phi, complex);
    end
else % betadistr
    for ii = 1:num
        % generate random variates with beta distribution
        gam_a = randg((g+1)./2);
        gam_b = randg((h+1)./2);
        y = gam_a ./ (gam_a + gam_b);
        
        phi = asin(sqrt(y));

        % bernoulli distributed random variates (coin flip)
        bino = rand(1, phi_n) <= .5;
        
        t = (h == 0 & bino);
        phi(t) = pi - phi(t);

        A(:,:,ii) = tau(ii) .* makeA(dim, phi, complex);
    end
end
end


function A = makeA(dim, phi, complex)
T = zeros(dim);
x = zeros(length(phi)+1, 1);

l = 1;
for idx = 1:length(phi)
    x(idx) = l * cos(phi(idx));
    l = l * sin(phi(idx));
end
x(idx+1) = l;

if complex
    idx = 0;
    for m = 1:dim
        idx2 = idx + 2*m - 1;
        T(1:m,m) = x(idx+1:2:idx2) + 1i*[x(idx+2:2:idx2); 0];
        idx = idx2;
    end
else
    idx = 0;
    for m = 1:dim
        idx2 = idx + m;
        T(1:m,m) = x(idx+1:idx2);
        idx = idx2;
    end
end

A = T' * T;
end


function [c,s] = zf_real(n)

c = zeros(1,n*(n+1)/2-1);
a = c;
b = c;

k = 1:(n-1);

c(k.*(k+1)/2) = n+1-k;

for ii=k
    m = 0:ii;
    l = ii*(ii+1)/2 + m;
    a(l) = ii-1;
    b(l) = ii+1+m;
end

s = n^2 - a*n - b;

end


function [c,s] = zf_complex(n)

c = zeros(1,n^2-1);
a = c;
b = c;

k = 1:(n-1);

c(k.^2) = 2*(n-(1:(n-1)))+1;

for ii=k
    m = 0:2*ii;
    l = ii^2 + m;
    a(l) = n-ii-1;
    b(l) = (ii-1)*n+1+m;
end

s = n^2 + a*n - b;

end


function tau = rndtrace(dim, lb, ub, num, type)
switch lower(type)
    case 'complex'
        a = dim^2;
    case 'real'
        a = (dim^2 + dim) / 2;
    otherwise
        error('Unknown type ''%s''. Use either ''real'' or ''complex''.', type);
end

if lb < 0
    error('Trace may not be negative.');
end

if ub < lb
    error('Upper bound must be greater than lower bound!');
end

u = rand(1,num);

if lb == 0
    tau = ub .* u.^(1/a);
else
    tau = lb .* (((ub/lb)^a - 1) .* u + 1).^(1/a);
end
end

% vim: sw=4 ts=4 expandtab
