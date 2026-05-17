function [X, baseline] = gmsEDA(signal)
%%%%%%%%% INPUT %%%%%%%%%%%%
%%% signal: n by m matrix, each column is an EDA signal

%%%%%%%%% OUTPUT %%%%%%%%%%%%
% 5/16/2026
% by Xuemei Chen


num_of_chops = 5;
overlap = 0.5; 
Y = chop_data(signal, num_of_chops, overlap);
[N, K] = size(Y); % N is new signal length after chopping

%%% parameter for generating filter matrix H
tau1 = 2; 
tau2 = 0.75; 
H = make_H(N, tau1, tau2, 'tril');

%%% decompose via GMS-P
para.rho_outer = 1; 
para.lasso_rho = 1; 
para.max_iter = 100;
para.lasso_max_iter = 50;
para.tol_outer = 1e-8;
para.lasso_tol = 1e-5;
para.lasso_decomp = "svd";      
para.lasso_method = 'FISTA';
lam = 1/sqrt(max(N, K))*3; 


output = gen_matrix_sep_con(Y, H, lam, para);
S = output.S;
S(S<0) = 0;
B = Y - H*S;
X = reconstr_data(S, size(signal), num_of_chops, overlap);
baseline = reconstr_data(B, size(signal), num_of_chops, overlap);


end


%%%%%%%%%%%% helper functions %%%%%%%%%%%%%%
function Y = chop_data(data, num_of_chops, overlap)
% chop_data slices each column of input into segments with optional overlap
% This is processed per column, as different col represents different participant.
% Inputs:
% - data := input matrix 
% - num_of_chops := number of segments per column if overlap = 0
% - overlap := overlap ratio, at most 0.5
%
% Returns:
% - Y: matrix of chopped data 
% rewritten on 8/13/2025 as overlap was not working properly

arguments
    data
    num_of_chops = 2
    overlap = 0
end

[num_rows, num_cols] = size(data);

% num_of_chops will determine the length of columns in Y
Y_nrows = ceil(num_rows/num_of_chops);

overlap_length = floor(Y_nrows * overlap);
nonlap_length = Y_nrows - overlap_length;
temp = floor((num_rows - Y_nrows)/nonlap_length);
if temp < (num_rows - Y_nrows)/nonlap_length
    multiplier = temp + 2;
else
    multiplier = temp + 1;
end

Y = zeros(Y_nrows, num_cols*multiplier);
for col = 1:num_cols
    y = zeros(Y_nrows, multiplier);
    if temp < (num_rows - Y_nrows)/nonlap_length
        for i = 1:(multiplier-1)
            start_idx = nonlap_length*(i-1)+1;
            end_idx = nonlap_length*(i-1)+ Y_nrows;
            y(:,i) = data(start_idx:end_idx, col);
        end
        % last col of y is different
        y(:,multiplier) = data((num_rows-Y_nrows+1):num_rows, col);
    else
        for i = 1:(multiplier)
            start_idx = nonlap_length*(i-1)+1;
            end_idx = nonlap_length*(i-1)+ Y_nrows;
            y(:,i) = data(start_idx:end_idx, col);
        end
    end
    s_idx = (col-1)*multiplier+1;
    e_idx = (col)*multiplier;
    Y(:, s_idx:e_idx) = y;

   
    
end

end

function H = make_H(n, tau1, tau2, mode)
%%%%%%%%%% Inputs 
% n: signal length, need to be bigger than 160
% tau1 > tau2 > 0

%%%%%%%%%% output 
% H: n by n matrix

% Using model from Compressed Sensing paper, construct H
% updated 12/18/2025
% Xuemei Chen

T = 40; % duration in seconds
scale = 1;
TT = T/scale;
fs = 4; % signal sampled at every 1/fs second
f_len = TT*fs;
u = linspace(0, TT, f_len);

uu = u*scale;
f = 2 * ((exp(-uu / tau1)) - (exp(-uu / tau2))); % 160 coordinates

ncols_main_H = n + f_len - 1;
h = zeros(n + f_len - 1, 1);
h(1:f_len) = f;

% Construct the Toeplitz matrix
hflip = fliplr(h');
r = circshift(hflip, 1);
H_temp = toeplitz(h, r); 

main_H = H_temp(:,1:n); % n+f_len-1 by n

switch mode
    case 'top'
        H = main_H(1:n,:);
    case 'tril'
        H = main_H(1:n,:);
        H = tril(H);
    case 'middle'
        center_idx = floor(ncols_main_H/2);
        half_n = floor(n/2);
        H = main_H((center_idx - half_n):(n + center_idx - half_n -1), :);
end

end

function output = gen_matrix_sep_con(M, H, lam, para)
% preconditioned generalized matrix separation
% min||L||_* + lam||S||_1, subject to, L + HS = M
% input: M: given matrix: m x n
% input: H: given matrix: m x p, not necessarily circulant
% input: lam
% input: para has fields: 
%        - rho_outer 
%        - rho_inner
%        - N_outer
%        - N_inner
%        - tol_outer
%        - tol_inner
%        - decomp: 'svd' or 'chol'

% output has fields
% - L: m x n
% - S: p x n
% - count_outer
% - para
% created on 7/18/2025, Xuemei Chen
% updated on 8/5/2025


[u, s, v] = svd(H);
k = rank(H);
uk = u(:,1:k); % m by k
vk = v(:,1:k); % p by k

s = diag(s); % a vector
s = s(1:k);
s_inv = 1./s;

C = uk*diag(s_inv)*uk';

CH = uk*vk';
para.preconditioned = true;
para.V_H = v;

[~, p] = size(H);
ss = zeros(p, 1);
ss(1:k) = 1;
para.S2_H = ss;
outputC = gen_matrix_sep(C*M, CH, lam, para);


output.S = outputC.S;
output.L = M - H*outputC.S;
output.count_outer = outputC.count_outer;
%output.isCirc = outputC.isCirc;
output.para = outputC.para;

end

function output = gen_matrix_sep(M, H, lam, para)
% generalized matrix separation
% use ADMM to solve the following
% min||L||_* + lam||S||_1, subject to, L + HS = M
% input: 
%   M: given matrix: m x n
%   H: given matrix: m x p, not necessarily circulant
%   lam: positive scalar
%   para has fields: 
%       .rho_outer: step size for ADMM, provided by user
%       .max_iter: max number of iterations for ADMM
%       .tol_outer: convergence tolerance for ADMM
%       .preconditioned: true or false. If true, H has been preconditioned
%       with all singular values to be 1
%       .lasso_decomp: {'svd', 'chol'}
%       .lasso_method: {'ADMM', 'FISTA'}
%       .lasso_rho: step size for lasso if lasso_method is ADMM
%       .lasso_max_iter: max num of iterations for lasso
%       .lasso_tol: convergence tolerance for lasso

% para.preconditioned needs to be provided by user!

% output has fields
% - L: m x n
% - S: p x n
% - count_outer
% - para: similar to input with possibly more fields
% updated on 4/9/2025
% updated on 4/21/2025
% updated on 7/5/2025: precompute A^Tb outside of loop
% updated on 8/5/2025: use my_lasso() and take advantage of V if
% preconditioning H

[m, n] = size(M);
[~, p] = size(H);
%% default parameter values
rho_outer_default = 1;
max_iter_default = 200;
tol_outer_default = 1e-7;
lasso_method_default = 'FISTA';
lasso_decomp_default = 'svd';
lasso_rho_default = 1;
lasso_tol_default = 1e-5;
%% pass the parameters
if ~isfield(para, 'preconditioned')
    error('Please specify "para.preconditioned" value.')
end

if ~isfield(para,'rho_outer')
    para.rho_outer = rho_outer_default;        
end

if ~isfield(para, 'max_iter')
    para.max_iter = max_iter_default;    
end

if ~isfield(para, 'tol_outer')
    para.tol_outer = tol_outer_default;   
end

if ~ isfield(para, 'lasso_method')
    para.lasso_method = lasso_method_default;   
end

if ~ isfield(para, 'lasso_decomp')
    para.lasso_decomp = lasso_decomp_default; 
end

if ~isfield(para, 'lasso_tol')
    para.lasso_tol = lasso_tol_default;
end

rho_outer = para.rho_outer;
N_outer = para.max_iter;
tol_outer = para.tol_outer;
lasso_method = para.lasso_method;
lasso_decomp = para.lasso_decomp;
%% parameters for my_lasso()
para_lasso.max_iter = para.lasso_max_iter;
para_lasso.tol = para.lasso_tol;
para_lasso.method = lasso_method;
para_lasso.decomp = lasso_decomp; % default value

% para_lasso.isCirculant is determined later

if strcmp(lasso_method, 'FISTA')
    % FISTA
    % H circulant does not matter
    isCirculant = [];
    if para.preconditioned
        para_lasso.L = 1; % this is biggest singular value of H'H
    else
        s = svd(H'*H);
        s = diag(s);
        para_lasso.L = s(1);
    end
elseif strcmp(lasso_method, 'ADMM')
    % ADMM
    if ~isfield(para, 'lasso_rho')
        % if para.lasso_rho not given
        para.lasso_rho = lasso_rho_default;
    end
    para_lasso.rho = para.lasso_rho; 
    if is_circulant(H)
        isCirculant = true; 
        %para_lasso.isCirculant = isCirculant; 
        d = abs(fft(H(:,1)));
        para_lasso.coef = d.^2;
    else
        % more common case for ADMM
        isCirculant = false;
        %para_lasso.isCirculant = isCirculant;
        if (para_lasso.decomp == "svd") && para.preconditioned
            % singular values of H'H
            para_lasso.S2_A = para.S2_H; % p by 1
            % singular vectors of H or H'H
            para_lasso.V_A = para.V_H; % p by p
        elseif (para_lasso.decomp == "svd") && (~para.preconditioned)
            [para_lasso.V_A, para_lasso.S2_A] = pref(H);
        elseif para_lasso.decomp == "chol"
            L_H = chol(H'*H + para.lasso_rho*eye(p), 'lower');
            para_lasso.L = L_H;
        end

            
    end
end
para_lasso.isCirculant = isCirculant;
para.isCirculant = isCirculant;

%% initialization and main loop
L = zeros(m, n);
S = zeros(p, n);
U = zeros(m, n);

RelChg = 1;
count_outer = 0;
while RelChg > tol_outer && count_outer < N_outer
    Slast = S;
    Llast = L;
    L = SVT(M - H*S - U, 1/rho_outer);
    %L(L<0) = 0; %%%%%%%%%%%%%%%%%%%%%%%% this line is added on 10/12/2025
    [S, ~] = my_lasso(H, M - U - L, lam/rho_outer, para_lasso);  
    %S(S<0) = 0; %%%%%%%%%%%%%%%%%%%%%%%% this line is added on 8/12/2025

    U = U + L + H*S - M;
    count_outer = count_outer + 1;

    % Check convergence
    % finding the relative error
    Ldn = norm(L - Llast, 'fro');
    Sdn = norm(S - Slast, 'fro');
    Ln = norm(Llast, 'fro');
    Sn = norm(Slast, 'fro');

    % updating stopping critera
    RelChg = (Ldn^2 + Sdn^2)^0.5 / ((Ln^2 + Sn^2)^0.5 + 1);
end

output.L = L;
output.S = S;
output.count_outer = count_outer;
output.para = para;
end

function B = SVT(A,a)
[U,S,V] = svd(A,'econ');
S2 = SoftThresh(S,a);
B = U*S2*V';
end

function z = SoftThresh(x, kappa)

% soft thresholding operator, works for vectors or matrices
% z(i)=x(i) - kappa if x(i)>kappa
% z(i)=x(i) + kappa if x(i)<-kappa
% z(i)=0 otherwise

    z = max( 0, x - kappa ) - max( 0, -x - kappa );
end

function [x, num_iter] = my_lasso(A, b, lam, para)
% returns argmin_x {0.5|Ax-b|^2 + lam*|x|_1}
% input:
%   A: m x p
%   b: m x k
%   lam: positive scalar
%   para:
%       .method: 'ADMM', 'ISTA', 'FISTA'
%       .tol: stopping criteria
%       .max_iter
%       .V_A
%       .S2_A
%       .L_A
%       .rho: step size if method is ADMM
%       .L: lipschitz constant if method is FISTA
%       .decomp: 'svd' or 'chol'
%       .isCirculant: whether A is circulant. only used when gen_matrix_sep
%       is called.
% output:
%   .x: (p x k) 
%   .num_iter: number of iterations ran

%%%%%%%% notes %%%%%%%%%%%
% for para.method: choose ADMM or FISTA. ISTA is slow in general.

% reference for FISTA: Beck, Amir, and Marc Teboulle. 
% "A fast iterative shrinkage-thresholding algorithm for linear inverse problems." 
% SIAM journal on imaging sciences 2.1 (2009): 183-202.

% written by Xuemei Chen 8/4/2025

[~, p] = size(A);
[~, k] = size(b);
%% pass parameter
% default
max_iter_default  = 100;
tol_default       = 1e-5;
method_default    = 'FISTA';
decomp_default = 'svd';
% for ADMM, can prefactor A'A using svd or 
% prefactor A'A + rho using cholesky
rho_default = 1;

if ~isfield(para, 'max_iter')   
    para.max_iter = max_iter_default; 
end

if ~isfield(para,'tol')   
    para.tol = tol_default ; 
end

if ~isfield(para, 'method')
    para.method = method_default;
end


max_iter = para.max_iter;
tol = para.tol;
method = para.method;


if strcmp(method, 'ADMM')
    if ~isfield(para, 'rho')
        para.rho = rho_default;
    end
    
    if ~isfield(para, 'decomp')
        para.decomp = decomp_default;
    end
    if ~isfield(para, 'isCirculant')
        % if not provided, evaluate
        para.isCirculant = is_circulant(A);
    end
    rho = para.rho;
    decomp = para.decomp;
    if para.isCirculant && (~isfield(para, 'coef'))
        % if A is circulant and para.coef not provided, evaluate
        d = abs(fft(A(:,1)));
        para.coef = d.^2;
    end

    if para.isCirculant
        coef = para.coef;
    else % not circulant
        if strcmp(decomp, 'svd') && (~isfield(para, 'V_A'))
        % if use svd and V_A not provided, evaluate
        [para.V_A, para.S2_A] = pref(A);
        end
        if strcmp(decomp, 'svd')
            % if use svd
            V = para.V_A;
            sig = para.S2_A;
        end
        if strcmp(decomp, 'chol') && (~isfield(para, 'L_A'))
            % if use chol and L_A not provided, evaluate
            para.L_A = chol(A'*A + rho*eye(p), 'lower');
        end
        if strcmp(decomp, 'chol')
            % if use chol
            L_A = para.L_A;            
        end

    end
    
elseif strcmp(method, 'FISTA')
    if ~isfield(para, 'L')
        % if para.L not provided, evaluate
        s = svd(A'*A);
        s = diag(s);
        para.L = s(1);        
    end
    L = para.L;
end


%% main part
x = zeros(p, k);

switch method
    case 'ADMM'
        Atb = A'*b; % m by k, Atb precomputed        
        z = zeros(p, k);
        u = zeros(p, k);
        for j = 1:max_iter              
            xlast = x;

            % update x
            rhs = Atb + rho*(z - u); 
            if para.isCirculant                 
              % [m,n]./[m,1] will divide each col
                x = ifft(fft(rhs)./(coef + rho));
                x = real(x);
            elseif strcmp(decomp, 'svd')
                x = (V'*rhs)./(sig + rho);
                x = V*x;
            elseif strcmp(decomp, 'chol')
                y = L_A \ rhs;  
                x = L_A' \ y; 
            end

            % update z
            z = SoftThresh(x + u, lam/rho);

            % update u
            u = u + x - z;

            % check convergence
            x_change = norm(x - xlast, 'fro')/(norm(xlast, 'fro') + 1);
            if x_change < tol
                break
            end
              
        end
        num_iter = j;
    case 'ISTA'
        
        s = svd(A'*A);
        s = diag(s);
        s = s(1);
        t = 1/s;
        for j = 1:max_iter
            xlast = x;
            y = x - t*A'*(A*x - b);
            x = SoftThresh(y, lam*t);
            x_change = norm(x - xlast, 'fro')/(norm(xlast, 'fro') + 1);
            if x_change < tol   
                break
            end
        end
        num_iter = j;
    case 'FISTA'
        
        
        t = 1;
        y = x;
        for j = 1:max_iter
            xlast = x;    
            tlast = t;
            x = SoftThresh(y - (1/L)*A'*(A*y - b), lam/L);
            
            t = (1+sqrt(1+4*t^2))/2;
            y = x + (tlast - 1)/t*(x - xlast);
            x_change = norm(x - xlast, 'fro')/(norm(xlast, 'fro') + 1);
            if x_change < tol   
                break
            end
        end
        num_iter = j;
        

end

end

%% helper functions
function [V, sig] = pref(A)
% [~, sig, V] = svd(A'*A)

[p, m] = size(A);

if m < p
   [~, sig, V] = svd(A'*A); % V is m by m       
   sig = diag(sig); % m by 1
            
else
            [~, s, V] = svd(A); % V is m by m 
            s = diag(s); % p by 1
            sig = zeros(m, 1);
            sig(1:length(s)) = (abs(s)).^2;
end
end

function flag = is_circulant(A)
    [m, n] = size(A);
    
    % we need A to be square if circulant
    if m ~= n
        flag = false;
        return;
    end

    first_row = A(1, :);
    
    % check if each row is a shift of the first row
    for i = 2:m
        expected_row = circshift(first_row, [0, i-1]);
        if any(A(i, :) ~= expected_row)
            flag = false;
            return;
        end
    end
    flag = true;
end

function data = reconstr_data(Y, data_size, num_of_chops, overlap)
% this is the reverse of chop_data()
% take the average at overlap area
% created on 8/16/2025 by Xuemei Chen

%[num_rows, num_cols] = data_size;
num_rows = data_size(1);
num_cols = data_size(2);

[Y_nrows, Y_ncols] = size(Y);

%{
overlap_length = floor(Y_nrows*overlap);
nonlap_length = Y_nrows - overlap_length;
temp = floor((num_rows - Y_nrows)/nonlap_length);
if temp < (num_rows - Y_nrows)/nonlap_length
    multiplier = temp + 2;
else
    multiplier = temp + 1;
end
%}

multiplier = ceil(Y_ncols/num_cols);

idx = chop_data((1:num_rows)', num_of_chops, overlap);

data = zeros(data_size);

%DIFF = idx(1,2) - idx(1,1);
idx_diff = idx(1,:) - idx(1, 1);

for ci = 1:num_cols
    j = 1;
    
    nonoverlap_idx = setdiff(idx(:,j), idx(:,j+1));
    data(nonoverlap_idx, ci) = Y(nonoverlap_idx,j+(ci-1)*multiplier);

    overlap_idx = intersect(idx(:,j), idx(:,j+1));
    data(overlap_idx, ci) = 0.5*(Y(overlap_idx-idx_diff(j),j+(ci-1)*multiplier) + Y(overlap_idx-idx_diff(j+1),j+1+(ci-1)*multiplier));

    for j = 2:(multiplier-1)
        nonoverlap_idx = setdiff(idx(:,j), idx(:,j+1));
        % discard if already covered by previous j
        nonoverlap_idx2 = setdiff(nonoverlap_idx, idx(:,j-1));
        data(nonoverlap_idx2, ci) = Y(nonoverlap_idx2-idx_diff(j),j+(ci-1)*multiplier);

        overlap_idx = intersect(idx(:,j), idx(:,j+1));
        %keyboard
        data(overlap_idx, ci) = 0.5*(Y(overlap_idx-idx_diff(j),j+(ci-1)*multiplier) + Y(overlap_idx-idx_diff(j+1),j+1+(ci-1)*multiplier));
    
    end
    j = multiplier;
    need_idx = setdiff(idx(:,j), idx(:,j-1));

    
    data(need_idx, ci) = Y(need_idx-idx_diff(j),j+(ci-1)*multiplier);
    
end
end









