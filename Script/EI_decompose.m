function res = EI_decompose(y, X_exc, X_inh)
% TWO-BLOCK PLS/CCA + NNLS, returning per-receptor loadings and E/I metrics
% 
% INPUT
%   y      : N x 1 response
%   X_exc  : N x Pe (excitatory receptor maps)
%   X_inh  : N x Pi (inhibitory receptor maps)
%
% OUTPUT (struct res):
%   res.t_exc, res.t_inh           : N x 1 block scores (zscored)
%   res.r_exc_weights, res.r_inh_weights : per-receptor weights mapping receptor -> block score
%   res.alpha, res.beta            : Stage-2 coefficients (alpha >=0, beta <=0)
%   res.per_receptor_contrib_exc   : Pe x 1, excitatory receptor contribution coefficients to y
%   res.per_receptor_contrib_inh   : Pi x 1, inhibitory receptor contribution coefficients to y
%   res.EI_ratio_coef              : alpha / |beta|
%   res.EI_ratio_score             : (alpha*||t_exc||) / (|beta|*||t_inh||)
%   res.EI_balance_score           : alpha*std(t_exc) - |beta|*std(t_inh)
%   res.EI_ratio_sum               : sum(per_receptor_contrib_exc)/sum(abs(per_receptor_contrib_inh))
%   res.y_pred, res.r              : prediction and Pearson r
%
% Note: requires Statistics Toolbox (plsregress, canoncorr) and Optimization Toolbox (lsqnonneg).
num_exc = size(X_exc,2);
num_inh = size(X_inh,2);

X = [X_exc, X_inh];
w =lsqminnorm(X, y);
y_pred = X * w;
r_recon = corr(y,y_pred);
vars = explained_variance(y, y_pred);

idx_exc = 1:num_exc;
idx_inh = num_exc+1:num_exc+num_inh;

w_exc = w(idx_exc);
w_inh = w(idx_inh);

c_exc = abs(w_exc / sum(abs(w_exc)));
c_inh = abs(w_inh / sum(abs(w_inh)));

EI_coef = sum(abs(w_exc)) / sum(abs(w_inh));
E_map = X_exc * w_exc;
I_map = X_inh * w_inh;

res.w = w;
res.r_recon = r_recon;
res.vars = vars;
res.w_exc = w_exc;
res.w_inh = w_inh;
res.c_exc = c_exc;
res.c_inh = c_inh;
res.EI_coef = EI_coef;
res.E_map = E_map;
res.I_map = I_map;
res.y_pred = y_pred;
