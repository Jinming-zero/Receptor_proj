function res = EEG_recep_decomp(X, E_recep,I_recep)

[num_roi, num_ts] = size(X);
num_exc = size(E_recep,2);
num_inh = size(I_recep,2);

%H = zeros(num_recept, num_ts, num_subj);
recon_r = zeros(num_ts,1);
w_exc = zeros(num_exc,num_ts);
w_inh = zeros(num_inh,num_ts);
c_exc = zeros(num_exc,num_ts);
c_inh = zeros(num_inh,num_ts);

E_map = zeros(num_roi, num_ts);
I_map = zeros(num_roi, num_ts);
E_map_norm = zeros(num_roi, 1);
I_map_norm = zeros(num_roi, 1);
X_recon = zeros(num_roi, num_ts);

for j = 1:num_ts
    temp = EI_decompose(X(:,j), E_recep, I_recep);
    recon_r(j,1) = temp.r_recon;
    w_exc(:,j) = temp.w_exc;
    w_inh(:,j) = temp.w_inh;
    c_exc(:,j) = temp.c_exc;
    c_inh(:,j) = temp.c_inh;
    E_map(:,j) = temp.E_map;
    I_map(:,j) = temp.I_map;
    X_recon(:,j) = temp.y_pred;
end

for r = 1:num_roi
    E_map_norm(r,1) = norm(E_map(r,:));
    I_map_norm(r,1) = norm(I_map(r,:));
end
    
E_I_ratio = E_map_norm ./ I_map_norm;
c_exc_average = mean(c_exc, 2);
c_inh_average = mean(c_inh, 2);

res.recon_r = recon_r;
res.w_exc = w_exc;
res.w_inh = w_inh;
res.c_exc = c_exc;
res.c_inh = c_inh;
res.E_map = E_map;
res.I_map = I_map;
res.E_map_norm = E_map_norm;
res.I_map_norm = I_map_norm;
res.E_I_ratio = E_I_ratio;

res.c_exc_average = c_exc_average;
res.c_inh_average = c_inh_average;
res.X_recon = X_recon;
end

