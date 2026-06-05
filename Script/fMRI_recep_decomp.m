function res = fMRI_recep_decomp(X,E_recep,I_recep)
[num_roi, num_ts, num_subj] = size(X);
num_exc = size(E_recep,2);
num_inh = size(I_recep,2);

%H = zeros(num_recept, num_ts, num_subj);
recon_r = zeros(num_ts, num_subj);
X_recon = zeros(num_roi, num_ts, num_subj);
w_exc = zeros(num_exc,num_ts, num_subj);
w_inh = zeros(num_inh,num_ts, num_subj);
c_exc = zeros(num_exc,num_ts, num_subj);
c_inh = zeros(num_inh,num_ts, num_subj);

E_map = zeros(num_roi, num_ts, num_subj);
I_map = zeros(num_roi, num_ts, num_subj);
E_map_norm = zeros(num_roi, num_subj);
I_map_norm = zeros(num_roi, num_subj);

for i = 1:num_subj
    X_sub = X(:,:,i);
    %H_sub = zeros(num_recept, num_ts);
    for j = 1:num_ts
        X_ts = X_sub(:,j);
        temp = EI_decompose(X_ts, E_recep, I_recep);
        recon_r(j,i) = temp.r_recon;
        c_exc(:,j,i) = temp.c_exc;
        c_inh(:,j,i) = temp.c_inh;
        E_map(:,j,i) = temp.E_map;
        I_map(:,j,i) = temp.I_map;
        X_recon(:,j,i) = temp.y_pred;
    end
    
    for r = 1:num_roi
        E_map_norm(r,i) = norm(E_map(r,:,i));
        I_map_norm(r,i) = norm(I_map(r,:,i));
    end
end

E_I_ratio = E_map_norm ./ I_map_norm;

c_exc_subj = squeeze(mean(c_exc, 2));
c_exc_average = mean(c_exc_subj,2);
c_inh_subj = squeeze(mean(c_inh, 2));
c_inh_average = mean(c_inh_subj,2);
E_I_ratio_average = mean(E_I_ratio,2);


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


res.c_exc_subj = c_exc_subj;
res.c_exc_average = c_exc_average;
res.c_inh_subj = c_inh_subj;
res.c_inh_average = c_inh_average;
res.E_I_ratio_average = E_I_ratio_average;
res.X_recon = X_recon;




end



