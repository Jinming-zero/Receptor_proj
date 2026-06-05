clear all;clc

%%  Pre data
clear all; clc
Pre_Data = importdata('D:\Project\Recepter-GSP\data\fMRI-TMS\Pre.mat');

num_sub = length(Pre_Data);
Pre_Ts = [];
for i = 1:num_sub
    sub_data = Pre_Data{i, 1};
    sub_data = zscore(sub_data');
    Pre_Ts = cat(3,Pre_Ts,sub_data');
end
sc = importdata('D:\DATA\TMS-EEG\important\my_process\Structure_conn\structure_net_raw.mat');
sc_spa = sc.sc;

[X_c, X_d, stats, stats_sub] = GSP(sc_spa, Pre_Ts);
save('D:\Project\Recepter-GSP\process\fMRI_TMS\GSP\GSP_results.mat','X_c','X_d','stats','stats_sub');

load('D:\Project\Recepter-GSP\process\receptor\recepteor_data_200_E_I.mat');
E_recep = zscore(E_recep');
I_recep = zscore(I_recep');

Pre_res_d = fMRI_recep_decomp(X_d,E_recep,I_recep);
Pre_res_c = fMRI_recep_decomp(X_c,E_recep,I_recep);

save_path = 'D:\Project\Recepter-GSP\process\fMRI_TMS\Pre_GSP\Pre_GSP.mat';
save(save_path,'Pre_res_d','Pre_res_c');

load(save_path);
EI_map_c = Pre_res_c.E_I_ratio_average;
EI_map_d = Pre_res_d.E_I_ratio_average;

save_path = 'D:\Project\Recepter-GSP\process\fMRI_TMS\FIg\Pre\EI_map';
temporal_file = 'D:\DATA\HCP_DATA\template_and_method\Schaefer2018_LocalGlobal_new\Schaefer2018_200Parcels_7Networks_order_neo.dlabel.nii';

roi2surf(EI_map_c, temporal_file, fullfile(save_path,'EI_map_c.dscalar.nii'));
roi2surf(EI_map_d, temporal_file, fullfile(save_path,'EI_map_d.dscalar.nii'));

roi2surf(zscore(EI_map_c), temporal_file, fullfile(save_path,'EI_map_c_Z.dscalar.nii'));
roi2surf(zscore(EI_map_d), temporal_file, fullfile(save_path,'EI_map_d_Z.dscalar.nii'));


%%  CTR data
clear all; clc
CTR_Data = importdata('D:\Project\Recepter-GSP\data\fMRI-TMS\CTR-post.mat');

num_sub = length(CTR_Data);
CTR_Ts = [];
for i = 1:num_sub
    sub_data = CTR_Data{i, 1};
    sub_data = zscore(sub_data');
    CTR_Ts = cat(3,CTR_Ts,sub_data');
end
sc = importdata('D:\DATA\TMS-EEG\important\my_process\Structure_conn\structure_net_raw.mat');
sc_spa = sc.sc;

[X_c, X_d, stats, stats_sub] = GSP(sc_spa, CTR_Ts);
save('D:\Project\Recepter-GSP\process\fMRI_TMS\GSP\CTR_GSP.mat','X_c','X_d','stats','stats_sub');

load('D:\Project\Recepter-GSP\process\receptor\recepteor_data_200_E_I.mat');
E_recep = zscore(E_recep');
I_recep = zscore(I_recep');

CTR_res_d = fMRI_recep_decomp(X_d,E_recep,I_recep);
CTR_res_c = fMRI_recep_decomp(X_c,E_recep,I_recep);

save_path = 'D:\Project\Recepter-GSP\process\fMRI_TMS\Results\CTR_results.mat';
save(save_path,'CTR_res_d','CTR_res_c');


%%  FRO data
clear all; clc
FRO_Data = importdata('D:\Project\Recepter-GSP\data\fMRI-TMS\FRO-post.mat');

num_sub = length(FRO_Data);
FRO_Ts = [];
for i = 1:num_sub
    sub_data = FRO_Data{i, 1};
    sub_data = zscore(sub_data');
    FRO_Ts = cat(3,FRO_Ts,sub_data');
end
sc = importdata('D:\DATA\TMS-EEG\important\my_process\Structure_conn\structure_net_raw.mat');
sc_spa = sc.sc;

[X_c, X_d, stats, stats_sub] = GSP(sc_spa, FRO_Ts);
save('D:\Project\Recepter-GSP\process\fMRI_TMS\GSP\FRO_GSP.mat','X_c','X_d','stats','stats_sub');

load('D:\Project\Recepter-GSP\process\receptor\recepteor_data_200_E_I.mat');
E_recep = zscore(E_recep');
I_recep = zscore(I_recep');

FRO_res_d = fMRI_recep_decomp(X_d,E_recep,I_recep);
FRO_res_c = fMRI_recep_decomp(X_c,E_recep,I_recep);

save_path = 'D:\Project\Recepter-GSP\process\fMRI_TMS\Results\FRO_results.mat';
save(save_path,'FRO_res_d','FRO_res_c');

%%  OCC data
clear all; clc
OCC_Data = importdata('D:\Project\Recepter-GSP\data\fMRI-TMS\OCC-post.mat');

num_sub = length(OCC_Data);
OCC_Ts = [];
for i = 1:num_sub
    sub_data = OCC_Data{i, 1};
    sub_data = zscore(sub_data');
    OCC_Ts = cat(3,OCC_Ts,sub_data');
end
sc = importdata('D:\DATA\TMS-EEG\important\my_process\Structure_conn\structure_net_raw.mat');
sc_spa = sc.sc;

[X_c, X_d, stats, stats_sub] = GSP(sc_spa, OCC_Ts);
save('D:\Project\Recepter-GSP\process\fMRI_TMS\GSP\OCC_GSP.mat','X_c','X_d','stats','stats_sub');

load('D:\Project\Recepter-GSP\process\receptor\recepteor_data_200_E_I.mat');
E_recep = zscore(E_recep');
I_recep = zscore(I_recep');

OCC_res_d = fMRI_recep_decomp(X_d,E_recep,I_recep);
OCC_res_c = fMRI_recep_decomp(X_c,E_recep,I_recep);

save_path = 'D:\Project\Recepter-GSP\process\fMRI_TMS\Results\OCC_results.mat';
save(save_path,'OCC_res_d','OCC_res_c');



%% CTR vs Pre
clear all;clc
Pre_res = importdata('D:\Project\Recepter-GSP\process\fMRI_TMS\Results\Pre_results.mat');
CTR_res = importdata('D:\Project\Recepter-GSP\process\fMRI_TMS\Results\CTR_results.mat');

Pre_res_c = Pre_res.Pre_res_c;
Pre_res_d = Pre_res.Pre_res_d;

CTR_res_c = CTR_res.CTR_res_c;
CTR_res_d = CTR_res.CTR_res_d;

stat_CTR_d = group_EI_analysis(CTR_res_d, Pre_res_d);
stat_CTR_c = group_EI_analysis(CTR_res_c, Pre_res_c);

save_path = 'D:\Project\Recepter-GSP\process\fMRI_TMS\Results\CTR_vs_Pre_stats.mat';
save(save_path, 'stat_CTR_d','stat_CTR_c');

load(save_path);
CTR_c_cohend = stat_CTR_c.stat_E_I_ratio(:,3);
CTR_d_cohend = stat_CTR_d.stat_E_I_ratio(:,3);

cifti_path = 'D:\Project\Recepter-GSP\process\fMRI_TMS\Results\cifiti\CTR';
temporal_file = 'D:\DATA\HCP_DATA\template_and_method\Schaefer2018_LocalGlobal_new\Schaefer2018_200Parcels_7Networks_order_neo.dlabel.nii';

roi2surf_cohen_thr(CTR_c_cohend, temporal_file, 0.8, fullfile(cifti_path,'CTR_EI_c_cohend.dscalar.nii'),  fullfile(cifti_path,'CTR_EI_c_cohend.dlabel.nii'));
roi2surf_cohen_thr(CTR_d_cohend, temporal_file, 0.8, fullfile(cifti_path,'CTR_EI_d_cohend.dscalar.nii'),  fullfile(cifti_path,'CTR_EI_d_cohend.dlabel.nii'));



%% FRO vs Pre
clear all;clc
Pre_res = importdata('D:\Project\Recepter-GSP\process\fMRI_TMS\Results\Pre_results.mat');
FRO_res = importdata('D:\Project\Recepter-GSP\process\fMRI_TMS\Results\FRO_results.mat');

Pre_res_c = Pre_res.Pre_res_c;
Pre_res_d = Pre_res.Pre_res_d;

FRO_res_c = FRO_res.FRO_res_c;
FRO_res_d = FRO_res.FRO_res_d;

stat_FRO_d = group_EI_analysis(FRO_res_d, Pre_res_d);
stat_FRO_c = group_EI_analysis(FRO_res_c, Pre_res_c);

save_path = 'D:\Project\Recepter-GSP\process\fMRI_TMS\Results\FRO_vs_Pre_stats.mat';
save(save_path, 'stat_FRO_d','stat_FRO_c');

load(save_path);
FRO_c_cohend = stat_FRO_c.stat_E_I_ratio(:,3);
FRO_d_cohend = stat_FRO_d.stat_E_I_ratio(:,3);

cifti_path = 'D:\Project\Recepter-GSP\process\fMRI_TMS\Results\cifiti\FRO';
temporal_file = 'D:\DATA\HCP_DATA\template_and_method\Schaefer2018_LocalGlobal_new\Schaefer2018_200Parcels_7Networks_order_neo.dlabel.nii';
roi2surf_cohen_thr(FRO_c_cohend, temporal_file, 0.8, fullfile(cifti_path,'FRO_EI_c_cohend.dscalar.nii'),  fullfile(cifti_path,'FRO_EI_c_cohend.dlabel.nii'));
roi2surf_cohen_thr(FRO_d_cohend, temporal_file, 0.8, fullfile(cifti_path,'FRO_EI_d_cohend.dscalar.nii'),  fullfile(cifti_path,'FRO_EI_d_cohend.dlabel.nii'));


%% OCC vs Pre
clear all;clc
Pre_res = importdata('D:\Project\Recepter-GSP\process\fMRI_TMS\Results\Pre_results.mat');
OCC_res = importdata('D:\Project\Recepter-GSP\process\fMRI_TMS\Results\OCC_results.mat');

Pre_res_c = Pre_res.Pre_res_c;
Pre_res_d = Pre_res.Pre_res_d;

OCC_res_c = OCC_res.OCC_res_c;
OCC_res_d = OCC_res.OCC_res_d;

stat_OCC_d = group_EI_analysis(OCC_res_d, Pre_res_d);
stat_OCC_c = group_EI_analysis(OCC_res_c, Pre_res_c);

save_path = 'D:\Project\Recepter-GSP\process\fMRI_TMS\Results\OCC_vs_Pre_stats.mat';
save(save_path, 'stat_OCC_d','stat_OCC_c');

load(save_path);
OCC_c_cohend = stat_OCC_c.stat_E_I_ratio(:,3);
OCC_d_cohend = stat_OCC_d.stat_E_I_ratio(:,3);

cifti_path = 'D:\Project\Recepter-GSP\process\fMRI_TMS\Results\cifiti\OCC';
temporal_file = 'D:\DATA\HCP_DATA\template_and_method\Schaefer2018_LocalGlobal_new\Schaefer2018_200Parcels_7Networks_order_neo.dlabel.nii';
roi2surf_cohen_thr(OCC_c_cohend, temporal_file, 0.8, fullfile(cifti_path,'OCC_EI_c_cohend.dscalar.nii'),  fullfile(cifti_path,'OCC_EI_c_cohend.dlabel.nii'));
roi2surf_cohen_thr(OCC_d_cohend, temporal_file, 0.8, fullfile(cifti_path,'OCC_EI_d_cohend.dscalar.nii'),  fullfile(cifti_path,'OCC_EI_d_cohend.dlabel.nii'));






































%%

X_delta = Pre_Ts(:,2:end,:) - Pre_Ts(:,1:end-1,:);
[X_c_delta, X_d_delta, stats, stats_sub] = GSP(sc_spa, X_delta);

load('D:\Project\Recepter-GSP\process\receptor\recepteor_data_200_E_I.mat');
recept_mat = [zscore(E_recep'), -zscore(I_recep'),ones(200,1)];
recept_mat = [zscore(E_recep'), -zscore(I_recep')];
recept_mat = [E_recep', I_recep'];
recept_mat = [E_recep', I_recep',ones(200,1)];
E_recep = zscore(E_recep');
I_recep = zscore(I_recep');

X_d_delta = X_d(:,2:end,:) - X_d(:,1:end-1,:);
X_c_delta = X_c(:,2:end,:) - X_c(:,1:end-1,:);

[H_d_e, H_d_i, recon_spatial_d] = dynamic_recep_decomp(X_d,recept_mat);
[H_c_e, H_c_i, recon_spatial_c] = dynamic_recep_decomp(X_c,recept_mat);

recon_spatial_d(isnan(recon_spatial_d)) = 0;
mean(mean(recon_spatial_d))
mean(mean(recon_spatial_c))

hist(recon_spatial_d)
hist(recon_spatial_c)

recon_spatial_d = fMRI_recep_decomp(X_d,E_recep,I_recep);
recon_spatial_c = fMRI_recep_decomp(X_c,E_recep,I_recep);





temporal_file = 'D:\DATA\HCP_DATA\template_and_method\Schaefer2018_LocalGlobal_new\Schaefer2018_200Parcels_7Networks_order_neo.dlabel.nii';
roi2surf(stats.mean_c, temporal_file, fullfile(save_path, 'Mean_c.dscalar.nii'));
roi2surf(zscore(stats.mean_c), temporal_file, fullfile(save_path,'Z_Mean_c.dscalar.nii'));
roi2surf(stats.mean_d, temporal_file, fullfile(save_path, 'Mean_d.dscalar.nii'));
roi2surf(zscore(stats.mean_d), temporal_file, fullfile(save_path,'Z_Mean_d.dscalar.nii'));


mean_c = stats.mean_c;
mean_d = stats.mean_d;

%mean_c = normalize(mean_c, 'range', [0, 1]);
%mean_d = normalize(mean_d, 'range', [0, 1]);

load('D:\Project\Recepter-GSP\process\receptor\recepteor_data_200_E_I.mat');
recept_mat = [E_recep', -I_recep',ones(200,1)];

H_d = lsqnonneg(recept_mat, mean_d);
recon_d = recept_mat * H_d;
recon_spatial_corr_d = diag(corr(recon_d, mean_d));

e_idx = [1:8];
i_idx = [9:15];
E_d = H_d(e_idx,:);
I_d = H_d(i_idx,:);


H_c = lsqnonneg(recept_mat, mean_c);
recon_c = recept_mat * H_c;
recon_spatial_corr_c = diag(corr(recon_c, mean_c));
E_c = H_c(e_idx,:);
I_c = H_c(i_idx,:);



