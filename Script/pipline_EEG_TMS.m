
%% GSP for TMS-EEG
clear all;clc

data = importdata('D:\DATA\TMS-EEG\important\my_process\Net_Ts_path\tms-eeg-activity.mat');
evoke = data.evoke;
evoke_ts = evoke(:,1000:1300);
plot(evoke_ts')

sc = importdata('D:\DATA\TMS-EEG\important\my_process\Structure_conn\structure_net_raw.mat');
sc_spa = sc.sc;
imagesc(sc_spa)

[X_c, X_d, stats] = GSP(sc_spa, evoke_ts);
subplot(3,1,1), plot(evoke_ts');
subplot(3,1,2), plot(X_c');
subplot(3,1,3), plot(X_d');

subplot(2,1,1), scatter(1:200,stats.mean_c);
subplot(2,1,2), scatter(1:200,stats.mean_d);

save_path = 'D:\Project\Recepter-GSP\process\EEG_TMS\GSP';
temporal_file = 'D:\DATA\HCP_DATA\template_and_method\Schaefer2018_LocalGlobal_new\Schaefer2018_200Parcels_7Networks_order_neo.dlabel.nii';

roi2surf(stats.mean_c, temporal_file, fullfile(save_path,'GSP_Low.dscalar.nii'));
roi2surf(stats.mean_d, temporal_file, fullfile(save_path,'GSP_High.dscalar.nii'));

roi2surf(zscore(stats.mean_c), temporal_file, fullfile(save_path,'GSP_Low_Z.dscalar.nii'));
roi2surf(zscore(stats.mean_d), temporal_file, fullfile(save_path,'GSP_High_Z.dscalar.nii'));
save('D:\Project\Recepter-GSP\process\EEG_TMS\GSP\gsp_result.mat','X_c','X_d', 'stats','sc_spa','evoke_ts');

%% receptor model fitting
clear all;clc
load('D:\Project\Recepter-GSP\process\TMS\evoke_SC_raw\gsp_result.mat');
load('D:\Project\Recepter-GSP\process\receptor\recepteor_data_200_E_I.mat');
E_recep = E_recep';
I_recep = I_recep';

EEG_res_d = EEG_recep_decomp(X_d,E_recep,I_recep);
EEG_res_c = EEG_recep_decomp(X_c,E_recep,I_recep);

recon_d = mean(mean(EEG_res_d.recon_r));
recon_c = mean(mean(EEG_res_c.recon_r));

save_path = 'D:\Project\Recepter-GSP\process\EEG_TMS\Fit\EEG_result.mat';
save(save_path, 'EEG_res_d', 'EEG_res_c');





%%














