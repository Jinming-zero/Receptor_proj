function [X_c, X_d, stats, stats_sub] = GSP(W, X_RS)
% W: Input network (ROI x ROI)
% X_RS: Input TS (ROI x Time_step x Subject)

%% GSP_laplacian
% data prepare
%zX_RS=zscore(X_RS,0,2);
zX_RS = X_RS;
% Number of regions
n_ROI = size(W,1);
% number of subjects
nsubjs_RS=size(zX_RS,3);

% calculate laplacian
D = diag(sum(W,2));
Wsymm = D^(-1/2) * W * D^(-1/2);
Wnew = Wsymm;
L = eye(n_ROI) - Wnew;

% laplacian decomposition
[U,LambdaL] = eig(L);
[LambdaL, IndL]=sort(diag(LambdaL));
U=U(:,IndL);

% compute cut-off frequency
for s = 1:nsubjs_RS
    X_hat_L = U' * zX_RS(:,:,s);
end

pow = abs(X_hat_L).^2;
PSD = squeeze(mean(pow,2));
mPSD = mean(PSD,2);
AUCTOT = trapz(mPSD(1:n_ROI));


i=0;
AUC=0;
while AUC<AUCTOT/2
    AUC=trapz(mPSD(1:i));
    i=i+1;
end
NN=i-1; %CUTOFF FREQUENCY C : number of low frequency eigenvalues to consider in order to have the same energy as the high freq ones
NNL=n_ROI-NN;

% split structural harmonics in high/low frequency
M=fliplr(U); %Laplacian eigenvectors flipped in order (high frequencies first)
Vlow=zeros(size(M));
Vhigh=zeros(size(M));

Vhigh(:,1:NNL)=M(:,1:NNL);%high frequencies= decoupled 
Vlow(:,end-NN+1:end)=M(:,end-NN+1:end);%low frequencies = coupled 


%% Graph Signal Analysis

data=zX_RS;

%% compute fMRI HF/LF portions
clear X_hat X_c X_d X_m N_c N_d X_all 
for s=1:size(data,3)
    X_hat(:,:,s)=M'*data(:,:,s);
    X_c(:,:,s)=Vlow*X_hat(:,:,s);
    X_d(:,:,s)=Vhigh*X_hat(:,:,s);
    X_all(:,:,s)=M*X_hat(:,:,s); % reconstruct back full signal to check norm
    
    %% norms  of the weights
    for r=1:n_ROI
        N_c(r,s)=norm(X_c(r,:,s));
        N_d(r,s)=norm(X_d(r,:,s));
    end
end

% mean across subjects
mean_c=mean(N_c,2); %average coupling
mean_d=mean(N_d,2); %average decoupling

stats.mean_c = mean_c;
stats.mean_d = mean_d;

stats_sub.mean_c = N_c;
stats_sub.mean_d = N_d;
