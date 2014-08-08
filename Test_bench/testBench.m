%testBench.m creates numerical accuracy plots
%author: Alice Bates
%date: July 2014

% Install SHT 
install_sht


L = 5; 

%% Generate random signal
% flmt = zeros(1,L^2);
% index = [1,0,5,0,17];%[1,0,5,0,17,0,37];
% for l = 1:2:5
%      flmt(index(l):(index(l) +(2*l-2) )) =  rand(1,(2*l-1)) + 1i*rand(1,(2*l-1));
% end
% save('Rand_Antipodal_Sig_L5.mat','flmt');

load('Rand_Antipodal_Sig_L5.mat');

%% inverse transform 
ft  = nsht_inverse(flmt,L);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TESTING FORWARD

% Install SHT 
%install_sht


%L = 5; 


%ft = rand(1,(15)) + 1i*rand(1,(15));

%forward
[flmr T_mat_inv] = nsht_forward(ft,L);


%% Absolute error
Error_max = max(abs(flmt-flmr));

%%MEAN error
Error_mean = sum(abs(flmt-flmr))/L^2;




