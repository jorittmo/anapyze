figure();

formula_adnimem_mmse = 'ADNI_MEM ~ TimeFromBaseliney*(Group) + BL_Age + Sex + Education + MMSE + (TimeFromBaseliney|Subject)';
formula_adnief_mmse = 'ADNI_EF ~ TimeFromBaseliney*(Group) +  BL_Age + Sex + Education + MMSE + (TimeFromBaseliney|Subject)';
formula_adnidelta_mmse = 'DMEMEF ~ TimeFromBaseliney*(Group) +  BL_Age + Sex + Education + MMSE + (TimeFromBaseliney|Subject)';
formula_mmse = 'MMSE ~ TimeFromBaseliney*(Group) +  BL_Age + Sex + Education + (TimeFromBaseliney|Subject)';

%MMSE
lme_mmse_mci = fitlme(CLINICALMMSELongitMCI,formula_mmse);
lme_mmse_ad = fitlme(CLINICALMMSELongitAD,formula_mmse);

%ADNI_MEM
lme_adnimem_mci = fitlme(CLINICALCOMPLongitMCI,formula_adnimem_mmse);
lme_adnimem_ad = fitlme(CLINICALCOMPLongitAD,formula_adnimem_mmse);

%ADNI_EF
lme_adnief_mci = fitlme(CLINICALCOMPLongitMCI,formula_adnief_mmse);
lme_adnief_ad = fitlme(CLINICALCOMPLongitAD,formula_adnief_mmse);

%ADNI_VS
lme_adnivs_mci = fitlme(CLINICALCOMPLongitMCI,formula_adnidelta_mmse);
lme_adnivs_ad = fitlme(CLINICALCOMPLongitAD,formula_adnidelta_mmse);

% Predictions for MCI

subplot(2,4,1);

% MMSE MCI

% Set AD-like
Dummy_MCI(:,'Group') = {'A+T+Asyn-'};
% Predict
[ypred,yCI,DF] = predict(lme_mmse_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','r')

% Set Mixed
Dummy_MCI(:,'Group') = {'A+T+Asyn+'};
% Predict
[ypred,yCI,DF] = predict(lme_mmse_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','g')

% Set LB
Dummy_MCI(:,'Group') = {'Asyn+'};
% Predict
[ypred,yCI,DF] = predict(lme_mmse_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','b')

Dummy_MCI(:,'Group') = {'A-T-Asyn-'};
% Predict
[ypred,yCI,DF] = predict(lme_mmse_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','m')

mmse_mci_lb = ypred
ylabel('MMSE')
xlabel('Time from baseline (y)')

subplot(2,4,2);

%ADNI_MEM MCI

% Set AD-like
Dummy_MCI(:,'Group') = {'A+T+Asyn-'};

% Predict
[ypred,yCI,DF] = predict(lme_adnimem_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','r')

Dummy_MCI(:,'Group') = {'A+T+Asyn+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnimem_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','g')

% Set sLB
Dummy_MCI(:,'Group') = {'Asyn+'};

% Predict
[ypred,yCI,DF] = predict(lme_adnimem_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','b')

% Set Negative
Dummy_MCI(:,'Group') = {'A-T-Asyn-'};

% Predict
[ypred,yCI,DF] = predict(lme_adnimem_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','m')

ylabel('ADNI MEM')
xlabel('Time from baseline (y)')

subplot(2,4,3);

% ADNI_EF MCI

% Set AD-like
Dummy_MCI(:,'Group') = {'A+T+Asyn-'};

% Predict
[ypred,yCI,DF] = predict(lme_adnief_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','r')

Dummy_MCI(:,'Group') = {'A+T+Asyn+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnief_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','g')

% Set sLB
Dummy_MCI(:,'Group') = {'Asyn+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnief_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','b')

% Set Negative
Dummy_MCI(:,'Group') = {'A-T-Asyn-'};
% Predict
[ypred,yCI,DF] = predict(lme_adnief_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','m')

ylabel('ADNI EF')
xlabel('Time from baseline (y)')

subplot(2,4,4);

% D(MEM_EF) MCI

% Set AD-like
Dummy_MCI(:,'Group') = {'A+T+Asyn-'};
% Predict
[ypred,yCI,DF] = predict(lme_adnivs_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','r')

Dummy_MCI(:,'Group') = {'A+T+Asyn+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnivs_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','g')

Dummy_MCI(:,'Group') = {'Asyn+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnivs_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','b')

Dummy_MCI(:,'Group') = {'A-T-Asyn-'};
% Predict
[ypred,yCI,DF] = predict(lme_adnivs_mci,Dummy_MCI);
% Plot
boundedline(Dummy_MCI.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','m')

ylabel('D(MEM-EF)')
xlabel('Time from baseline (y)')

% Predictions for ADD

subplot(2,4,5);

%MMSE AD

% Set AD-like
Dummy_AD(:,'Group') = {'A+T+Asyn-'};
% Predict
[ypred,yCI,DF] = predict(lme_mmse_ad,Dummy_AD);
% Plot
boundedline(Dummy_AD.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','r')


% Set Mixed
Dummy_AD(:,'Group') = {'A+T+Asyn+'};
% Predict
[ypred,yCI,DF] = predict(lme_mmse_ad,Dummy_AD);
% Plot
boundedline(Dummy_AD.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','g')


% Set LB
Dummy_AD(:,'Group') = {'Asyn+'};
% Predict
[ypred,yCI,DF] = predict(lme_mmse_ad,Dummy_AD);
% Plot
boundedline(Dummy_AD.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','b')

mmse_ad_lb = ypred

ylabel('MMSE')
xlabel('Time from baseline (y)')

subplot(2,4,6);

% Set AD-like
Dummy_AD(:,'Group') = {'A+T+Asyn-'};
% Predict
[ypred,yCI,DF] = predict(lme_adnimem_ad,Dummy_AD);
% Plot
boundedline(Dummy_AD.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','r')

% Set Mixed
Dummy_AD(:,'Group') = {'A+T+Asyn+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnimem_ad,Dummy_AD);
% Plot
boundedline(Dummy_AD.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','g')

% Set LB
Dummy_AD(:,'Group') = {'Asyn+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnimem_ad,Dummy_AD);
% Plot
boundedline(Dummy_AD.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','b')

ylabel('ADNI MEM')
xlabel('Time from baseline (y)')

subplot(2,4,7);

% Set AD-like
Dummy_AD(:,'Group') = {'A+T+Asyn-'};
% Predict
[ypred,yCI,DF] = predict(lme_adnief_ad,Dummy_AD);
% Plot
boundedline(Dummy_AD.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','r')

% Set Mixed
Dummy_AD(:,'Group') = {'A+T+Asyn+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnief_ad,Dummy_AD);
% Plot
boundedline(Dummy_AD.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','g')

% Set Mixed
Dummy_AD(:,'Group') = {'Asyn+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnief_ad,Dummy_AD);
% Plot
boundedline(Dummy_AD.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','b')

ylabel('ADNI EF')
xlabel('Time from baseline (y)')

subplot(2,4,8);

% Set AD-like
Dummy_AD(:,'Group') = {'A+T+Asyn-'};
% Predict
[ypred,yCI,DF] = predict(lme_adnivs_ad,Dummy_AD);
% Plot
boundedline(Dummy_AD.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','r')

% Set Mixed
Dummy_AD(:,'Group') = {'A+T+Asyn+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnivs_ad,Dummy_AD);
% Plot
boundedline(Dummy_AD.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','g')

% Set LB
Dummy_AD(:,'Group') = {'Asyn+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnivs_ad,Dummy_AD);
% Plot
boundedline(Dummy_AD.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','b')

ylabel('D(MEM-EF)')
xlabel('Time from baseline (y)')


