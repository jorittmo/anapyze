figure();

formula_adnimem_mmse = 'ADNI_MEM ~ TimeFromBaseliney*(Group) + BL_Age + Sex + Education + MMSE + (TimeFromBaseliney|Subject)';
formula_adnief_mmse = 'ADNI_EF ~ TimeFromBaseliney*(Group) +  BL_Age + Sex + Education + MMSE + (TimeFromBaseliney|Subject)';
formula_adnidelta_mmse = 'DMEMEF ~ TimeFromBaseliney*(Group) +  BL_Age + Sex + Education + MMSE + (TimeFromBaseliney|Subject)';
formula_mmse = 'MMSE ~ TimeFromBaseliney*(Group) +  BL_Age + Sex + Education + (TimeFromBaseliney|Subject)';

%MMSE
lme_mmse = fitlme(CLINICALMMSELongit,formula_mmse);

%ADNI_MEM
lme_adnimem = fitlme(CLINICALCOMPLongit,formula_adnimem_mmse);

%ADNI_EF
lme_adnief = fitlme(CLINICALCOMPLongit,formula_adnief_mmse);

%ADNI_DMEMEF
lme_adnidelta = fitlme(CLINICALCOMPLongit,formula_adnidelta_mmse);

%-----------------------------------------------------------------------

% Predict and Plot

subplot(2,2,1);

% MMSE MCI

% Set AD+LB-
Dummy_Patient(:,'Group') = {'AD+LB-'};
% Predict
[ypred,yCI,DF] = predict(lme_mmse,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','r')

% Set AD-LB+
Dummy_Patient(:,'Group') = {'AD-LB+'};
% Predict
[ypred,yCI,DF] = predict(lme_mmse_mci,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','g')

% Set AD+LB+
Dummy_Patient(:,'Group') = {'AD+LB+'};
% Predict
[ypred,yCI,DF] = predict(lme_mmse_mci,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','b')

% Set AD-LB-
Dummy_Patient(:,'Group') = {'AD-LB-'};
% Predict
[ypred,yCI,DF] = predict(lme_mmse_mci,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','m')


ylabel('MMSE')
xlabel('Time from baseline (y)')

subplot(2,2,2);

% ADNI_MEM

% Set AD+LB-
Dummy_Patient(:,'Group') = {'AD+LB-'};
% Predict
[ypred,yCI,DF] = predict(lme_adnimem,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','r')

% Set AD-LB+
Dummy_Patient(:,'Group') = {'AD-LB+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnimem,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','g')

% Set AD+LB+
Dummy_Patient(:,'Group') = {'AD+LB+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnimem,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','b')

% Set AD-LB-
Dummy_Patient(:,'Group') = {'AD-LB-'};
% Predict
[ypred,yCI,DF] = predict(lme_adnimem,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','m')

ylabel('ADNI MEM')
xlabel('Time from baseline (y)')

subplot(2,2,3);

% ADNI_EF

% Set AD+LB-
Dummy_Patient(:,'Group') = {'AD+LB-'};
% Predict
[ypred,yCI,DF] = predict(lme_adnief,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','r')

% Set AD-LB+
Dummy_Patient(:,'Group') = {'AD-LB+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnief,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','g')

% Set AD+LB+
Dummy_Patient(:,'Group') = {'AD+LB+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnief,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','b')

% Set AD-LB-
Dummy_Patient(:,'Group') = {'AD-LB-'};
% Predict
[ypred,yCI,DF] = predict(lme_adnief,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','m')

ylabel('ADNI EF')
xlabel('Time from baseline (y)')

subplot(2,2,4);

% D(MEM_EF)

% Set AD+LB-
Dummy_Patient(:,'Group') = {'AD+LB-'};
% Predict
[ypred,yCI,DF] = predict(lme_adnivs_mci,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','r')

% Set AD-LB+
Dummy_Patient(:,'Group') = {'AD-LB+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnivs_mci,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','g')

% Set AD+LB+
Dummy_Patient(:,'Group') = {'AD+LB+'};
% Predict
[ypred,yCI,DF] = predict(lme_adnivs_mci,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','b')

% Set AD-LB-
Dummy_Patient(:,'Group') = {'AD-LB-'};
% Predict
[ypred,yCI,DF] = predict(lme_adnivs_mci,Dummy_Patient);
% Plot
boundedline(Dummy_Patient.TimeFromBaseliney, ypred, [ypred-yCI(:,1) yCI(:,2)-ypred],'alpha','m')

ylabel('D(MEM-EF)')
xlabel('Time from baseline (y)')

