clc
clear 
close all
addpath('pca_ica/pca_ica')
plotData = false;

% GUIforInput();
% 
% return
%% Program to unfold the data and pick the ERP at the given event

X = load('Subject05_s2.mat');
[stackedERPdata, stackedERPidx, stackedNotERPdata, stackedNotERPidx] = exploreOpen(X,6);
erpDataCombined = [];
NoterpDataCombined = [];
for i=1:10
    erpDataCombined = [erpDataCombined; stackedERPdata{i}];
    NoterpDataCombined = [NoterpDataCombined; stackedNotERPdata{i}];
end


target = [ones(1,size(erpDataCombined,1)) zeros(1,size(NoterpDataCombined,1));...
          zeros(1,size(erpDataCombined,1)) ones(1,size(NoterpDataCombined,1))];


if(plotData)
    dat = (X.run{1, 1}.eeg)';
    pos = X.run{1, 1}.header.EVENT.POS;  
    typ = X.run{1, 1}.header.EVENT.TYP;        
    erpDataStPoint = pos(find(typ==5));
    erpDataEndPoint = erpDataStPoint+512;
    temp = min(dat(1,:));
    x = [erpDataStPoint erpDataStPoint]';
    y = [zeros(1,length(erpDataStPoint)); temp*ones(1,length(erpDataStPoint))];
    x1 = [erpDataEndPoint erpDataEndPoint]';
    plot(dat(1,:));
    hold on
    line(x,y,'Color','red');
    hold on
    line(x1,y,'Color','green');
    %ylim([max(dat(1,:)) min(dat(1,:))])
    hold off
end

%% Analysis for ERP data
Z1 = erpDataCombined';
[Zw1, T] = whitenRows(Z1); % Whitening rows
r = 64;
[Zpca1, U, mu, eigVecs] = PCA(Zw1,r); %PCA
[Zica1, W, T, mu] = kICA(Zpca1,r); %ICA

%% Analysis for Non ERP data
Z2 = NoterpDataCombined';
[Zw2, T] = whitenRows(Z2); % Whitening rows
r = 64;
[Zpca2, U, mu, eigVecs] = PCA(Zw2,r); %PCA
[Zica2, W, T, mu] = kICA(Zpca2,r); %ICA

DataTrain = [Zica1 Zica2];

%% Testing 



function [stackedERPdata, stackedERPidx, stackedNotERPdata, stackedNotERPidx] = exploreOpen(X,eveTyp)   
    X = X.run;
    l = length(X);    
    for i=1:l
        stackERPidx = [];
        dat = X{1,i};
        eeg = dat.eeg;
        stackNotERPidx = 1:size(eeg,1);
        header = dat.header;
        SampleRate = header.SampleRate;
        %label = header.Label;
        eventType = header.EVENT.TYP;
        eventPos = header.EVENT.POS;       
        erpDataStPoint = eventPos(find(eventType==eveTyp));
        numERP = length(erpDataStPoint);
        for k=1:numERP            
            stackERPidx = [stackERPidx erpDataStPoint(k):erpDataStPoint(k)+SampleRate-1];           
        end
        stackNotERPidx(stackERPidx) = [];
        eegERP = eeg(stackERPidx,:);
        eegNotERP = eeg(stackNotERPidx,:);
        stackedERPdata{i} = eegERP;
        stackedERPidx{i} = stackERPidx;        
        stackedNotERPdata{i} = eegNotERP;
        stackedNotERPidx{i} = stackNotERPidx;       
    end
end


