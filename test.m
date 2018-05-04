clc
clear
close all

% load('dataERP6.mat');
% load('dataNotERP6.mat');

X = load('Subject05_s1.mat');
[stackedERPdata, stackedERPidx, stackedNotERPdata, stackedNotERPidx] = exploreOpen(X,6);
temp = stackedNotERPdata{1,1};
% x1 = dataERP(1:512,:)';

for i=1:1
    
    x2 = temp(1:512,:)';
    
    Z1 = x2;
    [Zw1, T] = whitenRows(Z1); % Whitening rows
    r = 64;
    [Zpca1, U, mu, eigVecs] = PCA(Zw1,r); %PCA
    [Zica1, W, T, mu] = kICA(Zpca1,r); %ICA
    
    y1 = myNeuralNetworkFunction(Zica1);
    sum(y1')
end
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
