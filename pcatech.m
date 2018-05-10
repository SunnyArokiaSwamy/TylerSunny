clc
clear
close all

addpath('/Users/sunnyarokiaswamybellary/Documents/EEGLAB_Code/ERP_Dataset/Datasets');
addpath('pca_ica')
electrode = 1;

% Training the Network
Subjects = [1];
Events = [1];
EveNum = [6 9];
[outputDataERP,outputDataNotERP]= stackERP_NotERP(Subjects, Events, EveNum);
outputDataNotERP = outputDataNotERP(:,1:size(outputDataERP,2));

erpTemp = reshape(outputDataERP(electrode,:),[512 size(outputDataERP,2)/512])';
NotTemp = reshape(outputDataNotERP(electrode,:),[512 size(outputDataNotERP,2)/512])';
Data = [erpTemp' NotTemp'];
targets1 = [ones(1,size(erpTemp,1)) zeros(1,size(NotTemp,1)); zeros(1,size(erpTemp,1)) ones(1,size(NotTemp,1))];

% Analysis for ERP data
Z1 = Data;
[Zw1, T] = whitenRows(Z1); % Whitening rows
r = 64;
[Zpca1, U, mu, eigVecs] = PCA(Zw1,r); %PCA
[Zica1, W, T, mu] = kICA(Zpca1,r); %ICA


%% Testing
% Training the Network
Subjects = [1];
Events = [2];
EveNum = [6 9];
[outputDataERP,outputDataNotERP]= stackERP_NotERP(Subjects, Events, EveNum);
outputDataNotERP = outputDataNotERP(:,1:size(outputDataERP,2));

erpTemp = reshape(outputDataERP(electrode,:),[512 size(outputDataERP,2)/512])';
NotTemp = reshape(outputDataNotERP(electrode,:),[512 size(outputDataNotERP,2)/512])';
Data = [erpTemp' NotTemp'];
targets2 = [ones(1,size(erpTemp,1)) zeros(1,size(NotTemp,1)); zeros(1,size(erpTemp,1)) ones(1,size(NotTemp,1))];
Z1 = Data;
[Zw1, T] = whitenRows(Z1); % Whitening rows
r = 64;
[Zpca1, U, mu, eigVecs] = PCA(Zw1,r); %PCA
[Zica2, W, T, mu] = kICA(Zpca1,r); %ICA



% % Analysis for Non ERP data
% Z2 = NoterpDataCombined';
% [Zw2, T] = whitenRows(Z2); % Whitening rows
% r = 64;
% [Zpca2, U, mu, eigVecs] = PCA(Zw2,r); %PCA
% [Zica2, W, T, mu] = kICA(Zpca2,r); %ICA
% 
% DataTrain = [Zica1 Zica2];

function [outputDataERP,outputDataNotERP]= stackERP_NotERP(Subjects, Events, EveNum)
% Stack as ERP and Not ERP
stk=[];
tk=[];
outputDataERP = [];
outputDataNotERP = [];
for i=1:length(Subjects)
    subNum = Subjects(i);
    for j=1:length(Events)
        eventNum = Events(j);
        str = strcat('Subject0',num2str(subNum),'_s',num2str(eventNum));
        X = load(str);
        for eve=1:length(EveNum)
            [stackedERPdata, stackedERPidx, stackedNotERPdata, stackedNotERPidx] = exploreOpen(X,EveNum(eve));
            stackedERPdata = cat(1,stackedERPdata{:})';
            stackedNotERPdata = cat(1,stackedNotERPdata{:})';
            outputDataERP = [outputDataERP; stackedERPdata'];
            outputDataNotERP = [outputDataNotERP; stackedNotERPdata'];
        end     
    end
end

outputDataERP = outputDataERP';
outputDataNotERP = outputDataNotERP';

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

