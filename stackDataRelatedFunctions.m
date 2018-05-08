clc
clear
close all

addpath('W:\Documents\MATLAB\Thesis\Database');
Subjects = [1 2 3 4 5 6];
Events = [1];
EveNum = [6 9];


% Concatenates all the data
if(0)
    stk=[];
    tk=[];
    for i=1:length(Subjects)
    subNum = Subjects(i);
    for j=1:length(Events)
        eventNum = Events(j);
        str = strcat('Subject0',num2str(subNum),'_s',num2str(eventNum));
        load(str);
        totalTimesPerSub = length(run); % Can be controlled by subtracting
        for k=1:totalTimesPerSub
           temp = run{1,k}.eeg';
           tk = [tk; size(run{1,k}.eeg')];
           stk = [stk temp]; 
        end        
    end
end
end

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





