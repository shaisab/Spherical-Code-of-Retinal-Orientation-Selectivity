function []= Shufflepori()
% figure;

filedir=cd;
%% Load files
[FileName,PathName] = uigetfile('*.mat','Select pooledmap file for all clusters');

load(FileName);


Shufflepori.dir=[0,45,90,135,180,225,270,315]; %stimulus direction

%ON SUS - 5, ON TRANS - 2, ON-OFF SUS - 3, ON-OFF TRANS - 8
            
for i=1:67 %choosing half the cells only from the first retina (cells were duplicated)
    data=pooledmap.meanr_all{i,:}; % mean response with reps
    d=Shufflepori.dir;
    
    for j=1:100
    ind=randsample(4,4,'true'); %random indexing of trial reps (4 reps)
    data_shuffled=[data(ind(1),:);data(ind(2),:);data(ind(3),:);data(ind(4),:)]; % mean response with shuffled reps
    a=mean(data_shuffled,1); %mean response for all shuffled reps 
    
    [oriFWHM,u,k,gof2,x,xq,Rq]=doubleGaussianFitOSI(j,d,a); %doubleGaussianFit for mean shuffled response
    
    Shufflepori.pori(i).pori_Shuffle(j)=u;
    Shufflepori.pori(i).R2(j)=gof2.rsquare;

    end
    
    a= mean(pooledmap.meanr_all{i,:}); %mean response for all original reps
    [oriFWHM,u,k,gof2,x,xq,Rq]=doubleGaussianFitOSI(j,d,a); %doubleGaussianFit for mean response for all original reps
    Shufflepori.pori_real(i)=u;
    Shufflepori.R2_real(i)=gof2.rsquare;
    Shufflepori.OSI(i)=pooledmap.OSI(i);
    Shufflepori.meanr_all(i)=pooledmap.meanr_all(i);
    Shufflepori.typeClust(i)=pooledmap.typeClust(i);
    
end
 

    ind=strfind(filedir,'\');
    filename=filedir;
    filename(ind(end):length(filename))=[];
    filename(1:ind(end-1))=[];
    save(['Shufflepori_',filename,'.mat'],'Shufflepori');
    
end
