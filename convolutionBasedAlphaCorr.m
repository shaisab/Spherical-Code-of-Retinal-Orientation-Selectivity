
function []= convolutionBasedAlphaCorr(resolution)

%% calculate the average alpha-beta map to be used as a reference for the 1-D convolution
% load all allPossibleFields files
[FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files to avarege','MultiSelect','on');

if ~iscell(FileName)
    alphaBetaMap{1}=load(FileName);
else
    for i=1:size(FileName,2)    
        alphaBetaMap{i}=load(FileName{i});
    end
end

% extract alpha beta values for the given cutoff        
c=10;
cutoff=(['c',num2str(c)]);
k=1;
h=1;
        for t=1:size(alphaBetaMap,2)
            k=1;
            h=1;
            for alpha=0:resolution:360
                for beta=0:resolution:180
                    mat180{t}(k,h)=alphaBetaMap{1,t}.allPossibleFields.matchind180(k,h).(cutoff);
                    k=k+1;
                end
                k=1;
                h=h+1;
            end
        end
% calculate the mean matrix
b=cellfun(@(x) [x],mat180,'un',0)
meanMat180=mean(cat(3,b{:}),3)

figure;
alpha=0:resolution:360;
beta=0:resolution:180;
contourf(alpha,beta,meanMat180);


%% calculate alpha correction for all retinas
filedir=cd;
namesMap=getFileList(filedir,'allPossibleFields_',0,'anywhere');

for i=1:size(namesMap,2)
    ind=[];
    Map{1,i}=load(namesMap{i});
        
% extract alpha beta values for the given cutoff        
k=1;
h=1;               
for alpha=0:resolution:360
    for beta=0:resolution:180
        ABmap{1,i}(k,h)=Map{1,i}.allPossibleFields.matchind180(k,h).(cutoff);
        k=k+1;
    end
    k=1;
    h=h+1;
end
                 
figure;
subplot(1,2,1);
alpha=0:resolution:360;
beta=0:resolution:180;
for shift=1:size(meanMat180,2)
d=circshift(ABmap{1,i},[0 shift]);
shiftProduct{shift,i}=d.*meanMat180;
shiftSum(shift,i)=sum(sum(shiftProduct{shift,i}));
contourf(alpha,beta,shiftProduct{shift,i});
drawnow expose
end

subplot(1,2,2);
plot(shiftSum(:,i));
hold off
[~,shiftMax(i)]=max(shiftSum(:,i));

shiftMaxS=shiftMax-size(meanMat180,2);
 for s=1:size(shiftMax,2)
     if shiftMaxS(s)<-floor(size(meanMat180,2)/2)
     alphaCorr(s)=resolution*(shiftMaxS(s)+size(meanMat180,2));
 elseif shiftMaxS(s)>=-floor(size(meanMat180,2)/2)    
     alphaCorr(s)=resolution*shiftMaxS(s);
     end
 end

ind=strfind(namesMap{i},'_');
filename{i}=namesMap{i}(ind(1)+1:ind(4)-1);

end

figure;
scatter(1:size(alphaCorr,2),alphaCorr)

% save file alphaCorr.mat
ind=strfind(namesMap{1},'_');
inddot=strfind(namesMap{1},'.');
try
    type=namesMap{1}(ind(14)+1:inddot(1)-1);
catch
    type='';
end
currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');

save(['alphaCorr_',num2str(type),'_',currenttime,'.mat'],'filename','shiftProduct','shiftSum','shiftMax','alphaCorr');

end


