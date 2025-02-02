
function []= combineAllpossibleMaps()


%% load and combine all PossibleFields files
clear all;

resolution=10;

[FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of the four single-channel maps','MultiSelect','on');
        
        % sort file order by ascending alpha
        for k=1:size(FileName,2)
        indUnder=find(FileName{k}=='_');
        order(k)=str2num(FileName{k}(indUnder(4)+1:indUnder(5)-1));
        end
        [~,I]=sort(order);
        FileNameSorted=FileName(I);
        
        for i=1:size(FileName,2)
            transField{i}=load(FileNameSorted{i});
        end
allOutput={transField{1,1}.allPossibleFields.meanOutput,transField{1,2}.allPossibleFields.meanOutput,...
    transField{1,3}.allPossibleFields.meanOutput,transField{1,4}.allPossibleFields.meanOutput};     

figure;

%% plot selection over four translational channels
clear tmp

%%%%%%%%%%%%%%%%%%
% input parameters
alpha=':';
beta=20;
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%

if alpha==':'
    a=':';
    b=(beta/resolution)+1;
x=0:10:360;
tmp=[allOutput{1,1}(b,a);allOutput{1,2}(b,a);allOutput{1,3}(b,a);allOutput{1,4}(b,a)];
plot(x,tmp)
xlim([0 360]);
ylim([0 1]);
elseif beta==':'
    b=':';
    a=(alpha/resolution)+1;
x=0:10:180;
tmp=[allOutput{1,1}(b,a),allOutput{1,2}(b,a),allOutput{1,3}(b,a),allOutput{1,4}(b,a)];
plot(x,tmp)
xlim([0 180]);
ylim([0 1]);
end

%% calculate pair-wise euclidian distance between four translation channels

X=[reshape(allOutput{1,1},1,[]);reshape(allOutput{1,2},1,[]);reshape(allOutput{1,3},1,[]);reshape(allOutput{1,4},1,[])];
for j=1:size(X,2)
eDistance(j,:)=pdist(X(:,j),'euclidean');
end
meaneDistance=mean(eDistance,2);
meaneDistance=reshape(meaneDistance,size(allOutput{1,1},1),[]);

%% plot meanOutput map
f1=figure;
alpha=0:resolution:360;
beta=0:resolution:180;
contourf(alpha,beta,meaneDistance);

end