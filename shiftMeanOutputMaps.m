
function []=shiftMeanOutputMaps(elements)

        resolution=10;

        [FileName,PathName] = uigetfile('*.*','Select a allPossibleFields file');

        if strfind(FileName,'allPossibleFields')
            load(FileName);
           data=allPossibleFields.meanOutput; 
        end
    

% [~,I]=max(mat180tmp(:));
%         [Irow,Icol]=ind2sub(size(mat180tmp),I);        
        

f1=figure;
alpha=0:resolution:360;
beta=0:resolution:180;
contourf(alpha,beta,data);


d=circshift(data,[0 elements]);

f2=figure;
alpha=0:resolution:360;
beta=0:resolution:180;
contourf(alpha,beta,d);

end


