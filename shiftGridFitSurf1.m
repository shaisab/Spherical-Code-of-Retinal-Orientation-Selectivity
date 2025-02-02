
function []=shiftGridFitSurf1(elements)

        resolution=10;

        [FileName,PathName] = uigetfile('*.*','Select a allPossibleFields file');

        if strfind(FileName,'allPossibleFields')
            load(FileName);
            
            c=10;
            cutoff=(['c',num2str(c)]);
            k=1;
            h=1;
            for Alpha=0:resolution:360
                for beta=0:resolution:180
                    mat180(k,h)=allPossibleFields.matchind180(k,h).(cutoff);
                    mat90(k,h)=allPossibleFields.matchind90(k,h).(cutoff);
                    k=k+1;
                end
                k=1;
                h=h+1;
            end   
            data=mat180;
        end
    


f1=figure;
alpha=0:resolution:360;
beta=0:resolution:180;
contourf(alpha,beta,data);




% a=1:20;
% b=repmat(a,10,1);
d=circshift(data,[0 elements]);

f2=figure;
alpha=0:resolution:360;
beta=0:resolution:180;
contourf(alpha,beta,d);

end


