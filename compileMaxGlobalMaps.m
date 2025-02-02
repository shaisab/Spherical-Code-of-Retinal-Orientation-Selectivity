function []=compileMaxGlobalMaps()


%% Load files and extract meanOutput matrices
filedir=cd;
names=getFileList(filedir,'.mat',0,'anywhere');     % detect all .mat files
Right=[];
type={'translation','rotation'};
for t=1:2
    typestr=type{t};
    for i=1:size(names,2)
        if strfind(names{i},typestr)
            if strfind(names{i},'Right')
                if strfind(names{i},'sing1')
                    load(num2str(names{i}));
                    Right.(typestr){1,1}=[allPossibleFields.maxAziR,allPossibleFields.maxElevR];
                elseif strfind(names{i},'sing2')
                    load(num2str(names{i}));
                    Right.(typestr){1,2}=[allPossibleFields.maxAziR,allPossibleFields.maxElevR];
                elseif strfind(names{i},'sing3')
                    load(num2str(names{i}));
                    Right.(typestr){1,3}=[allPossibleFields.maxAziR,allPossibleFields.maxElevR];
                elseif strfind(names{i},'sing4')
                    load(num2str(names{i}));
                    Right.(typestr){1,4}=[allPossibleFields.maxAziR,allPossibleFields.maxElevR];
                end
            elseif strfind(names{i},'Left')
                if strfind(names{i},'sing1')
                    load(num2str(names{i}));
                    Left.(typestr){1,1}=[allPossibleFields.maxAziL,allPossibleFields.maxElevL];
                elseif strfind(names{i},'sing2')
                    load(num2str(names{i}));
                    Left.(typestr){1,2}=[allPossibleFields.maxAziL,allPossibleFields.maxElevL];
                elseif strfind(names{i},'sing3')
                    load(num2str(names{i}));
                    Left.(typestr){1,3}=[allPossibleFields.maxAziL,allPossibleFields.maxElevL];
                elseif strfind(names{i},'sing4')
                    load(num2str(names{i}));
                    Left.(typestr){1,4}=[allPossibleFields.maxAziL,allPossibleFields.maxElevL];
                end
            end
        end
    end
    
end
save('singularities of global maps','Left','Right');

end
