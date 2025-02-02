% clear all
% close all

load('pooledMapAlphaCorr_OS data_40_OSI0.20_1.00_SI_0.00_created_Dec_01_2022_14_56.mat');

% x=abs(asind(sind(allCells.poriSz)));
% Orientation_alignment_index=sum(x>60|x<30)/length(x);

for i=1:size(pooledmap.cellID,1)
    for j=1:size(pooledmap.clusterdata,2)
        
        try Orientation_alignment=abs(asind(sind(pooledmap.clusterdata{1, j}.poriSz)));
        pooledmap.clusterdata{1, j}.Orientation_alignment_index=sum(Orientation_alignment>60|Orientation_alignment<30)/length(Orientation_alignment);
        
        catch pooledmap.clusterdata{1, j}.Orientation_alignment_index=0;
            
        end
    end 
end

% pori = asind(pooledmap.UVvecComp.V);

% for i=1:size(pooledmap.cellID,1)
%     for j=1:size(pooledmap.clusterdata,2)
%         for k=1:size(pooledmap.clusterdata{1, j}.cellID)
%             index=find(strcmp(pooledmap.cellID, pooledmap.clusterdata{1, j}.cellID{k, 1} ),1);
%             pori=abs(asind(pooledmap.UVvecComp.V(index)));
%             pooledmap.clusterdata{1, j}.pori(k, 1)=pori;
%         end
%         try pooledmap.clusterdata{1, j}.Orientation_alignment_index=sum(pooledmap.clusterdata{1, j}.pori<5)/sum(pooledmap.clusterdata{1, j}.pori>0);
%         
%         catch pooledmap.clusterdata{1, j}.Orientation_alignment_index=0;
%             
%         end
%     end 
% end