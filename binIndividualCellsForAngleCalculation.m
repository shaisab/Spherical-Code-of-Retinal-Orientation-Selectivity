

j=1;
for i=1:size(pooledmap.UVStvecComp.USt,1)
    ind=intersect(find(pooledmap.UVStvecComp.USt==pooledmap.UVStvecComp.USt(i)),find(pooledmap.UVStvecComp.VSt==pooledmap.UVStvecComp.VSt(i)));
    if ~isempty(ind)
        for k=1:size(ind,1)
            if pooledmap.Itrans(ind(k))==1
                PU(j,1)=pooledmap.UVStvecComp.USt(ind(k));
                PV(j,1)=pooledmap.UVStvecComp.VSt(ind(k));
                PvecU(j,1)=pooledmap.UVStvecComp.VecUSt(ind(k));
                PvecV(j,1)=pooledmap.UVStvecComp.VecVSt(ind(k));
            elseif pooledmap.Itrans(ind(k))==2
                PU(j,2)=pooledmap.UVStvecComp.USt(ind(k));
                PV(j,2)=pooledmap.UVStvecComp.VSt(ind(k));
                PvecU(j,2)=pooledmap.UVStvecComp.VecUSt(ind(k));
                PvecV(j,2)=pooledmap.UVStvecComp.VecVSt(ind(k));
            elseif pooledmap.Itrans(ind(k))==3
                PU(j,3)=pooledmap.UVStvecComp.USt(ind(k));
                PV(j,3)=pooledmap.UVStvecComp.VSt(ind(k));
                PvecU(j,3)=pooledmap.UVStvecComp.VecUSt(ind(k));
                PvecV(j,3)=pooledmap.UVStvecComp.VecVSt(ind(k));
            elseif pooledmap.Itrans(ind(k))==4
                PU(j,4)=pooledmap.UVStvecComp.USt(ind(k));
                PV(j,4)=pooledmap.UVStvecComp.VSt(ind(k));
                PvecU(j,4)=pooledmap.UVStvecComp.VecUSt(ind(k));
                PvecV(j,4)=pooledmap.UVStvecComp.VecVSt(ind(k));  
            end 
        end
        j=j+1;
    end
end