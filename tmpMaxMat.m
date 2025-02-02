
resolution=10;


%%% code for generating matchind matrices and plotting contour plots
c=10;
cutoff=(['c',num2str(c)]);
k=1;
h=1;
for alpha=0:resolution:360
for beta=0:resolution:180
    mat180(k,h)=allPossibleFields.matchind180(k,h).(cutoff);
    mat90(k,h)=allPossibleFields.matchind90(k,h).(cutoff);
    k=k+1;
end
k=1;
h=h+1;
end


%%% find max singularities
for i=1:4
    if i==1
        [~,I]=max(mat180(:));
        [Irow,Icol]=ind2sub(size(mat180),I);
    else
        Icol=Icolmax(i-1)+(90/resolution);
        if Icol>size(mat180,2)
            Icol=Icol-size(mat180,2);
        end
        [~,Irow]=max(mat180(:,Icol)) ;    
    end
Icolmax(i)=Icol;
Irowmax(i)=Irow;
end
alpha=0:resolution:360;
beta=0:resolution:180;
alphasing(i)=alpha(Icol);
betasing(i)=beta(Irow);
