


%% generate matchind matrix
resolution=10;
c=10;
cutoff=(['c',num2str(c)]);
k=1;
h=1;
for alpha=0:resolution:360
for beta=0:resolution:180
    mat180(k,h)=allPossibleFields.matchind180(k,h).(cutoff);    
    k=k+1;
end
k=1;
h=h+1;
end
tmp=reshape(mat180,[],1);
figure;
histogram(tmp,30)


y=quantile(tmp,[0.25 0.75]) ;
bandwidth=y(2)-y(1)

%median(mat180(:))
% mean(mat180(:))
% std(mat180(:))