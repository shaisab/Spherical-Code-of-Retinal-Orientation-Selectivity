
function mat180 = generateMatchindMatrixMix(allPossibleFields)

resolution=10;
c=10;
cutoff=(['c',num2str(c)]);
k=1;
h=1;
for alpha=120:resolution:300
for beta=0:resolution:90
    mat180(k,h)=allPossibleFields.matchindMix(k,h).(cutoff);
    k=k+1;
end
k=1;
h=h+1;
end

end