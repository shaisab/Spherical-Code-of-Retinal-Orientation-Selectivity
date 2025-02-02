
mat=mat180;

figure;
resolution=5;
alpha=0:resolution:360;
beta=0:resolution:180;
contourf(alpha,beta,mat);


matS=smooth2a(mat,3,3);

figure
contourf(alpha,beta,matS);

