

a1=dlmread('Sec1VecField_50_ASC');
U1=a1(:,1:50);
V1=a1(:,51:100);
vecU1=a1(:,101:150);
vecV1=a1(:,151:200);
a2=dlmread('Sec2VecField_50_ASC');
U2=a2(:,1:50);
V2=a2(:,51:100);
vecU2=a2(:,101:150);
vecV2=a2(:,151:200);
a3=dlmread('Sec3VecField_50_ASC');
U3=a3(:,1:50);
V3=a3(:,51:100);
vecU3=a3(:,101:150);
vecV3=a3(:,151:200);
a4=dlmread('Sec4VecField_50_ASC');
U4=a4(:,1:50);
V4=a4(:,51:100);
vecU4=a4(:,101:150);
vecV4=a4(:,151:200);



figure
for i=1:50
plot(vecU1(:,i),vecV1(:,i))
hold on;
end

for i=1:50
plot(vecU2(:,i),vecV2(:,i))
hold on;
end

for i=1:50
plot(vecU3(i,:),vecV3(i,:))
hold on;
end

for i=1:2:50
plot(vecU4(:,i),vecV4(:,i))
hold on;
end


U=(1:10)';
V=(10:20)';
vecU=(20:30)';
vecV=(30:40)';

X=[];
Y=[];
    for i=1:50
        X=[X;[U1(i,j);vecU1(i,j)]];
        Y=[Y;[V1(i,j);vecV1(i,j)]];
    end


figure;
plot(X,Y);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

