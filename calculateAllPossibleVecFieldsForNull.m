function calculateAllPossibleVecFieldsForNull(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,numIter,polarity2,mapType)

load(FileName);     % loads pooledmap file

polarity=1;
discPoints=40;
vizAngle=180;
interpType=2;       % find closest grid point

%% loads single-channel maps
switch mapType
    case 0      % translation - singularities based on real data
        
        [FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of single-channel translation maps','MultiSelect','on');
        
        % sort file order by ascending alpha
        for k=1:size(FileName,2)
            indUnder=find(FileName{k}=='_');
            order(k)=str2num(FileName{k}(indUnder(4)+1:indUnder(5)-1));
        end
        [~,I]=sort(order);
        FileNameSorted=FileName(I);
        
        for i=1:size(FileName,2)
            load(FileNameSorted{i});
            transField{i}=generateMatchindMatrix(allPossibleFields);
        end
        
        [FileName,PathName] = uigetfile('*.mat','Select allPossibleFields translation file of the real data','MultiSelect','on');
        if strfind(FileName,'_real')
            load(FileName);
            real=generateMatchindMatrix(allPossibleFields);
            alphasing=allPossibleFields.alphasing;
            betasing=allPossibleFields.betasing;
        end
        
        
        case 2      % mixed orthogonal
        
        [FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of single-channel mixed orthogonal maps','MultiSelect','on');
        
            load(FileName);
            transField=generateMatchindMatrixMix(allPossibleFields);
        
        
        [FileName,PathName] = uigetfile('*.mat','Select allPossibleFields mixed orthogonal file of the real data','MultiSelect','on');
        if strfind(FileName,'_real')
            load(FileName);
            real=generateMatchindMatrix(allPossibleFields);
            alphasing=allPossibleFields.alphasing;
            betasing=allPossibleFields.betasing;
        end
        
        
    case 1  % rotation
        
switch polarity2
             case 1
                [FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of the three single-channel rotation maps and the real data','MultiSelect','on');
                for i=1:size(FileName,2)
                    if strfind(FileName{i},'_ASC')
                        load(FileName{i});
                        matASC=generateMatchindMatrix(allPossibleFields);
                    elseif strfind(FileName{i},'_PSC')
                        load(FileName{i});
                        matPSC=generateMatchindMatrix(allPossibleFields);
                    elseif strfind(FileName{i},'_LSC')
                        load(FileName{i});
                        matLSC=generateMatchindMatrix(allPossibleFields);
                    else
                        load(FileName{i});
                        real=generateMatchindMatrix(allPossibleFields);
                    end
                end
                
            case 2  
                [FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of the six single-channel rotation maps and the real data','MultiSelect','on');
                for i=1:size(FileName,2)
                    if strfind(FileName{i},'_ASC')
                        load(FileName{i});
                        matASC=generateMatchindMatrix(allPossibleFields);
                    elseif strfind(FileName{i},'_PSC')
                        load(FileName{i});
                        matPSC=generateMatchindMatrix(allPossibleFields);
                    elseif strfind(FileName{i},'_LSC')
                        load(FileName{i});
                        matLSC=generateMatchindMatrix(allPossibleFields);
                    elseif strfind(FileName{i},'_iASC')
                        load(FileName{i});
                        matiASC=generateMatchindMatrix(allPossibleFields);
                    elseif strfind(FileName{i},'_iPSC')
                        load(FileName{i});
                        matiPSC=generateMatchindMatrix(allPossibleFields);
                    elseif strfind(FileName{i},'_iLSC')
                        load(FileName{i});
                        matiLSC=generateMatchindMatrix(allPossibleFields);
                    elseif strfind(FileName{i},'_real')
                        load(FileName{i});
                        real=generateMatchindMatrix(allPossibleFields);
                    end
                end     
end
end
%% generate alpha-beta maps for random vectors

if DATAMODELstate
    data=pooledmap.clusterdata{1,CLUSTNstate};
elseif ~DATAMODELstate
    data=pooledmap;
end


for i=1:numIter
    clear allPossibleFields
    allPossibleFields.field={};
    k=1;
    h=1;
    maxVec=[];  % for compatability with VecComparisonPooled4
    DISTstate=[];  % for compatability with VecComparisonPooled4
    locDist=[];  % for compatability with VecComparisonPooled4
    widthDist=[];  % for compatability with VecComparisonPooled4
    
    a=deg2rad(360*rand(size(data.XYZvecComp.VecXcompLSC,1),1));
    a=a';
    for alpha=120:resolution:300
        for beta=0:resolution:90
            if mapType==2
                fieldType=['alpha_',num2str(alpha),'_','beta_',num2str(beta),'_','iteration_',num2str(i)];
                allPossibleFields=VecComparisonPooled4Trans(data,a,maxVec,deg2rad(alpha),deg2rad(beta),polarity,fieldType,vizAngle,allPossibleFields,k,h,interpType,histType);     
                k=k+1;
            else
                fieldType=['alpha_',num2str(alpha),'_','beta_',num2str(beta),'_','iteration_',num2str(i)];
                allPossibleFields=VecComparisonPooled4(data,a,maxVec,deg2rad(alpha),deg2rad(beta),polarity,fieldType,vizAngle,allPossibleFields,k,h,interpType,histType,DISTstate,locDist,widthDist);
                k=k+1;
            end
        end
        k=1;
        h=h+1;
        
    end
    
    
    % currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
    % switch ALPHACORRstate
    %     case 1
    %         switch histType
    %             case 10
    %                     save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_rand_',num2str(i),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
    %         end
    %
    %     case 0
    %         switch histType
    %             case 10
    %                     save(['allPossibleFields_',num2str(filename),'_rand_',num2str(i),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
    %         end
    % end
    
    
    
    %%% code for generating matchind matrices
    c=10;
    cutoff=(['c',num2str(c)]);
    k=1;
    h=1;
    for alpha=120:resolution:300
        for beta=0:resolution:90
            if mapType==2
                null.randmat(i).mat180(k,h)=allPossibleFields.matchindMix(k,h).(cutoff);
            else
                null.randmat(i).mat180(k,h)=allPossibleFields.matchind180(k,h).(cutoff);
            end
            k=k+1;
        end
        k=1;
        h=h+1;
    end
    
    matrand= null.randmat(i).mat180;      % extract the current map to be used for the fit below
    
    %% process mat180 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch mapType
        case 0  % translation
            
            C1=reshape(transField{1},[],1);
            C2=reshape(transField{2},[],1);
            Cr=reshape(matrand,[],1);
            
            %%% fit four translation channels
            % non linear optimization
            guess=[0.5,0.5];
            lowerBound=[0,0];
            upperBound=[1,1];
            regressor=[C1,C2];
            
            options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
            [x,resnorm,residual,exitflag,output]=lsqnonlin(@(guess)myfunSurfDynamic(guess,regressor,Cr),guess,lowerBound,upperBound,options)
            
            % multiple linear regression
            [b,bint,r,rint,stats]=regress(Cr,[regressor]);
            
            % multiple regrression without a constant from my dimentionality of color
            % vision work based on Maloney's paper.
            X=regressor;
            Y=Cr;
            B=X\Y;          % this the same as: B=inv(X'*X)*X'*Y;
            SSR=B'*X'*Y;
            SST=Y'*Y;
            R2=SSR/SST;
            
            % reconstruct the contour using the estimated coefficients
            Cfit=x(1).*C1+x(2).*C2;
            CfitRs=reshape(Cfit,size(real,1),[]);
            CresRs=reshape(residual,size(real,1),[]);
            
            % normalize canal coefficients
            xn=x./sum(x);
            
            
            case 2  % mixed orthogonal
            
            C1=reshape(transField,[],1);
            Cr=reshape(matrand,[],1);
            
            %%% fit four translation channels
            % non linear optimization
            guess=0.5;
            lowerBound=0;
            upperBound=1;
            regressor=C1;
            
            options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
            [x,resnorm,residual,exitflag,output]=lsqnonlin(@(guess)myfunSurfDynamic(guess,regressor,Cr),guess,lowerBound,upperBound,options)
            
            % multiple linear regression
            [b,bint,r,rint,stats]=regress(Cr,[regressor]);
            
            % multiple regrression without a constant from my dimentionality of color
            % vision work based on Maloney's paper.
            X=regressor;
            Y=Cr;
            B=X\Y;          % this the same as: B=inv(X'*X)*X'*Y;
            SSR=B'*X'*Y;
            SST=Y'*Y;
            R2=SSR/SST;
            
            % reconstruct the contour using the estimated coefficients
            Cfit=x(1).*C1;
            CfitRs=reshape(Cfit,size(real,1),[]);
            CresRs=reshape(residual,size(real,1),[]);
            
            % normalize canal coefficients
            xn=x./sum(x);
            
            
            
            
        case 1 % rotation
            
            C1=reshape(mat180ASC,[],1);
            C2=reshape(mat180PSC,[],1);
            C3=reshape(mat180LSC,[],1);
            Cr=reshape(mat180rand,[],1);
            try Ci1=reshape(mat180iASC,[],1); end
            try Ci2=reshape(mat180iPSC,[],1); end
            try Ci3=reshape(mat180iLSC,[],1); end
            
            switch polarity2
                case 1
                    % non linear optimization
                    guess=[0.5,0.5,0.5];
                    options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
                    [x180,resnorm180,residual180,exitflag180,output180]=lsqnonlin(@(guess) myfunSurf(guess,C1,C2,C3,Cr),guess,[0,0,0],[1,1,1],options)
                    
                    % multiple linear regression
                    [b180,bint180,r180,rint180,stats180]=regress(Cr,[C1,C2,C3]);
                    
                    % multiple regrression without a constant from my dimentionality of color
                    % vision work based on Maloney's paper.
                    X=[C1,C2,C3];
                    Y=Cr;
                    B=X\Y;          % this the same as: B=inv(X'*X)*X'*Y;
                    SSR=B'*X'*Y;
                    SST=Y'*Y;
                    R2_180=SSR/SST;
                    
                    % reconstruct the contour using the estimated coefficients
                    Cfit=x180(1).*C1+x180(2).*C2+x180(3).*C3;
                    %Cfit=b180(1).*C1+b180(2).*C2+b180(3).*C3;
                    CfitRs180=reshape(Cfit,size(mat180ASC,1),[]);
                    
                    % normalize canal coefficients
                    xn180=x180./sum(x180);
                    
                case 2
                    
                    % non linear optimization
                    guess=[0.5,0.5,0.5,0.5,0.5,0.5];
                    options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
                    [x180,resnorm180,residual180,exitflag180,output180]=lsqnonlin(@(guess) myfunSurf6c(guess,C1,C2,C3,Ci1,Ci2,Ci3,Cr),guess,[0,0,0,0,0,0],[1,1,1,1,1,1],options)
                    
                    % multiple linear regression
                    [b180,bint180,r180,rint180,stats180]=regress(Cr,[C1,C2,C3,Ci1,Ci2,Ci3]);
                    
                    % multiple regrression without a constant from my dimentionality of color
                    % vision work based on Maloney's paper.
                    X=[C1,C2,C3,Ci1,Ci2,Ci3];
                    Y=Cr;
                    B=X\Y;          % this the same as: B=inv(X'*X)*X'*Y;
                    SSR=B'*X'*Y;
                    SST=Y'*Y;
                    R2_180=SSR/SST;
                    
                    % reconstruct the contour using the estimated coefficients
                    Cfit=x180(1).*C1+x180(2).*C2+x180(3).*C3+x180(4).*Ci1+x180(5).*Ci2+x180(6).*Ci3;
                    %Cfit=b180(1).*C1+b180(2).*C2+b180(3).*C3;
                    CfitRs180=reshape(Cfit,size(mat180ASC,1),[]);
                    
                    % normalize canal coefficients
                    xn180=x180./sum(x180);
                    
            end
    end
    null.coeff(i,:)=xn;
    null.R2(i,1)=R2;
    null.b(i,:)=b';
    null.x(i,:)=x;
    null.FileName{i,1}=FileName;
    %fitRandSummary{i}=ws2struct();
    
    f1=figure;
    alpha=120:resolution:300;
    beta=0:resolution:90;
    contourf(alpha,beta,matrand);
    %caxis([1 17]);
    
end         % loop over all iterations

%% save all data files

filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];

currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
switch ALPHACORRstate
    case 1
        switch histType
            case 10
                save(['nullAlphaCorr_',num2str(filename),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'null');
                %save(['fitRandSummaryAlphaCorr_',num2str(filename),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'fitRandSummary');
        end
        
    case 0
        switch histType
            case 10
                save(['null_',num2str(filename),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'null');
                %save(['fitRandSummary_',num2str(filename),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'fitRandSummary');
        end
end

end

