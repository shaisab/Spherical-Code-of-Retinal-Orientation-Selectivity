function []= fitSurf2mixed(channelID,resolution,mapType)
%% import data
% channelID is equivalent to polarity2
switch mapType
    case 1          % rotation canal-based
        switch channelID
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
                
            case 3
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
                
            case 4      % rotation - singularities based on real data
                
                [FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of single-channel rotation maps','MultiSelect','on');
                
                % sort file order by ascending alpha
                for k=1:size(FileName,2)
                    indUnder=find(FileName{k}=='_');
                    order(k)=str2num(FileName{k}(indUnder(6)+1:indUnder(7)-1));     % for individual retinas
                    %order(k)=str2num(FileName{k}(indUnder(4)+1:indUnder(5)-1));     % for all retinas pooled
                    %order(k)=str2num(FileName{k}(indUnder(3)+1:indUnder(4)-1));     % use when the text '_Calcium imaging data' was omitted from file name 
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
                
                
        end
        
    case 0      % translation - singularities based on real data
        
        [FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of single-channel translation or rotation maps','MultiSelect','on');
        
        % sort file order by ascending alpha
        for k=1:size(FileName,2)
            indUnder=find(FileName{k}=='_');
            order(k)=str2num(FileName{k}(indUnder(6)+1:indUnder(7)-1));     % for individual retinas
            %order(k)=str2num(FileName{k}(indUnder(4)+1:indUnder(5)-1));     % for all retinas pooled
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
        
        
        case 2      % mixed orthogonal - singularities based on real data
        
        [FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of single-channel translation or mixed orthogonal maps','MultiSelect','on');
        load(FileName);
        transField=generateMatchindMatrixMix(allPossibleFields);       
       
        
        [FileName,PathName] = uigetfile('*.mat','Select allPossibleFields file of translation or mixed orthogonal of the real data','MultiSelect','on');
        if strfind(FileName,'_real')
            load(FileName);
            real=generateMatchindMatrixMix(allPossibleFields);
            alphasing=allPossibleFields.alphasing;
            betasing=allPossibleFields.betasing;    
        end
        
        
    case 3      % mixed vertical translation + horizontal rotation (singularities basde on real data)
        
        [names,PathName] = uigetfile('*.mat','Select allPossibleFields vertical translation horizontal rotation maps','MultiSelect','on');
        
        % sort file order by ascending alpha
        for k=1:size(names,2)
            indUnder=find(names{k}=='_');
            if strfind(names{k},'single channel probed')
                order(k)=str2num(names{k}(indUnder(13)+1:indUnder(14)-1));
            else
                order(k)=str2num(names{k}(indUnder(4)+1:indUnder(5)-1));
            end
        end
        [~,I]=sort(order);
        nameSorted=names(I);
        
        transn=1;
        rotn=1;
        for i=1:size(nameSorted,2)
            if strfind(nameSorted{i},'translation_alpha')                               % two vertical translation channels
                load(num2str(nameSorted{i}));
                trans{1,transn}=generateMatchindMatrix(allPossibleFields);
                transn=transn+1;
            elseif strfind(nameSorted{i},'single channel probed with_translation')      % two horizontal rotation channels
                load(num2str(nameSorted{i}));
                rot{1,rotn}=generateMatchindMatrix(allPossibleFields);
                rotn=rotn+1;
            end
        end
        
        [FileName,PathName] = uigetfile('*.mat','Select allPossibleFields translation file of the real data','MultiSelect','on');
        if strfind(FileName,'_real')
            load(FileName);
            real=generateMatchindMatrix(allPossibleFields);
            alphasing=allPossibleFields.alphasing;
            betasing=allPossibleFields.betasing;    
        end
        
end

%% perform fit
switch mapType
    case 1          % rotation
                
        switch channelID
            case 1      % sense polarity for each of three canalcentric channels (3 channels) -  canal-based singularities
                C1=reshape(matASC,[],1);
                C2=reshape(matPSC,[],1); 
                C3=reshape(matLSC,[],1);
                Cr=reshape(matreal,[],1);

                % non linear optimization
                guess=[0.5,0.5,0.5];
                options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
                [x,resnorm,residual,exitflag,output]=lsqnonlin(@(guess) myfunSurf(guess,C1,C2,C3,Cr),guess,[0,0,0],[1,1,1],options)
                
                % multiple linear regression
                [b,bint,r,rint,stats]=regress(Cr,[C1,C2,C3]);
                
                % multiple regrression without a constant from my dimentionality of color
                % vision work based on Maloney's paper.
                X=[C1,C2,C3];
                Y=Cr;
                B=X\Y;          % this the same as: B=inv(X'*X)*X'*Y;
                SSR=B'*X'*Y;
                SST=Y'*Y;
                R2=SSR/SST;
                
                % reconstruct the contour using the estimated coefficients
                Cfit=x(1).*C1+x(2).*C2+x(3).*C3;
                %Cfit=b180(1).*C1+b180(2).*C2+b180(3).*C3;
                CfitRs=reshape(Cfit,size(matASC,1),[]);
                
                % normalize canal coefficients
                xn=x./sum(x);
                
            case 2          % sense and antisense polarity for each of three canalcentric channels (6 channels) -  canal-based singularities
                
            C1=reshape(matASC,[],1);
            C2=reshape(matPSC,[],1);
            C3=reshape(matLSC,[],1);
            Cr=reshape(matreal,[],1);
            Ci1=reshape(matiASC,[],1); 
            Ci2=reshape(matiPSC,[],1);
            Ci3=reshape(matiLSC,[],1);
        
                % non linear optimization
                guess=[0.5,0.5,0.5,0.5,0.5,0.5];
                options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
                [x,resnorm,residual,exitflag,output]=lsqnonlin(@(guess) myfunSurf6c(guess,C1,C2,C3,Ci1,Ci2,Ci3,Cr),guess,[0,0,0,0,0,0],[1,1,1,1,1,1],options)
                
                % multiple linear regression
                [b,bint,r,rint,stats]=regress(Cr,[C1,C2,C3,Ci1,Ci2,Ci3]);
                
                % multiple regrression without a constant from my dimentionality of color
                % vision work based on Maloney's paper.
                X=[C1,C2,C3,Ci1,Ci2,Ci3];
                Y=Cr;
                B=X\Y;          % this the same as: B=inv(X'*X)*X'*Y;
                SSR=B'*X'*Y;
                SST=Y'*Y;
                R2=SSR/SST;
                
                % reconstruct the contour using the estimated coefficients
                Cfit=x(1).*C1+x(2).*C2+x(3).*C3+x(4).*Ci1+x(5).*Ci2+x(6).*Ci3;
                %Cfit=b180(1).*C1+b180(2).*C2+b180(3).*C3;
                CfitRs=reshape(Cfit,size(matASC,1),[]);
                
                % normalize canal coefficients
                xn=x./sum(x);
                
                
            case 3          % sense and antisense polarity for ASC and LSC channels (4 channels) -  canal-based singularities
                
            C1=reshape(matASC,[],1);
            C2=reshape(matPSC,[],1);
            C3=reshape(matLSC,[],1);
            Cr=reshape(matreal,[],1);
            Ci1=reshape(matiASC,[],1); 
            Ci2=reshape(matiPSC,[],1);
            Ci3=reshape(matiLSC,[],1);
   
                % non linear optimization
                guess=[0.5,0.5,0.5,0.5];
                lowerBound=[0,0,0,0];
                upperBound=[1,1,1,1];
                regressor=[C1,C3,Ci1,Ci3];
                
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
                Cfit=x(1).*C1+x(2).*C3+x(3).*Ci1+x(4).*Ci3;
                CfitRs=reshape(Cfit,size(matASC,1),[]);
                CresRs=reshape(residual,size(matASC,1),[]);
                
                % normalize canal coefficients
                xn=x./sum(x);
                
                % non linear optimization on residuals using the two PSC channels
                guess=[0.5,0.5];
                lowerBound=[0,0];
                upperBound=[1,1];
                regressor=[C2,Ci2];
                
                options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
                [PSCx,PSCresnorm,PSCresidual,PSCexitflag,PSCoutput]=lsqnonlin(@(guess)myfunSurfDynamic(guess,regressor,Cr),guess,lowerBound,upperBound,options)
                
                % multiple regrression without a constant from my dimentionality of color
                % vision work based on Maloney's paper.
                X=regressor;
                Y=Cr;
                B=X\Y;          % this the same as: B=inv(X'*X)*X'*Y;
                SSR=B'*X'*Y;
                SST=Y'*Y;
                PSC_R2=SSR/SST;
                
                % reconstruct the contour using the estimated coefficients
                PSC_Cfit=x(1).*C2+x(2).*Ci2;
                PSC_CfitRs=reshape(PSC_Cfit,size(matASC,1),[]);
                PSC_CresRs=reshape(PSCresidual,size(matASC,1),[]);
                
                % normalize canal coefficients
                PSCxn=PSCx./sum(PSCx);
                
            case 4          % four rotation channels based on real singularities
                
                C1=reshape(transField{1},[],1);
                C2=reshape(transField{2},[],1);
                Cr=reshape(real,[],1);
                
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
                
                
        end
        
        
        
        
        
        
        
    case 0          % translation
        
        C1=reshape(transField{1},[],1);
        C2=reshape(transField{2},[],1);
        C3=reshape(transField{3},[],1);
        C4=reshape(transField{4},[],1);
        Cr=reshape(real,[],1);
        
        %%% fit four translation channels
        % non linear optimization
        guess=[0.5,0.5,0.5,0.5];
        lowerBound=[0,0,0,0];
        upperBound=[1,1,1,1];
        regressor=[C1,C2,C3,C4];
        
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
        
        
        
        case 2          % mixed orthogonal
        
        C1=reshape(transField,[],1);
        Cr=reshape(real,[],1);
        
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
        R2_adjusted=1-((1-R2).*((length(Cr)-1)/(length(Cr)-1)));
        
        % reconstruct the contour using the estimated coefficients
        Cfit=x(1).*C1;
        CfitRs=reshape(Cfit,size(real,1),[]);
        CresRs=reshape(residual,size(real,1),[]);
        
        % normalize canal coefficients
        xn=x./sum(x);
        
        
    case 3          % mixed vertical translation + horizontal rotation
        
        C1=reshape(trans{1,1},[],1);
        C2=reshape(trans{1,2},[],1);
        C3=reshape(rot{1,1},[],1);
        C4=reshape(rot{1,2},[],1);
        Cr=reshape(real,[],1);
        
        %%% fit four translation channels
        % non linear optimization
        guess=[0.5,0.5,0.5,0.5];
        lowerBound=[0,0,0,0];
        upperBound=[1,1,1,1];
        regressor=[C1,C2,C3,C4];
        
        options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
        %options=optimset('MaxFunEvals',1e+20,'MaxIter',5000,'TolFun',1e-10,'TolX',1e-10,'FinDiffType','central','Display','iter');
        [x,resnorm,residual,exitflag,output]=lsqnonlin(@(guess)myfunSurfDynamic(guess,regressor,Cr),guess,lowerBound,upperBound,options)
        
        % multiple linear regression
        %[b,bint,r,rint,stats]=regress(Cr,[regressor]);
        
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
        
        
end


%% calculate alpha and beta of canals (rotation) OR of observed hot points (translation)

if mapType==1           % rotation
    
    normalASC=[0.592,0.764,0.247];          %ASC
    [alphaASC,betaASC]=AB(normalASC);
    alphaASC=rad2deg(alphaASC);
    betaASC=rad2deg(betaASC);
    betaiASC=90-abs(90-betaASC);
    if alphaASC<=180
        alphaiASC=180+alphaASC;
    else alphaASC>180
        alphaiASC=alphaASC-180;
    end
    
    normalPSC=[0.470,-0.764,0.438];         %PSC
    [alphaPSC,betaPSC]=AB(normalPSC);
    alphaPSC=rad2deg(alphaPSC);
    betaPSC=rad2deg(betaPSC);
    betaiPSC=90-abs(90-betaPSC);
    if alphaPSC<=180
        alphaiPSC=180+alphaPSC;
    else alphaPSC>180
        alphaiPSC=alphaPSC-180;
    end
    
    normalLSC=[0.421,-0.056,-0.901];        %LSC
    [alphaLSC,betaLSC]=AB(normalLSC);
    alphaLSC=rad2deg(alphaLSC);
    betaLSC=rad2deg(betaLSC);
    betaiLSC=90-abs(90-betaLSC);
    if alphaLSC<=180
        alphaiLSC=180+alphaLSC;
    else alphaLSC>180
        alphaiLSC=alphaLSC-180;
    end
    
end

%% plot fit, real, and residual contours

% plot model contour
f1=figure;
alpha=120:resolution:300;
beta=0:resolution:90;
subplot(2,2,1);
contourf(alpha,beta,CfitRs,10);
%contourf(alpha,beta,CfitRs180,10);

switch mapType
    case 1
        switch channelID
            case 1              % rotation
                title(['model  ','ASC:',num2str(xn(1),2),...
                    ',  PSC:',num2str(xn(2),2),',  LSC:',num2str(xn(3),2),',  R^2: ',num2str(R2,2)]);
            case 2
                title(['model  ','ASC:',num2str(xn(1),2),...
                    ',  PSC:',num2str(xn(2),2),',  LSC:',num2str(xn(3),2),'  iASC:',num2str(xn(4),2),...
                    ',  iPSC:',num2str(xn(5),2),'  iLSC:',num2str(xn(6),2),',  R^2: ',num2str(R2,2)]);
            case 3
                title(['model  ','ASC:',num2str(xn(1),2),...
                    ',  LSC:',num2str(xn(2),2),'  iASC:',num2str(xn(3),2),...
                    '  iLSC:',num2str(xn(4),2),',  R^2: ',num2str(R2,2)]);
            case 4
                title(['model  ','dorsal: ',num2str(xn(1),2),...
            ',  temporal: ',num2str(xn(2),2),'  ventral: ',num2str(xn(3),2),...
            '  nasal: ',num2str(xn(4),2),',  R^2: ',num2str(R2,2)]);
        hold on;
        scatter(alphasing,betasing,'bo','LineWidth',3);
        hold off;
        
        end
        hold on;
        if channelID~=4
        scatter(alphaASC,betaASC,'bx','LineWidth',3);
        scatter(alphaiASC,betaiASC,'bo','LineWidth',3);
        scatter(alphaPSC,betaPSC,'gx','LineWidth',3);
        scatter(alphaiPSC,betaiPSC,'go','LineWidth',3);
        scatter(alphaLSC,betaLSC,'rx','LineWidth',3);
        scatter(alphaiLSC,betaiLSC,'ro','LineWidth',3);
        end
        hold off;
        
    case 0              % translation
        
        title(['model  ','nasal: ',num2str(xn(1),2),...
            ',  dorsal: ',num2str(xn(2),2),'  temporal: ',num2str(xn(3),2),...
            '  ventral: ',num2str(xn(4),2),',  R^2: ',num2str(R2,2)]);
        hold on;
        scatter(alphasing,betasing,'bo','LineWidth',3);
        hold off;
        
        
        case 2              % mixed orthogonal
        
        title(['model  ','nasal: ',num2str(xn(1),2),',  R^2: ',num2str(R2_adjusted,2)]);
        hold on;
        scatter(alphasing,betasing,'bo','LineWidth',3);
        hold off;
        
    case 3              % mixed translation-rotation
        
        title(['model  ','trans field 1: ',num2str(xn(1),2),...
            ',  trans field 2: ',num2str(xn(2),2),'  rot field 3: ',num2str(xn(3),2),...
            '  rot field 4: ',num2str(xn(4),2),',  R^2: ',num2str(R2,2)]);
        hold on;
        scatter(alphasing,betasing,'bo','LineWidth',3);
        hold off;
        
end

% plot real contour
subplot(2,2,3);
contourf(alpha,beta,real);
title('real');
hold on;
switch mapType
    case 1          % rotation
        scatter(alphaASC,betaASC,'bx','LineWidth',3);
        scatter(alphaiASC,betaiASC,'bo','LineWidth',3);
        scatter(alphaPSC,betaPSC,'gx','LineWidth',3);
        scatter(alphaiPSC,betaiPSC,'go','LineWidth',3);
        scatter(alphaLSC,betaLSC,'rx','LineWidth',3);
        scatter(alphaiLSC,betaiLSC,'ro','LineWidth',3);
    case 0          % translation
        scatter(alphasing,betasing,'bo','LineWidth',3);
    case 2          % mixed orthogonal
        scatter(alphasing,betasing,'bo','LineWidth',3);
    case 3          % mixed old
        scatter(alphasing,betasing,'bo','LineWidth',3);
        hold off;
end

% plot bar plot of fit coefficients
subplot(2,2,2);
ax = gca;
switch mapType
    case 1          % rotation
        switch channelID
            case 1
                bar([xn(1),xn(2),xn(3)],'k');
                ax.XTickLabel={'ASC','PSC','LSC'};
                xlim([0.5 3.5]);
            case 2
                bar([xn(1),xn(4),xn(2),xn(5),xn(3),xn(6)],'k');
                ax.XTickLabel={'+ASC','-ASC','+PSC','-PSC','+LSC','-LSC'};
                xlim([0.5 6.5]);
            case 3
                bar([xn(1),xn(2),xn(3),xn(4)],'k');
                ax.XTickLabel={'+ASC','-ASC','+LSC','-LSC'};
                xlim([0.5 4.5]);
            case 4
            bar([xn(1),xn(2),xn(3),xn(4)],'k');     %for ON-OFF rotation fit where the alpha of the nasal axis is 350 deg
            bar([xn(2),xn(3),xn(4),xn(1)],'k');     %for ON rotation fit where the alpha of the nasal axis is 5 deg
            %ax.XTickLabel={'dorsal','temporal','ventral','nasal'};
            ax.XTickLabel={'nasal','dorsal','temporal','ventral'};              
            
            xlim([0.5 4.5]);
        end
    case 0          % translation
        bar([xn(1),xn(2)],0.8,'grouped');
        xlim([0.65 2.4]);
        
    case 2          % mixed orthogonal
        bar(xn(1),0.8,'grouped');
        xlim([0.65 2.4]);
        
    case 3          % mixed translation-rotation
        bar([xn(1),xn(2),xn(3),xn(4)],'k');
        ax.XTickLabel={'trans field 1','trans field 2','rot field 3','rot field 4'};
        %ax.XTickLabel={'nasal','dorsal','temporal','ventral'};
        xlim([0.5 4.5]);
end
set(gca,'FontSize',16)
xlabel('Channel','FontSize',20);ylabel('Relative contribution','FontSize',20);
set(ax,'box','off');
box off
title(['\itR\rm^2 = ',num2str(R2,2)]);


if channelID==3 || mapType==0 || mapType==2
    % plot residuals
    subplot(2,2,4);
    contourf(alpha,beta,CresRs);
    title('residuals');
    hold on;
    switch mapType
        case 1          % rotation
            scatter(alphaASC,betaASC,'bx','LineWidth',3);
            scatter(alphaiASC,betaiASC,'bo','LineWidth',3);
            scatter(alphaPSC,betaPSC,'gx','LineWidth',3);
            scatter(alphaiPSC,betaiPSC,'go','LineWidth',3);
            scatter(alphaLSC,betaLSC,'rx','LineWidth',3);
            scatter(alphaiLSC,betaiLSC,'ro','LineWidth',3);
        case 0          % translation
            scatter(alphasing,betasing,'bo','LineWidth',3);
        case 2          % mixed orthogonal
            scatter(alphasing,betasing,'bo','LineWidth',3);
    end
    hold off;
end

if channelID==3
    % plot PSC model contour
    f2=figure;
    alpha=0:resolution:360;
    beta=0:resolution:180;
    subplot(2,2,1);
    contourf(alpha,beta,PSC_CfitRs);
    title(['model  ','PSC:',num2str(PSCxn(1),2),...
        ',  iPSC:',num2str(xn(2),2),',  R^2: ',num2str(PSC_R2,2)]);
    hold on;
    scatter(alphaASC,betaASC,'bx','LineWidth',3);
    scatter(alphaiASC,betaiASC,'bo','LineWidth',3);
    scatter(alphaPSC,betaPSC,'gx','LineWidth',3);
    scatter(alphaiPSC,betaiPSC,'go','LineWidth',3);
    scatter(alphaLSC,betaLSC,'rx','LineWidth',3);
    scatter(alphaiLSC,betaiLSC,'ro','LineWidth',3);
    hold off;
    
    
    % plot residuals after fitting the PSC
    subplot(2,2,2);
    contourf(alpha,beta,PSC_CresRs);
    title('residuals after fitting PSC');
    hold on;
    scatter(alphaASC,betaASC,'bx','LineWidth',3);
    scatter(alphaiASC,betaiASC,'bo','LineWidth',3);
    scatter(alphaPSC,betaPSC,'gx','LineWidth',3);
    scatter(alphaiPSC,betaiPSC,'go','LineWidth',3);
    scatter(alphaLSC,betaLSC,'rx','LineWidth',3);
    scatter(alphaiLSC,betaiLSC,'ro','LineWidth',3);
    hold off;
    
end
%% save file and figure
filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];

currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
save(['fitAlphaCorr_',num2str(filename),'_DS_created_',currenttime]);

%saveas(f1,['fitAlphaCorr_',num2str(filename),'_cutoff_',num2str(10),'_DS_created_',currenttime,'.fig']);



end

