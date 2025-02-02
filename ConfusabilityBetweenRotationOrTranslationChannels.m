
                                
% load the translational data (net spiking in each channel as a
% function of az/el of translation axis)                                    

M = zeros(181,361,8);
 
load('TransR1.mat');
M(:,:,1) = allPossibleFields.meanOutput;

load('TransR2.mat');
M(:,:,2) = allPossibleFields.meanOutput;

load('TransR3.mat');
M(:,:,3) = allPossibleFields.meanOutput;

load('TransR4.mat');
M(:,:,4) = allPossibleFields.meanOutput;

load('TransL1.mat');
M(:,:,5) = allPossibleFields.meanOutput;

load('TransL2.mat');
M(:,:,6) = allPossibleFields.meanOutput;

load('TransL3.mat');
M(:,:,7) = allPossibleFields.meanOutput;

load('TransL4.mat');
M(:,:,8) = allPossibleFields.meanOutput;


M(isnan(M))=0;  %convert all NaN to zero


% load the rotational data (net spiking in each channel as a
% function of az/el of rotational axis)        

N = zeros(181,361,8);
 
load('RotR1.mat');
N(:,:,1) = allPossibleFields.meanOutput;

load('RotR2.mat');
N(:,:,2) = allPossibleFields.meanOutput;

load('RotR3.mat');
N(:,:,3) = allPossibleFields.meanOutput;

load('RotR4.mat');
N(:,:,4) = allPossibleFields.meanOutput;

load('RotL1.mat');
N(:,:,5) = allPossibleFields.meanOutput;

load('RotL2.mat');
N(:,:,6) = allPossibleFields.meanOutput;

load('RotL3.mat');
N(:,:,7) = allPossibleFields.meanOutput;

load('RotL4.mat');
N(:,:,8) = allPossibleFields.meanOutput;

%N(isnan(N))=0;  %convert all NaN to zero


%% calculate confusability 

% CHANGE HERE
%%%%%%%%%%%%%%%%%%%%%%%
refmat=N;
testmat=M;
%%%%%%%%%%%%%%%%%%%%%%%

EuclidDist=zeros(181,361);
Confusion=zeros(181,361);
Confusion_lowres=zeros(18,36);

for C_elev=1:10:181       %restrict range of elevation to avoid glitches at 
    for C_az=1:10:361     %start at 2 instead of 1 to avoid glitch of zeros
        Ref_heat_values=refmat(C_elev,C_az,:); % get all 8 channel strengths for the reference azim-elev in the translation map 
                                                
        for W = 1:181
            for V = 1:361

                Test_heat_values=testmat(W,V,:);            % get the 8 channel strengths for all test azim-elevs in the rotation map 
                EuclidDist(W,V)=(sum((Ref_heat_values - Test_heat_values).^2))^0.5;
                
            end
        end
        Confusion(C_elev,C_az)=min(min(EuclidDist));
        Confusion_lowres(round(C_elev/10+1), round(C_az/10+1))=min(min(EuclidDist));  
    end
    C_elev    %monitor progress during execution
end

save('tmp.mat');


%% determine the theoretical max Euclidian distance for 8 dimensions
% (assumes range of possible outputs for any channel is 0 to 1)
%p = [0 0 0 0 0 0 0 0];
%q = [1 1 1 1 1 1 1 1];
%max_heat_theoretical = (sum((p - q).^2))^0.5; % Euclid dist for most disparate
                                        % values across all 6 channels
%max_heat = plot_heat_range*max_heat_theoretical;    % to set max for plots using user entered scale factor  
                                        %nothing ever approaches theoretical max 
                                        %so tone down but keep fixed across both plots

%plot the heat map for Euclidean distance for translational axes
%figure;
%imagesc(EuclidDist_Transl, [0 max_heat])
%title(strcat('TRANSLATION: Euclid distance from translation toward Az= ', num2str(Azimuth) ,'   El= ',num2str(Elevation)  ));

%plot the heat map for Euclidean distance for rotational axes
%figure;
%imagesc(EuclidDist_Rot, [0 max_heat])
%title(strcat('ROTATION: Euclid distance from translation toward Az= ', num2str(Azimuth),'   El= ',num2str(Elevation)  ));

% Confusion_lowres=circshift(Confusion_lowres,-1,1);
% Confusion_lowres=circshift(Confusion_lowres,-1,2);

figure;   %vertical flip of matrix needed to correct for inverted convention when using imagesc
x=[-179 179];
y=[-89 89];
imagesc(x,y,flipud(Confusion_lowres))
xlim([-180 180]); ylim([-90 90]);


Confusion( ~any(Confusion,2), : ) = [];  %remove rows where all elements are zeros
Confusion( :, ~any(Confusion,1) ) = [];  %remove columns where all elements are zeros

figure;   %vertical flip of matrix needed to correct for inverted convention when using imagesc
x=[-179 179];
y=[-89 89];
imagesc(x,y,flipud(Confusion));
xlim([-180 180]); ylim([-90 90]);
azimuthElevationPlotSettings

c=colorbar('northoutside');
c.Label.String='Discriminability';
c.FontSize=14;
caxis([0 0.8]);



title('discriminatability of 8 ch outputs to translation with the most similar rotational output')
disp ('Minimum Euclidian distance to any rotationally generated 8-channel output:  ')
min(min(Confusion_lowres))









