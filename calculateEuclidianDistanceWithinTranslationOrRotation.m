%      ####### ENTER REFERENCE AXIS HERE  ########
Azimuth = 0;
Elevation = 45;
plot_heat_range = 0.5; % fraction of theoretical max Euclid distance to 
                        %scale all plots to

RefAxis = [Elevation+90 Azimuth+180];  % note that order here is Elevation, Azimuth 
                                    %converting to form appropriate to
                                    %input matrices


                                    
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



N(isnan(N))=0;  %convert all NaN to zero

%create matrices to hold Euclidian distance results
EuclidDist_Transl = zeros(181,361);
EuclidDist_Rot = zeros(181,361);

Ref_heat_values = M(RefAxis(1),RefAxis(2),:); % get all 8 channel strengths 
                                                % for the reference azim-elev in the translation map

for W = 1:181
    for V = 1:361
        Test_heat_values_transl = M(W,V,:);   % get the 8 channel strengths for all test azim-elevs in the translation map                 
        EuclidDist_Transl(W,V) = (sum((Ref_heat_values - Test_heat_values_transl).^2))^0.5;  % put the Euclid distance value 
                                                                                             % for each other translational
                                                                                             % flow axis into the output array
                                                                                             % for plotting
    end
end

for W = 1:181
    for V = 1:361
       
        Test_heat_values_rot = N(W,V,:);            % get the 8 channel strengths for all test azim-elevs in the rotation map 
        EuclidDist_Rot(W,V) = (sum((Ref_heat_values - Test_heat_values_rot).^2))^0.5;  
        
            
    end
end

%determine the theoretical max Euclidian distance for 8 dimensions
% (assumes range of possible outputs for any channel is 0 to 1)
p = [0 0 0 0 0 0 0 0];
q = [1 1 1 1 1 1 1 1];
max_heat_theoretical = (sum((p - q).^2))^0.5; % Euclid dist for most disparate
                                        % values across all 6 channels
                                        % Actually this is just the square root of 8 
                                        
max_heat = plot_heat_range*max_heat_theoretical;    % to set max for plots using user entered scale factor  
                                        %nothing ever approaches theoretical max 
                                        %so tone down but keep fixed across both plots

%plot the heat map for Euclidean distance for translational axes
x=[-180 180];
y=[-90 90];

figure;
imagesc(x,y,flipud(EuclidDist_Transl), [0 max_heat])
azimuthElevationPlotSettings
title(strcat('TRANSLATION: Euclid distance from translation toward Az= ', num2str(Azimuth) ,'   El= ',num2str(Elevation)  ));

%plot the heat map for Euclidean distance for rotational axes; flip needed
%for using imageSC to square with Shai data
figure;
imagesc(x,y,flipud(EuclidDist_Rot), [0 max_heat])
azimuthElevationPlotSettings
title(strcat('ROTATION: Euclid distance from translation toward Az= ', num2str(Azimuth),'   El= ',num2str(Elevation)  ));










