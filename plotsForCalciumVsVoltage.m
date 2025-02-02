
% pllot for ONDS
V=physdata.rep{1, 1}.ordResponse{1, 2}(:,7);            % from June 22, 2015, fov#2, run 2, direction 7
C=pResponse;

figure;
subplot(2,1,1)
plot(Time+2,C);
xlim([0 20]);
subplot(2,1,2)
plot(physdata.rep{1, 1}.phystime(:,1),V);
xlim([0 20]);

% pllot for ONOFFDS
V=physdata.rep{1, 1}.ordResponse{1, 1}(:,3);            % from June 22, 2015, fov#2, run 2, direction 7
C=pResponse;

figure;
subplot(2,1,1)
plot(Time+2,C);
xlim([0 20]);
subplot(2,1,2)
plot(physdata.rep{1, 1}.phystime(:,1),V);
xlim([0 20]);


% figure;
% [ax,p1,p2] = plotyy(Time+2,C,physdata.rep{1, 1}.phystime(:,1),V,'plot','plot');