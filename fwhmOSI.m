function oriFWHM = fwhmOSI(xqTmp,RqTmp)

% Find the half max value.
halfMax = (min(RqTmp) + max(RqTmp)) / 2;
% Find where the data first drops below half the max.
index1 = find(RqTmp >= halfMax, 1, 'first');
% Find where the data last rises above half the max.
index2 = find(RqTmp >= halfMax, 1, 'last');
fwhm = index2-index1 + 1; % FWHM in indexes.
% OR, if you have an x vector
oriFWHM = xqTmp(index2) - xqTmp(index1);

end