

% plot histograma nd concordance indeces for several translation axes

% input alpha and beta coordinates of translation axis
b=30;
a=70;

figure;
histogram(allPossibleFields.anglediff180{b,a},50)
title(['Eccentricity =',num2str(b),' ,  Direction =',num2str(a), ' ,  Concordance index =',num2str(allPossibleFields.matchind180(b,a).c10)])

