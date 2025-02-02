a = (allCells.isRBPMS);
[C,ia,ic] = unique(a);
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];