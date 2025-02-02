% extract data from figure

 h = gco; %select line in the plot
 xdata=get(h,'XData');
 ydata=get(h,'YData');
 data=[xdata',ydata'];
 cdata=get(h,'CData');
 