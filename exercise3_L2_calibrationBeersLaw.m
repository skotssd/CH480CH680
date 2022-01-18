% for the Beer's law calibration determine the slope and uncertainty 
% in the slope estimation. Include a plot of the calibration curve with model line.

% for octave

% pkg load statistics

% input the data x=conc, y=Abs;

x=[0 1 2 4 5];
y=[0 0.41 0.78 1.64 1.98];

% determine the best fit slope

slope=(sum(y.*x))./(sum(x.^2))

% determine standard error

SSE_Best=sum((y-slope*x).^2);
P=1; N=length(y); nu=N-P; ssquared=SSE_Best/nu;
SE=sqrt(ssquared/sum(x.^2));

alpha=0.95; alphatable=alpha+(1-alpha)/2;
ttable=tinv(alphatable,nu);

CI95=SE*ttable
slopewitherror=[slope-CI95 slope slope+CI95]

% plot with the best answer

model=slope*x; modellow=(slope-CI95)*x; modelhigh=(slope+CI95)*x; 

figure(1)
h=plot(x,y,'ko',x,model,'k-',x,modelhigh,'k--',x,modellow,'k--');
set(h(1),'markersize',8,'markerfacecolor','b');
set(h(2),'linewidth',2)
legend('data','model','high model','low model','location','northwest')
set(gca,'linewidth',2,'fontsize',12)
xlabel('concentration (ppm)'); ylabel('absorbance')



