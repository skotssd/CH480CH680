% excercise #1 Lecture #1 CH480/680

figure(1); clf

%let's write a function to return  𝑓(𝑥)  for a given  𝑥 ,  𝜇  and  𝜎 
%then let's write a script to determine and plot distributions with mean 1.5 
% and  𝜎  0.1, 0.25 and 0.5 on the same graph
%finally, let's add an integration step to the script to show the area is always one

% define the normal function

f = @(x,mu,sigma) (1/(sigma*sqrt(2*pi)))*exp((-1/2)*((x-mu)./sigma).^2);

%define row vector of x values and mu and sigma

x=0:0.01:3; mu=1.5; sigma1=0.1; sigma2=0.25; sigma3=0.5;

% call the function, make sure variables are in the right order

y1=f(x,mu,sigma1); y2=f(x,mu,sigma2); y3=f(x,mu,sigma3);

% plot the 3 lines

plot(x,y1,x,y2,x,y3,'linewidth',2)
set(gca,'linewidth',2,'fontsize',11)
xlabel('x'); ylabel('frequency')
legend('\sigma=0.1', '\sigma=0.25','\sigma=0.5','location','northeast')

area1=trapz(x,y1)
area2=trapz(x,y2)
area3=trapz(x,y3) % didn't go far enough on the x-axis range.  but ok.
