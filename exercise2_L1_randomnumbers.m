% excercise #2 Lecture #1 CH480/680

figure(1); clf; figure(2); clf; figure(3); clf; figure(4); clf;figure(5); clf

%octave/matlab can generate random numbers with a normal distribution
%let's write a script to generate 2, 3, 10, 100, 1000 random numbers with
% a mean of 1.5 and a standard deviation of 0.25
% calculate the sample mean and standard deviation for these sets
% plot histograms for each data set and the theoretical distributions on the 
% same plot

% define the normal function

f = @(x,mu,sigma) (1/(sigma*sqrt(2*pi)))*exp((-1/2)*((x-mu)./sigma).^2);

% simulate for 2, 3, 10, 100 and 1000

simdata2=0.25*(randn(1,2))+1.5;
simdata3=0.25*(randn(1,3))+1.5;
simdata10=0.25*(randn(1,10))+1.5;
simdata100=0.25*(randn(1,100))+1.5;
simdata1000=0.25*(randn(1,1000))+1.5;

average2=mean(simdata2); average3=mean(simdata3); average10=mean(simdata10);
average100=mean(simdata100); average1000=mean(simdata1000);

stdev2=std(simdata2); stdev3=std(simdata3); stdev10=std(simdata10); 
stdev100=std(simdata100); stdev1000=std(simdata1000); 

display=[average2 average3 average10 average100 average1000
    stdev2 stdev3 stdev10 stdev100 stdev1000]

% now make the plots
% generate the theoretical curve  (area of 1)
x=0:0.01:3; mu=1.5; sigma=(0.25);

% call the function, make sure variables are in the right order
y=f(x,mu,sigma);

figure(1);
xlabel('x'); ylabel('frequency'); title('2 random numbers')
hold on; hist(simdata2,2); [counts,centers]=hist(simdata2,2);
binwidth=centers(2)-centers(1);
%[Mdata,I]=max(counts); Msim=max(y); factor=Mdata/Msim; y2=factor*y;
%plot(x,y2,'linewidth',2); set(gca,'linewidth',2,'fontsize',11)
plot(x,y*2*binwidth,'linewidth',2); set(gca,'linewidth',2,'fontsize',11)

figure(2);
xlabel('x'); ylabel('frequency'); title('3 random numbers')
hold on; hist(simdata3,3); [counts,centers]=hist(simdata3,3);
binwidth=centers(2)-centers(1);
%[Mdata,I]=max(counts); Msim=max(y); factor=Mdata/Msim; y3=factor*y;
%plot(x,y3,'linewidth',2); set(gca,'linewidth',2,'fontsize',11)
plot(x,y*3*binwidth,'linewidth',2); set(gca,'linewidth',2,'fontsize',11)

figure(3);
xlabel('x'); ylabel('frequency'); title('10 random numbers')
hold on; hist(simdata10,5); [counts,centers]=hist(simdata10,5);
binwidth=centers(2)-centers(1);
%[Mdata,I]=max(counts); Msim=max(y); factor=Mdata/Msim; y10=factor*y;
%plot(x,y10,'linewidth',2); set(gca,'linewidth',2,'fontsize',11)
plot(x,y*10*binwidth,'linewidth',2); set(gca,'linewidth',2,'fontsize',11)

figure(4);
xlabel('x'); ylabel('frequency'); title('100 random numbers')
hold on; hist(simdata100,10); [counts,centers]=hist(simdata100,10);
binwidth=centers(2)-centers(1);
%[Mdata,I]=max(counts); Msim=max(y); factor=Mdata/Msim; y100=factor*y;
%plot(x,y100,'linewidth',2); set(gca,'linewidth',2,'fontsize',11)
plot(x,y*100*binwidth,'linewidth',2); set(gca,'linewidth',2,'fontsize',11)

figure(5); % 1000 sims
xlabel('x'); ylabel('frequency'); title('1000 random numbers')
hold on; hist(simdata1000,10); [counts,centers]=hist(simdata1000,10); 
binwidth=centers(2)-centers(1);
%[Mdata,I]=max(counts); Msim=max(y); factor=Mdata/Msim; y1000=factor*y;
%plot(x,y1000,'linewidth',2); set(gca,'linewidth',2,'fontsize',11)
plot(x,y*1000*binwidth,'linewidth',2); set(gca,'linewidth',2,'fontsize',11)

% or using the built in function!

%for Octave use this line. for matlab comment it out.
%pkg load statistics

figure(6); histfit(simdata1000)
xlabel('x'); ylabel('frequency'); title('1000 random numbers')

