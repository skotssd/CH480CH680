figure(1); clf; figure(2); clf

x=[2 2.5 2.5 2.75 3 3 3];
y=[89 97 91 98 100 104 97];
X=[ones(size(x')) x']; 
y=y';
b=inv(X'*X)*X'*y
P=2; % two parameters
N=length(x); %number of observatrions
nu=N-P; % degrees of freedom.  no. obs-no.parameters
bestmodel=X*b;
s2=sum((bestmodel-y).^2)./nu; s=sqrt(s2)
tvalue=tinv(0.975,nu);
seB=s*sqrt(diag(inv(X'*X)))*tvalue
w=0:0.1:2.5*pi; % for plotting the elippse

[Q,R]=qr(X); R1 = R(1:2,1:2); invR1=inv(R1); Fvalue=finv(0.95,P,N-P);

for i=1:length(w)
	    scalar=sqrt(P*s2*Fvalue);
        BETA(:,i)=b+scalar*invR1*([cos(w(i)); sin(w(i))]);
    end

%plot(BETA(1,:),BETA(2,:),'linewidth',2)
%set(gca,'fontsize',11,'linewidth',2)
%xlabel('intercept'); ylabel('slope')

% mc test

residuals=y-bestmodel; sestimate=std(residuals);

% need to make simulated dataset

count = 1:length(y); % an integer vector with an entry in order per each point in the data

% replacing all data seems to work better

%noreplace=round(length(y)/3); % replace 1/3 of the datapoints
noreplace=round(length(y)/1); % replace all of the datapoints

indexreplace = randsample(count,noreplace); % pick at random wich of the 1/3 of the points to replace.

Ynew=y; % initialize the simulated data just as the original data.

for i=1:length(indexreplace)
    Ynew(indexreplace(i))=bestmodel(indexreplace(i))+randn(1,1)*sestimate;
end

%plot(x,y,'ko',x,Ynew,'k.')

% full MC


for j=1:1000

    Ynew=y; % initialize the simulated data just as the original data.

    for i=1:length(indexreplace)
        Ynew(indexreplace(i))=bestmodel(indexreplace(i))+randn(1,1)*sestimate;
    end
    
    MCbetas(:,j)=inv(X'*X)*X'*Ynew;
    fit(:,j)=X*MCbetas(:,j);
    
end

MCslope=MCbetas(2,:); MCintercepts=MCbetas(1,:);

mean(MCslope)
std(MCslope)

figure(1)

subplot(211); 
[counts,centres]=hist(MCslope,10);
hist(MCslope,10)
hold on; plot([b(2) b(2)],[0 max(counts)],'k','linewidth',2)
plot([b(2)-seB(2) b(2)-seB(2)],[0 max(counts)],'k--','linewidth',2)
plot([b(2)+seB(2) b(2)+seB(2)],[0 max(counts)],'k--','linewidth',2)
axis([b(2)-seB(2)*1.2 b(2)+seB(2)*1.2 0 max(counts)*1.1])

subplot(212)
[counts,centres]=hist(MCintercepts,10);
hist(MCintercepts,10)
hold on; plot([b(1) b(1)],[0 max(counts)],'k','linewidth',2)
plot([b(1)-seB(1) b(1)-seB(1)],[0 max(counts)],'k--','linewidth',2)
plot([b(1)+seB(1) b(1)+seB(1)],[0 max(counts)],'k--','linewidth',2)
axis([b(1)-seB(1)*1.2 b(1)+seB(1)*1.2 0 max(counts)*1.1])


% plot data and swarm of lines

figure(2)
plot(x,bestmodel,x,y,'ko')
hold on
plot(x,fit)

% add the dashed lines

Fvalue=finv(0.95,P,N-P);
xplot=1.8:0.1:3.2; xplot=xplot';
XPLOT=[ones(size(xplot)) xplot];
model=XPLOT*b;

plot(x,y,'ko','markersize',4,'markerfacecolor','b')
set(gca,'linewidth',2,'fontsize',11)
hold on

plot(xplot,model,'k','linewidth',2)

for i=1:length(xplot)
	xh=[1; xplot(i)];
	upper(i)=xh'*b+s*sqrt(xh'*inv(X'*X)*xh)*sqrt(P*Fvalue);
	lower(i)=xh'*b-s*sqrt(xh'*inv(X'*X)*xh)*sqrt(P*Fvalue);
end

plot(xplot,upper,'k--')
plot(xplot,lower,'k--')