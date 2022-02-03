function  exercise6_L4_equilib

figure(1); clf; figure(2); clf;

%define the equilibrium sytem (The tableau)
[KSOLUTION,ASOLUTION,SOLUTIONNAMES] = get_equilib_defn;

% initial guess (just Fe becasue pH is fixed)
Feguess=[-5.5]; Lguess=-6; guess=[10.^Feguess; 10.^Lguess];

%set the pH range and total
pH=2:0.1:12; FeT=1e-4; LT=4e-5; T=[FeT; LT];

%now for each pH solve
for i=1:length(pH)
% adjust for fixed pH
[Ksolution,Asolution]=get_equilib_fixed_pH(KSOLUTION,ASOLUTION,pH(i));
% calculate species using NR
[X,F,J,C] = nl_massbalancerrnosolid_NR(guess,Asolution,Ksolution,T);
species_summary(:,i)=C; err(:,i)=F;
end

for i=1:size(species_summary,1)
txt=[SOLUTIONNAMES(i,:),'=species_summary(i,:);'];
eval(txt);
end



figure(1);
h=plot(pH,Fe./T(1),pH,FeOH./T(1),pH,FeOH2./T(1),pH,FeOH4./T(1),pH,FeL./T(1));
set(gca,"fontsize",12); set(h,'linewidth',2); set(gca,'linewidth',2);
h=xlabel('pH'); set(h,'fontsize',12); h=ylabel('fraction of total Fe'); set(h,'fontsize',12);
legend('Fe','FeOH','FeOH2','FeOH4','FeL','Location','northeast','Orientation','vertical');
axis([min(pH) max(pH) 0 1])

figure(2)
plot(pH,err(1,:),'bo',pH,err(2,:),'r.')


end

%------------------ SUBFUNCTIONS --------------------

function [KSOLUTION,ASOLUTION,SOLUTIONNAMES] = get_equilib_defn 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tableau=[...
%H	  FeIII    L    logK            species
 1	  0	       0     0              {'H'}
 0	  1	       0     0              {'Fe'}
 0    0        1      0             {'L'}
 -1	  0	       0     -14            {'OH'}
-1    1        0     -2.19          {'FeOH'}
-2    1        0     -5.67          {'FeOH2'}
-4    1        0     -21.6          {'FeOH4'}
0     1        1    17.55           {'FeL'}
1     0        1     13.7           {'HL'}
2     0        1      16.67          {'H2L'}
];

n=size(Tableau,2);
ASOLUTION=cell2mat(Tableau(:,1:n-2));
KSOLUTION=cell2mat(Tableau(:,n-1));
SOLUTIONNAMES=strvcat(Tableau(:,n));

end

% ----------- for fixed pH ----------------

function [Ksolution,Asolution]=get_equilib_fixed_pH(KSOLUTION,ASOLUTION,pH)

    [N,M]=size(ASOLUTION);
    Ksolution=KSOLUTION-ASOLUTION(:,1)*pH;
    Asolution=[ASOLUTION(:,2:M)];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,F,J,C] = nl_massbalancerrnosolid_NR(X,Asolution,Ksolution,T)

[Nc,Nx]=size(Asolution); %Xsolution=X(1:Nx);
criteria=1e-15;

for i=1:1000

logC=(Ksolution)+Asolution*log10(X); C=10.^(logC); % calc species
R=Asolution'*C-T;

% Evaluate the Jacobian 
   z=zeros(Nx,Nx); 
for j=1:Nx 
	for k=1:Nx 
		for i=1:Nc; z(j,k)=z(j,k)+Asolution(i,j)*Asolution(i,k)*C(i)/X(k); end
   	end
end

J = z;

deltaX=z\(-1*R);
one_over_del=max([1, -1*deltaX'./(0.5*X')]);
del=1/one_over_del; X=X+del*deltaX;
    
tst=sum(abs(R));
if tst<=criteria; break; end

end

F=[R]; 

end
