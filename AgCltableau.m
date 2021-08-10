
function [Agp,Clm,AgClLsR,Nap,NO3m,AgCl,AgCl2m,AgCl3m2,AgCl4m3,MASSERR]=AgCltableauJ(pH,pe,AgT,ClT,flag1,flag2,database,acid)

% input tableau.  change this part % ----------------------------------------------

Tableau=[...
%H        e        Cl        Ag      Na      NO3-   logK     phase     species1 
1         0        0         0       0       0       0            0             {'H'}
0         1        0         0       0       0       0            0             {'e '}
0         0        1         0       0       0       0            0             {'Cl-'}
0         0        0         1       0       0       0            0             {'Ag+'}
0         0        0         0       1       0       0            0             {'Na+'}
0         0        0         0       0       1       0            0             {'NO3-'}
0         0        1         1       0       0       9.7453       1             {'AgCl(s)'}
0         0        1         1       0       0       3.31         0             {'AgCl'}
0         0        2         1       0       0       5.24         0             {'AgCl2-'}
0         0        3         1       0       0       5.2          0             {'AgCl3-2'}
0         0        4         1       0       0       5.32         0             {'AgCl4-3'}
];

NaT=ClT; NO3T=AgT;

if flag1==3; T=[ClT; AgT; NaT; NO3T]; end; % input totals
if flag1==1; T=[ClT; AgT]; Tableau=redtableau(Tableau); Nap=NaT; NO3m=NO3T; end
if flag1==2; T=[ClT; AgT]; Tableau=redtableau(Tableau); Nap=NaT; NO3m=NO3T; end

% end of tableau.  ------------------ % ----------------------------------------------

[KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES]=processtableau(Tableau,pH,pe);

[SPECIESCONCS,SPECIATIONNAMES,MASSERR,X]=returnspeciationRE(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2,database,acid);

% this will generate the outputs

for k=1:size(SPECIESCONCS,1)
      txt=[SPECIATIONNAMES(k,:),'=SPECIESCONCS(k);'];
      eval(txt)
end

end

function redtableau=redtableau(Tableau)

[N,M]=size(Tableau); redtableau=Tableau;

% M is the number of components. first two have to keep.  H and e
c=0;
for i=3:M-2
    tst=cell2mat(Tableau(:,i)); tstsum=sum(tst);
    if tstsum==1
        c=c+1;
        rowindex(c) = find(tst==1);
        columnindex(c)=i;
    end
end

C=0;
for i=1:c
    redtableau(rowindex(i)-C,:)=[];
    redtableau(:,columnindex(i)-C,:)=[];
    C=C+1;
end

end