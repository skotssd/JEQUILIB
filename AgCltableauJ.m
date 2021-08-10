
function [Agp,Clm,AgClLsR,Nap,NO3m,AgCl,AgCl2m,AgCl3m2,AgCl4m3,MASSERR]=AgCltableauJ(pH,pe,AgT,ClT,NaT,NO3T,flag1,flag2,database,acid)

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

T=[ClT; AgT; NaT; NO3T];

% end of tableau.  ------------------ % ----------------------------------------------

[KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES]=processtableau(Tableau,pH,pe);

[SPECIESCONCS,SPECIATIONNAMES,MASSERR,X]=returnspeciationRE(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2,database,acid);

% this will generate the outputs

for k=1:size(SPECIESCONCS,1)
      txt=[SPECIATIONNAMES(k,:),'=SPECIESCONCS(k);'];
      eval(txt)
end

end