
function [Eup3,PO4m3,EuPO4LsR,EuLOHR3LsR,MASSERR]=...
EuPsaltstableau(pH,pe,EuT,PT,CT,MgT,CaT,ST,NO3T,NaT,flag1,flag2,database,acid)

% input tableau.  change this part % ----------------------------------------------

pKa1=-log10(7.5e-3); pKa2=-log10(6.2e-8); pKa3=-log10(2.14e-13); khco3=10.3; kh2co3=6.3;

Tableau=[...
%H        e        Eu        PO4     Na      NO3-  CO3   Mg   Ca   SO4    logK                 phase            species1 
1         0        0         0       0       0      0    0    0     0      0                     0             {'H'}
0         1        0         0       0       0      0    0    0     0      0                     0             {'e'}
0         0        1         0       0       0      0    0    0     0      0                     0             {'Eu+3'}
0         0        0         1       0       0      0    0    0     0      0                     0             {'PO4-3'}
0         0        0         0       1       0      0    0    0     0      0                     0             {'Na+'}
0         0        0         0       0       1      0    0    0     0      0                     0             {'NO3-'}
0         0        0         0       0       0      1    0    0     0      0                     0             {'CO3-2'}
0         0        0         0       0       0      0    1    0     0      0                     0             {'Mg+2'}
0         0        0         0       0       0      0    0    1     0      0                     0             {'Ca+2'}
0         0        0         0       0       0      0    0    0     1      0                     0             {'SO4-2'}
0         0        1         1       0       0      0    0    0     0      24                    1             {'EuPO4(s)'}
1         0        0         1       0       0      0    0    0     0      pKa3                  0             {'HPO4-2'}
2         0        0         1       0       0      0    0    0     0      pKa3+pKa2             0             {'H2PO4-'}
3         0        0         1       0       0      0    0    0     0      pKa3+pKa2+pKa1        0             {'H3PO4'}
-2        0        2         0       0       0      0    0    0     0      -6.9182               0             {'Eu2(OH)2+4'}
-1        0        1         0       0       0      0    0    0     0      -7.9075               0             {'EuOH+2'}
-2        0        1         0       0       0      0    0    0     0      -16.337               0             {'EuO+'}
-3        0        1         0       0       0      0    0    0     0      -15.3481              1             {'Eu(OH)3(s)'}
1         0        0         0       0       0      1    0    0     0      khco3                 0             {'HCO3-'}
2         0        0         0       0       0      1    0    0     0      kh2co3                0             {'CO2'}
1         0        1         0       0       0      1    0    0     0      1.6258+khco3          0             {'EuHCO3+2'}
0         0        1         0       0       0      3    0    0     0      -16.8155+(3*khco3)    0             {'Eu(CO3)3-3'}
0         0        1         0       0       0      2    0    0     0      -8.3993+(2*khco3)     0             {'Eu(CO3)2-'}
0         0        1         0       0       0      1    0    0     0      -2.4057+khco3         0             {'EuCO3+'}
-1        0        1         0       0       0      1    0    0     0      -8.4941+khco3         1             {'EuOHCO3(s)'}
-1        0        1         0       0       0      2    0    0     0      -15.176+(2*khco3)     0             {'EuOH(CO3)2-2'}
-2        0        1         0       0       0      1    0    0     0      -17.8462+khco3        0             {'Eu(OH)2CO3-'}
0         0        0         1       0       0      0    1    0     0      -5.7328+(1/pKa3)      0             {'MgPO4-'}
1         0        0         1       0       0      0    1    0     0      2.9100+(1/pKa3)       0             {'MgHPO4'}
2         0        0         1       0       0      0    1    0     0      1.6600+(1/pKa3)       0             {'MgH2PO4+'}
0         0        0         0       0       0      1    1    0     0      -7.3499+khco3         1             {'Magnesite'}
1         0        0         0       0       0      1    1    0     0      1.0357+khco3          0             {'MgHCO3+'}
%-4        0        0         0       0       0      0    4    0     0      -39.75                0             {'Mg(OH)4-2'}
0         0        0         1       0       0      0    0    1     0      -5.8618+(1/pKa3)      0             {'CaPO4-'}
1         0        0         1       0       0      0    0    1     0      2.7400+(1/pKa3)       0             {'CaHPO4'}
2         0        0         1       0       0      0    0    1     0      1.4000+(1/pKa3)       0             {'CaH2PO4+'}
0         0        0         0       0       0      1    0    1     0      -7.0017+khco3         1             {'Calcite'}
1         0        0         0       0       0      1    0    1     0      1.0467+khco3          0             {'CaHCO3+'}
-1        0        0         0       0       0      0    0    1     0      -12.85                0             {'CaOH+'}
1         0        0         0       0       1      0    0    0     0      -1.3025               0             {'HNO3'}
0         0        1         0       0       1      0    0    0     0      0.8745                0             {'EuNO3+2'}
0         0        0         0       0       1      0    0    1     0      0.7000                0             {'CaNO3+'}
1         0        0         0       0       0      0    0    0     1      1.9791                0             {'HSO4-'}
0         0        1         0       0       0      0    0    0     1      3.6430                0             {'EuSO4+'}
0         0        1         0       0       0      0    0    0     2      5.4693                0             {'Eu(SO4)2-'}
0         0        0         0       0       0      0    1    0     1      2.4117                0             {'MgSO4'}
0         0        0         0       0       0      0    0    1     1      4.50                  1             {'Gypsum'}
];

T=[EuT; PT; NaT; NO3T; CT; MgT; CaT; ST]; 

% end of tableau.  ------------------ % ----------------------------------------------

[KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES]=processtableau(Tableau,pH,pe);

[SPECIESCONCS,SPECIATIONNAMES,MASSERR,X]=returnspeciationRE(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2,database,acid);

% this will generate the outputs

for k=1:size(SPECIESCONCS,1)
      txt=[SPECIATIONNAMES(k,:),'=SPECIESCONCS(k);'];
      eval(txt)
end

end