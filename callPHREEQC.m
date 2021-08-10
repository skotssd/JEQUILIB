
% return speciation as a function of pH, FeT, PT

function [solutionspeciesconcs, speciesnames, solidconcs, solidnames]=callPHREEQC(totalnames,totalvector,pH,pe,T,minerals,speciesexport,database,acid,flag2)
temp_val=T; pe_val=pe; pH_val=pH; NOOFSOLIDS=size(minerals,1);

% Construct the text file to run PHREEQC from MATLAB-------------------
fileID=fopen('runphreeqc.txt','w');
fclose(fileID);
fileID=fopen('runphreeqc.txt','a');
% add stuff to database
%fprintf(fileID,'PHASES\n');

% add stuff to database
fprintf(fileID,['PHASES \n']);
fprintf(fileID,['pH_Fix\n']);
fprintf(fileID,['H+ = H+\n']);
fprintf(fileID,['-log_k           0.0\n']);
fprintf(fileID,['pe_Fix\n']);
fprintf(fileID,['e- = e-\n']);
fprintf(fileID,['-log_k           0.0\n']);

% define solution  -------------------------------------------------------
fprintf(fileID,'SOLUTION \n');
fprintf(fileID,['       pe      ' num2str(pe_val), '\n']); % the redox value
fprintf(fileID,['       pH      ' num2str(pH_val), '\n']); % the pH condition
fprintf(fileID,['       temp      ' num2str(temp_val), '\n']); % the temperature
fprintf(fileID,'-units mol/kgw\n'); % the unit of the input; usually mol/L is used
% put in the totals
for i=1:size(totalnames,1)
   totaltxt=[cell2mat(totalnames(i)),' ', num2str(totalvector(i)), ' \n'];
   fprintf(fileID,totaltxt);
end

% define equilibrium phases (solids, fix pH, pe) -------------------------
fprintf(fileID,'EQUILIBRIUM_PHASES \n'); 
for i=1:NOOFSOLIDS
   phasestxt=[cell2mat(minerals(i)), '   0.0   0\n'];
   fprintf(fileID,phasestxt);
end

pHfixline=['       pH_Fix -',num2str(pH_val),'          ',acid,' 10.0\n']; 
fprintf(fileID,pHfixline);
fprintf(fileID,'-force_equality true\n');
pefixline=['       pe_Fix ',num2str(-1*pe_val),'          O2\n']; 
fprintf(fileID,pefixline);
fprintf(fileID,'-force_equality true\n');

% define outputs of model ------------------------------------------------
fprintf(fileID,'SELECTED_OUTPUT\n');
fprintf(fileID,'-file selected.out\n');
fprintf(fileID,'-selected_out true\n');
fprintf(fileID,'-user_punch true\n');
fprintf(fileID,'-high_precision true\n');
fprintf(fileID,'-reset false\n');
fprintf(fileID,'-simulation false\n');
fprintf(fileID,'-state false\n');
fprintf(fileID,'-distance false\n');
fprintf(fileID,'-time false\n');
fprintf(fileID,'-step false\n');
fprintf(fileID,'-ph true\n');
fprintf(fileID,'-pe true\n');
fprintf(fileID,'-reaction false\n');
fprintf(fileID,'-temperature false\n');
fprintf(fileID,'-alkalinity false\n');
fprintf(fileID,'-ionic_strength false\n');
fprintf(fileID,'-water false\n');
fprintf(fileID,'-charge_balance false\n');
fprintf(fileID,'-percent_error false\n');

%fprintf(fileID,'KNOBS\n');
%fprintf(fileID,'-iterations 1500\n');

% species to export
speciestxt=['-molalities '];
for i=1:size(speciesexport,1)
speciestxt=[speciestxt, cell2mat(speciesexport(i)), ' '];
end
speciestxt=[speciestxt, ' \n'];
fprintf(fileID,speciestxt);
equilibriumphases=['-equilibrium_phases '];
for i=1:NOOFSOLIDS
    equilibriumphases=[equilibriumphases, cell2mat(minerals(i)), '  '];
end
fprintf(fileID,equilibriumphases);
fclose(fileID);

% run the model -----------------------------------------------------------
str=['!phreeqc runphreeqc.txt out.txt ', database];
%system('export LD_PRELOAD=/usr/lib/libstdc++.so.6')
%system('phreeqc runphreeqc.txt out.txt llnl.dat')
%!phreeqc runphreeqc.txt out.txt llnl.dat
if flag2==2; eval(str); end % output to the screen
if flag2==1; evalc(str); end % so no screen output

% import the data from the run and prepare to export----------------
fid = fopen('selected.out','rt');
hdr = strtrim(regexp(fgetl(fid),'\t','split')); hdr=hdr';
mat = cell2mat(textscan(fid,repmat('%f',1,numel(hdr))));
fclose(fid);

out_PHREEQC=mat';
[n,m]=size(out_PHREEQC); hdr=hdr(1:n-1); out_PHREEQC=out_PHREEQC(1:n-1,:); n=n-1;
selectedconcs=out_PHREEQC(:,2);
solutionspeciesconcs=selectedconcs(1:n-2*NOOFSOLIDS);
solidconcs=selectedconcs(n-2*NOOFSOLIDS+1:n);
speciesnames=(hdr(1:n-2*NOOFSOLIDS,1));
solidnames=(hdr(n-2*NOOFSOLIDS+1:n,1));

end