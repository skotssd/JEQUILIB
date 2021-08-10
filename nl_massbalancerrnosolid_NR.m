

function [X,F,J,RSI,C] = nl_massbalancerrnosolid_NR(X,Asolution,Ksolution,Asolid,Ksolid,T,TYPX)

[Nc,Nx]=size(Asolution); %Xsolution=X(1:Nx);
criteria=1e-16;

for i=1:1000

logC=(Ksolution)+Asolution*log10(X); C=10.^(logC); % calc species
R=Asolution'*C-T; 

% Evaluate the Jacobian 
   z=zeros(Nx,Nx); 
for j=1:Nx; 
	for k=1:Nx; 
		for i=1:Nc; z(j,k)=z(j,k)+Asolution(i,j)*Asolution(i,k)*C(i)/X(k); end
   	end
end

    %[jac,err] = jacobianest(fun,X);
    %J=jac;
    J=z;
    Jstar = J./TYPX';
    %pause % scale Jacobian
    dx_star = -Jstar\(R); % faster than the alternatives it seems.  if warnings set to off.
    %dx_star =pinv(z)*(-1*R); % no scaling
    %dx_star =pinv(Jstar)*(-1*R);
    %deltaX=z\(-1*R);
    %deltaX = pinv(Jstar)*(-1*R);
    %dx_star = pinv(Jstar)*(-1*R);
    dx = dx_star.*TYPX; %reverse the scaling
    deltaX=dx;
    % rescale x

%deltaX=z\(-1*R);
%deltaX = pinv(z)*(-1*R);
one_over_del=max([1, -1*deltaX'./(0.5*X')]);
del=1/one_over_del; X=X+del*deltaX;
    
tst=sum(abs(R));
if tst<=criteria; break; end

end

Q=Asolid*log10(X); SI=(Q+Ksolid);
RSI=ones(size(SI))-SI;

F=[R]; 

end