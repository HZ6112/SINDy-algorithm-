function [Phi, Lambda, b] = DMD(X,Xprime,r)
[U,Sigma,V] = svd(X,'econ'); 
Ur = U(:,1:r);
Sigmar = Sigma(1:r,1:r);
Vr = V(:,1:r);
Atilde = Ur'*Xprime*Vr/Sigmar; [W,Lambda] = eig(Atilde);
Phi = Xprime*(Vr/Sigmar)*W; 
alpha1 = Sigmar*Vr(1,:)'; 
b = (W*Lambda)\alpha1;
end
X = VORTALL;
[Phi, Lambda, b] = DMD(X(:,1:end-1),X(:,2:end),21);