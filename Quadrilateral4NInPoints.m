function [weig,posgp,shapef,dershapef] = Quadrilateral4NInPoints

weig = [1 1 1 1];
posgp=1/sqrt(3)*[-1 -1; 1 -1; 1 1; -1 1]';
ndim=2; nnodeE=4;
ngaus=length(weig);
shapef=zeros(ngaus,nnodeE);
dershapef=zeros(ndim,nnodeE,ngaus);
for j=1:length(weig)
    xi=posgp(1,j);
    eta=posgp(2,j);
    [Ne, Bex]=Quadrilateral4N(xi,eta);
    shapef(j,:)=Ne;
    dershapef(:,:,j)=Bex;
end
end

function [N, B]=Quadrilateral4N(xi,eta)
N=0.25*[(1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)];

B=0.25*[-(1-eta) (1-eta) (1+eta) -(1+eta);
           -(1-xi) -(1+xi) (1+xi) (1-xi)];
end
