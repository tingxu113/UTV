function Dy=diffy(D)
[m,n]=size(D);
Dy=zeros(m-1,m);
        for j=1:m-1
            Dy(j,j)=1;
            Dy(j,j+1)=-1;
        end
end
