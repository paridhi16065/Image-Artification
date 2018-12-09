function [Xtilde]=irls(R,X,z)
    [tNc,Nij]=size(R);
    I= 5; %iterations
    Xk=X; 
    r=0.8;
    unsampled_pixs=double(~(sum(R,2)>0)); 
    for k=1:I
        
        A=unsampled_pixs;
        B=Xk.*unsampled_pixs;
        for i=1:Nij
            w=sum((Xk(logical(R(:,i)))-z(:,i)).^2 + 1e-10).^((r-2)/2);
            A=A+w*R(:,i);
            temp=R(:,i);
            temp(logical(temp))=z(:,i);
            B=B+w*temp;
        end
        Xk=(1./(A+1e-10)).*B;
    end
    Xtilde=Xk;
end
