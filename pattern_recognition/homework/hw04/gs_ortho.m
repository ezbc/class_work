
function Q = Gram_Schmidt(V)
%Gram_Schmidt process
    sz = size(V);
    vprime = zeros(sz);
    Q = zeros(sz);
    Q(:,1) = V(:,1)/norm(V(:,1));
    for i=2:sz(2)
        vprime(:,i) = zeros(sz(1),1);
        % sum the projections from all previous vectors
        for j=1:i-1
            vprime(:,i) = vprime(:,i) + (Q(:,j)'*V(:,i))*Q(:,j);
        end
        vprime(:,i) = V(:,i) - vprime(:,i);
        if any(vprime(:,i))
            Q(:,i) = vprime(:,i)/norm(vprime(:,i));
        else 
            Q(:,i) = vprime(:,i);
        end
    end
end