function [ R ] = cholupdatek( R , X , str)
%CHOLUPDATEK Rank-k Cholesky update helper function

    for i = 1:size(X,2)
        try
        R = cholupdate(R,X(:,i),str);
        catch
        end
    end
end

