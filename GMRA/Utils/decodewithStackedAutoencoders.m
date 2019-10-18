function X = decodewithStackedAutoencoders( feat, autoenc )

J = length(autoenc);

for j = J:-1:1
    if j<J
        X{j} = decode( autoenc{j}, X{j+1} );
    else
        X{j} = decode( autoenc{j}, feat{j} );
    end
end

return