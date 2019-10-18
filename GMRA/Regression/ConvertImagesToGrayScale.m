function grayImages = ConvertImagesToGrayScale( RGBImages, imSize )

grayImages = zeros(prod(imSize),size(RGBImages,2),class(RGBImages));

for k = 1:size(RGBImages,2)
    grayImages(:,k) = reshape(rgb2gray(reshape(RGBImages(:,k),[imSize(1),imSize(2),3])),[imSize(1)*imSize(2),1]);
end

return