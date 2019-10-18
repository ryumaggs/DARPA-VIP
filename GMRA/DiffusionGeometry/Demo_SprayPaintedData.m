%% Demo with some simple 2-dimensional "spray-painted" data sets
% New data sets may be created from any software that saves to bmp's (other formats may be easily added)

for k = 1:7
    ImageName = sprintf('DataSetFromImageEx_%.2d.bmp',k);
    GraphEigenFcns_FromImage
end