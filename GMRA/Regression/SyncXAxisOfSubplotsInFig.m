function SyncXAxisOfSubplotsInFig( figHandle, axesHandles )

if nargin<2 || isempty(axesHandles)
    axesHandles = get(figHandle,'Children');
end

for k = 1:length(axesHandles)
   classname = class(axesHandles(k));
   if ~isempty(strfind(classname,'Axes'))
       XLim(k,:) = get(axesHandles(k),'XLim');
   end
end

newXLim = [min(XLim(:,1)),max(XLim(:,2))];

for k = 1:length(axesHandles)
   classname = class(axesHandles(k));
   if ~isempty(strfind(classname,'Axes'))
       set(axesHandles(k),'XLim',newXLim);
   end
end

return