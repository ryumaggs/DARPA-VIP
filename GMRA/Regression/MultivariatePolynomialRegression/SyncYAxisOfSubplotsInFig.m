function SyncYAxisOfSubplotsInFig( figHandle, axesHandles )

if nargin<2 || isempty(axesHandles)
    axesHandles = get(figHandle,'Children');
end

for k = 1:length(axesHandles)
   classname = class(axesHandles(k));
   if ~isempty(strfind(classname,'Axes'))
       YLim(k,:) = get(axesHandles(k),'YLim');
   end
end

newYLim = [min(YLim(:,1)),max(YLim(:,2))];

for k = 1:length(axesHandles)
   classname = class(axesHandles(k));
   if ~isempty(strfind(classname,'Axes'))
       set(axesHandles(k),'YLim',newYLim);
   end
end

return