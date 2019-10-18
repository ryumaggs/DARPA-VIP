function MakeFigLargeAndNice( figHandle )

GOLDEN_RATIO = 1.618033988749894848204586834365638;

ss      = get(0,'screensize'); %The screen size
width   = ss(3);
height  = ss(4);

vert    = height/GOLDEN_RATIO;
horz    = width/GOLDEN_RATIO;

set(figHandle,'Position',[(width/2)-horz/2, (height/2)-vert/2, horz, vert]);

return