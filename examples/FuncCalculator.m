
classdef FuncCalculator < handle
properties
func, numEntries, xvals, yvals;
end
methods
  function obj=FuncCalculator(func, maxVals)
    obj.func = func;
    obj.numEntries = 0;
    if(nargin>1)
      maxvals = maxVals;
    else
      maxvals=100;
    end
    obj.xvals=zeros(maxvals,1);
    obj.yvals=zeros(maxvals,1);
  end
  function yVal=getVal(self,x)
  % call the function only if we don't have a saved value
  i = find(self.xvals==x);
  if(isempty(i) || i<=0)
    self.numEntries = self.numEntries + 1;
    i = self.numEntries;
    self.xvals(i) = x;
    yVal = self.func(x);
    self.yvals(i) = yVal;
  else
    yVal = self.yvals(i);
  end
  end
end
end