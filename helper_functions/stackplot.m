function stackplot(m)
%STACKPLOT Plot function replacement to show a stack of time traces.
%   The y axis is reversed, showing each of the time traces centered in its
%   index number. Only works for single input plot calls. 
%   
%   2014 Vicente Parot
%   Cohen Lab - Harvard University
%
plot(bsxfun(@plus,(1:size(m,2)),-1.3*m/(max(m(:))-min(m(:)))));
set(gca,'YDir','reverse')
