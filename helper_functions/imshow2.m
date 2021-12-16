function h = imshow2(varargin)
    varargin = [varargin, 'InitialMagnification', 'fit'];
    hh = imshow(varargin{:}); 
      
    if (nargout > 0)
    % Only return handle if caller requested it.
    h = hh;
    end
  
end