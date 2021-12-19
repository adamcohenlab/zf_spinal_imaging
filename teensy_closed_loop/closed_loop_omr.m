% Matlab script to communicate with the Teensy 3.6 when it is in closed
% loop OMR mode
%
% Calling closed_loop_omr without arguments opens the GUI and initializes
% the threshold of burst detectioin to 100
% 
% The GUI calls the subfunctions by passing arguments.
%
% Last modified 06/11/17.

function [] = closed_loop_omr(varargin)
global fig
if (nargin==0)                          % initialize
    
    fig = openfig('closed_loop_omr.fig','reuse');
    
    initializeteensy('COM3')            % replace COM3 with the name of
                                        % the USB port to which the Teensy
                                        % is connected

    zero;  
    
else                                    % feval switchyard
    
    if (nargout)
        [varargout{1:nargout}] = feval(varargin{:}); %#ok<NASGU>
    else
        feval(varargin{:});
    end
    
end


% -------------------------------------------------------------------------

function [] = initializeteensy(portName)
global teensy

teensy = serial(portName,'BaudRate',115200,'ByteOrder','littleEndian');
fopen(teensy);



% -------------------------------------------------------------------------

function [] = upload()
global teensy
global fig

thresh = get(findobj(fig, 'Tag','thresh'),'Value');

out = thresh;
out = typecast(single(out),'uint8'); 
fwrite(teensy,out);

% -------------------------------------------------------------------------

function [] = zero()
global fig

set(findobj(fig, 'Tag','thresh'),'Value',340);
getnumbers
upload

% -------------------------------------------------------------------------
function [] = getnumbers()
global fig
thresh = get(findobj(fig, 'Tag','thresh'),'Value');

set(findobj(fig, 'Tag','thresh_Num'),'String',num2str(thresh,'%0.2f'))



% -------------------------------------------------------------------------

function [] = closedcfigure %#ok<*DEFNU>
global teensy

delete(findobj('Tag','closed_loop_omr_figure'))
fclose(teensy);
delete(teensy)
clear teensy


