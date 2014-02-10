function movi16(seq,type,rate,minmax)

% function movi16(seq,type,rate,minmax)
% 
% Display movies of a 3-D uint16, uint8 or double dataset seq at max frame
% rate 'rate' fr/s (the actual rate can be lower due limited processor speed).
% seq is a 3-D uint16 or uint8 array of frames seq(:,:,k)
% type = 'hsv','jet','cool' etc define any of the Matlab colour maps.
% minmax is a 2-element vector defining the min and max grey-levels for display
% (defaults to min and max values in seq).
%
% Three orthogonal views of the dataset are displayed.
% GUIs on the figure allow the following controls (from top left):
%   Brightness and Contrast sliders (allow contrast to be enhanced so that
%     the full resolution of uint16 data can be displayed if required);
%   Selection (in green) of which view is the one to be incremented through the slices;
%   Buttons to control the slice position of the chosen view;
%   Exit button;
%   Slider to control the slice position;
%   Cross-hairs on/off toggle button.
%
% The routine uses a continuous loop for display, rather than being GUI button
% driven, so that it is not necessary to save a copy of the 3-D dataset seq
% in gcf.UserData, which could use up a lot of memory.
%
% Nick Kingsbury, Cambridge University, March 2000.

global ind di brk newsel 

if nargin < 2, type = []; end
if isempty(type), type = 'bone'; end

% Open figure and set up colour map for uint8 images.
figure(gcf)
colormap(eval([type '(256)']));

if nargin < 3, rate = []; end
if isempty(rate), rate = 10; end  % Default max frame rate.

% Find the initial range for image display (using min contrast).
if nargin < 4, minmax = []; end
if isempty(minmax),
	maxseq = double(max(seq(:)));
   minseq = double(min(seq(:)));
else
   maxseq = minmax(2);
   minseq = minmax(1);
end

fprintf('Fig %d: max = %d, min = %d\n',gcf,maxseq,minseq);

% Generate slice position control buttons with labels.
label = strvcat(' < | ','  <  ','  O  ','  >  ',' | > ',' Exit');
for u = 1:6,
   ui(u) = uicontrol(gcf,'str',label(u,:),'pos',[40*u 0 40 20]);
end

% Generate sliders for contrast, brightness and slice position.
ui(7) = uicontrol(gcf,'style','slider','value',0,'pos',[12 70 16 120]);
ui(8) = uicontrol(gcf,'style','slider','value',0.5,'pos',[12 210 16 120]);
ui(9) = uicontrol(gcf,'style','slider','pos',[280 0 160 20]);

% Generate buttons for image selection and cross-hairs on/off.
ui(10) = uicontrol(gcf,'pos',[0 30 30 30]);
ui(11) = uicontrol(gcf,'pos',[30 30 30 30]);
ui(12) = uicontrol(gcf,'pos',[0 0 30 30]);
ui(13) = uicontrol(gcf,'style','toggle','pos',[480 0 40 20],'str','X-hairs','value',1);

% Add labels.
ui(14) = uicontrol(gcf,'style','text','pos',[0 190 40 14],'str','Contrast');
ui(15) = uicontrol(gcf,'style','text','pos',[0 330 50 14],'str','Brightness');
ui(16) = uicontrol(gcf,'style','text','pos',[60 30 40 30],'str','Select Image');

% Define 'call' functions for position controls and image selection.
set(ui(1),'call','global ind newsel;di=0; ind(newsel)=ind(newsel)-1;')
set(ui(2),'call','global di;di=-1;')
set(ui(3),'call','global di;di=0;')
set(ui(4),'call','global di;di=1;')
set(ui(5),'call','global ind newsel;di=0; ind(newsel)=ind(newsel)+1;')
set(ui(6),'call','global brk;brk=1;')
set(ui(10),'call','global newsel;newsel=3;')
set(ui(11),'call','global newsel;newsel=2;')
set(ui(12),'call','global newsel;newsel=1;')

% Determine image parameters and plot the 3 image views.
contrast = get(ui(7),'value');
brightness = get(ui(8),'value');
imscale = exp(4*contrast) * 255/(maxseq - minseq + 1);
imoffset = 128 - imscale * (minseq + (1 - brightness)*(maxseq - minseq));
ss = size(seq);
ind = round(ss/2); % Initial slice position indices.
% Grey areas between subimages.
gap = 160;
w1 = uint8(gap*ones(ss(1),8)); 
w2 = uint8(gap*ones(8,ss(2)+8+ss(3)));
w3 = uint8(gap*ones(ss(3),8+ss(3)));
% Calculate views separately, using uint8 for speed.
im1 = uint8(double(seq(:,:,ind(3)))*imscale+imoffset);
im2 = uint8(double(squeeze(seq(:,ind(2),:)))*imscale+imoffset);
im3 = uint8(double(squeeze(seq(ind(1),:,:)).')*imscale+imoffset);
imh=image([im1 w1 im2; w2; im3 w3],'EraseMode','none');
axis image;axis off
set(gca,'position',[0.1 0.1 .85 .85]);

% Plot the cross-hairs, showing relative slice positions.
lh = line([[1;1]*[ind(2) ss(2)+8+ind(3)] [0 0;ss(2)+8+ss(3) ss(2)+8]], ...
   [[0 0;ss(1)+8+ss(3) ss(1)+8] [1;1]*[ind(1) ss(1)+8+ind(3)]], ...
   'color','r','EraseMode','none','linestyle',':');

dtime = 1/rate; % Frame period (sec).
t1 = clock;
t2 = dtime;
di = 0;

% Define selected view and colours for slice position control.
newsel = 3; sel = newsel;
colours = ['r';'r';'g';'r';'r']; 
newind = ind;
brk = 0;
xhairs = get(ui(13),'value'); % Turn cross-hairs on/off.
prev = [ind sel imscale imoffset xhairs]*0;

% This is a continuous loop until the 'Exit' button is pressed.
while 1,
   
   % Update slice position if slider has been used.
   if any(abs(newind - ind) > 1), ind = newind; end
   
   % Update slice position of selected view and slider.
   sel = newsel;
   if (ind(sel)+di) > ss(sel), di = -abs(di); % Code for continuous up/down sweeps.
   elseif (ind(sel)+di) < 1, di = abs(di);
   end
   ind(sel) = max(min(ind(sel) + di,ss(sel)),1);
   set(ui(9),'value',ind(sel)/ss(sel));
   set(ui(3),'str',sprintf('%d',ind(sel)));
   
   % Update colours and slice positions on 'Select Image' buttons.
   set(ui(10),'BackgroundColor',colours(sel,:),'str',sprintf('%d',ind(3)));
   set(ui(11),'BackgroundColor',colours(sel+1,:),'str',sprintf('%d',ind(2)));
   set(ui(12),'BackgroundColor',colours(sel+2,:),'str',sprintf('%d',ind(1)));
   
   % Update contrast and brightness parameters.
   contrast = get(ui(7),'value');
   brightness = get(ui(8),'value');
   imscale = exp(4*contrast) * 255/(maxseq - minseq + 1);
   imoffset = 128 - imscale * (minseq + (1 - brightness)*(maxseq - minseq));
   xhairs = get(ui(13),'value');
   
   % Update appropriate views if any parameters have changed.
   change = [ind sel imscale imoffset xhairs];
   ch = (change ~= prev); 
   if any(ch),     
      % Only update one view if possible, for max speed.
      if any(ch & [0 0 1 0 1 1 0]),
         im1 = uint8(double(seq(:,:,ind(3)))*imscale+imoffset); end
      if any(ch & [0 1 0 0 1 1 0]),
         im2 = uint8(double(squeeze(seq(:,ind(2),:)))*imscale+imoffset); end
      if any(ch & [1 0 0 0 1 1 0]),
         im3 = uint8(double(squeeze(seq(ind(1),:,:)).')*imscale+imoffset); end
      set(imh,'CData',[im1 w1 im2; w2; im3 w3]);
      % Image overwrites cross-hairs, so rewrite these if required.
      if (xhairs),
         r=rand(4,1)*1e-6; % This ensures that all cross hairs are changed slightly and hence updated.
         set(lh(1),'Xdata',[1;1]*ind(2)+r(1));
         set(lh(2),'Xdata',[1;1]*(ss(2)+8+ind(3))+r(2));        
         set(lh(3),'Ydata',[1;1]*ind(1)+r(3));
         set(lh(4),'Ydata',[1;1]*(ss(1)+8+ind(3))+r(4));
      end  
   else
      pause(dtime); % This frees the computer for other tasks when no changes are required.
   end
   prev = change;
   drawnow
   newind(sel) = round(get(ui(9),'value')*ss(sel));
   % Wait for elapsed time to expire before next frame.
   while etime(clock,t1) < dtime, dummy = 1; end
   t1 = clock;
   if brk, break; end % Exit from loop and terminate display.
end

% Remove GUI buttons.
delete(ui);

% Display position of slices and axes on final display.
xlabel(['Slices:' sprintf('  %d',ind)])
axis on

return



set(lh(1),'Xdata',[[1;1]*[ind(2) ss(2)+8+ind(3)] [0;ss(2)+8+ss(3)]*[1 1]], ...
   'Ydata',[[0;ss(1)+8+ss(3)]*[1 1] [1;1]*[ind(1) ss(1)+8+ind(3)]]);


lh(1) = line([1;1]*ind(2),[0;ss(1)+8+ss(3)],'color','r');
lh(2) = line([1;1]*(ss(2)+8+ind(3)),[0;ss(1)+8+ss(3)],'color','r');
lh(3) = line([0;ss(2)+8+ss(3)],[1;1]*ind(1),'color','r');
lh(4) = line([0;ss(2)+8+ss(3)],[1;1]*(ss(1)+8+ind(3)),'color','r');

        set(lh(1),'Xdata',[1;1]*ind(2),'EraseMode','none');
        set(lh(2),'Xdata',[1;1]*(ss(2)+8+ind(3)),'EraseMode','none');        
        set(lh(3),'Ydata',[1;1]*ind(1),'EraseMode','none');
        set(lh(4),'Ydata',[1;1]*(ss(1)+8+ind(3)),'EraseMode','none');        
        
        set(lh,'Xdata',[[1;1]*[ind(2) ss(2)+8+ind(3)] [0;ss(2)+8+ss(3)]*[1 1]], ...
           'Ydata',[[0;ss(1)+8+ss(3)]*[1 1] [1;1]*[ind(1) ss(1)+8+ind(3)]]);

