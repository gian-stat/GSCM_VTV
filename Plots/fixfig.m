%% fixfig
% FIXFIG(FIGNUM,VERBOSE)
%
% FIXFIG sets the font size and weight for all the items on a plot so that
%  they can be read more easily and used in presentations.
%
% Combined the original fixfig with paper settings 
%% Introduction
function fixfig(fignum,verbose)
% FIXFIG(FIGNUM,VERBOSE) fixes up a Matlab figure by increasing the font
%  sizes, making lines thicker, etc., so that the result is suitable for
%  presentation on a screen (e.g., Powerpoint). The font size for all text
%  is increased, lines are made thicker, markers are made larger,
%  background is set to white, etc. The values are set as constants near
%  the top of the file. You may want to change the sizes and fonts to suit
%  your system.
%
% FIGNUM    The number of the figure to modify. Default is the current
%           figure.
%
% VERBOSE   Level of status messages [0:2]. 0=none, 2=all. Default is 2.
%
% Example:
%
% Create a figure:
% 
% >> x = [1 2 3 4 5];
% >> a = rand(5,1);
% >> b = rand(5,1);
% >> plot(x,a,'.-r',x,b,'s-k');
% >> title('Random'); xlabel('x'); ylabel('a b');
% >> text(3,0.4,'The Point'); legend('Line1','Line2');
%
% The figure is plotted, with the default MATLAB fonts and sizes, which are
% quite small and impossible to read if used in a screen presentation.
% Now type:
%
% >> fixfig
%
% and the figure is transformed with larger fonts, bolder lines, etc.
%
%
% M. A. Hopcroft
%  mhopeng at ml1 dot net
%
% MH Apr2010
% v1.11  keep legend backgrounds opaque (thanks to SJ for feedback)
%        show axes tag
%        specify axes instead of using gca
%

% MH MAY2009
% v1.01 fix typo at line 96 "('mkr')"
% v1.0  added text and line objects
%       cleaned up for public release
%     
% MH MAR205
% v0.9 script for personal use
%

%% Set constants: line width, marker size, font

% How thick do we like our lines?
lnwidth = 2;
% How big do we like our markers?
mksize = 6;
mksizept = mksize * 3; % points ('.') are always smaller
% Which font do we like in our figures?
myfont = 'Times New Roman';
% How big do we like our fonts?
fontsizefact = 2; % font sizes are scaled by this value. 2-3 is typical.
% Font weight for axis and legends
fontwtaxlb = 'normal'; % {normal} | bold | light | demi
% Font weight for title
fontwttitl = 'normal'; % {normal} | bold | light | demi
% Font weight for xlabel
fontwtxlab = 'normal'; % {normal} | bold | light | demi
% Font weight for ylabel
fontwtylab = 'normal'; % {normal} | bold | light | demi
% Font weight for zlabel
fontwtzlab = 'normal'; % {normal} | bold | light | demi

% Note: set title and axis fonts near end of file
if nargin < 1
    fignum=gcf;
end
if nargin < 2
    verbose=2;
end


%% Modify the figure

% make the background white
set(fignum,'color','w');

% identify the subplots
a=get(fignum,'Children')';

% if verbose>=1, fprintf(1,'fixfig: Fixing Figure %d, with %d axes.\n',fignum,length(a)); end
% if verbose>=2, fprintf(1,'fixfig: Font: %s, Size Factor: %g.\n',myfont,fontsizefact); end

k=0;
for i=fliplr(a)
    k=k+1;
    %% Set the linewidths and marker sizes
    
    % find all the lines on this axis
    dataline = findobj(i,'Type','line');
%     if verbose>=2
%         fprintf(1,'fixfig: Axis %d has %d lines', k,length(dataline));
%         if ~isempty(get(i,'Tag'))
%             fprintf(1,' (%s).\n',get(i,'Tag'));
%         else
%             fprintf(1,'.\n');
%         end    
%     end

    % cycle through the lines on this axis
    for j = dataline'
        % set the linewidth
        set(j,'LineWidth',lnwidth);
        % set the marker size
        mkr = get(j,'Marker');
        if ~strcmp(mkr,'none')
            mkrsz = get(j,'MarkerSize');
            if mkrsz <= mksize, set(j,'MarkerSize',mksize); end
            if strcmp(mkr,'.')
                if mkrsz <= mksizept, set(j,'MarkerSize',mksizept); end
            end
        end
    end
    
    %% Set Font sizes
    
    % if there is text on the plot, find it and enlarge it
    datatext = findobj(i,'Type','text');
    for p = datatext'
        % set the font
        set(p,'FontSize',fontsizefact*7,'FontWeight','bold','FontName',myfont);
    end
    
    % axes(i); % specify the subplot
    set(fignum,'CurrentAxes',i);
    
    % axis values and legend text
    set(i,'FontSize',fontsizefact*8,'FontWeight',fontwtaxlb,'FontName',myfont)

    % title
    set(get(i,'Title'),'FontSize',fontsizefact*10,'FontWeight',fontwttitl,'FontName',myfont)

    % x,y,z axis label
    set(get(i,'XLabel'),'FontSize',fontsizefact*9,'FontWeight',fontwtxlab,'FontName',myfont)    
    set(get(i,'YLabel'),'FontSize',fontsizefact*9,'FontWeight',fontwtylab,'FontName',myfont)
    set(get(i,'ZLabel'),'FontSize',fontsizefact*9,'FontWeight',fontwtzlab,'FontName',myfont)

    % make background fill transparent so that legends, text, etc can be seen
    if ~strcmpi(get(i,'Tag'),'legend') % but keep legend opaque
        set(i,'Color','none')
    end
end
grid on
set(gcf,'Units','Pixels')
set(gcf,'Position',[100 100 830 580])
shg; pause(0.2)
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');   
set(gcf,'PaperPosition',[0 0 8.8 6.3]);  
set(gcf,'PaperSize',[8.5 6.0]);          
set(gcf,'PaperPosition',[0 0 8.5 6.0]);
set(gcf,'color', 'white');
set(gca,'LooseInset',get(gca,'TightInset'));
% legend(gca,'Location','SouthWest')
return