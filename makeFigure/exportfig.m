function exportfig(filename, format, m, type )
% EXPORTFIG - Export figure for articles, slides, etc.
% exportfig(filename,format, m)
% exportfig(filename,format, m, type )
%
% filename : filename without format name.
% format : image format including .eps, .pdf, .tiff, etc.
% m : size or aspect ratio of image.
% type: export format; 'TeX', 'TeX*', 'Slides' and 'Poster'.
%
% exsample
% exportfig('./fig', 'eps', [800 600]);
% exportfig('./fig', 'pdf', [4 3],'TeX');
%
% See also PRINT, SET, GET.

if isvector(m)
    if numel(m) ~= 2
        error('m must be [1 2] array');
    end
end

if nargin  < 4
    type = 'Manual';
end

% get figure handles
fig = gcf;
ax = findobj(fig,'Type','Axes');

if strcmpi(type,'Manual')
%         manual size
    set(fig,'PaperUnits','points');
    resizefig(fig,m);
    fsave(filename,format);
elseif  strcmpi(type,'TeX' )
%        export for TeX documents
    ratio = m(1)/m(2);
    w = 800;
    h = w/ratio;
    sfig(fig,ax,filename,format,[w h],'Times',24,30);
elseif  strcmpi(type,'TeX*' )
%        export for TeX documents as figure*
    ratio = m(1)/m(2);
    w = 800;
    h = w/ratio;
    sfig(fig,ax,filename,format,[w h],'Times',12,15);
elseif  strcmpi(type,'Slides');
%        export for Slides
    ratio = m(1)/m(2);
    h = 768;
    w = h*ratio;
    sfig(fig,ax,filename,format,[w h],'Helvetica',24,36);
elseif  strcmpi(type,'Poster');
%        export for Slides
    ratio = m(1)/m(2);
    w = 800;
    h = w/ratio;
    sfig(fig,ax,filename,format,[w h],'Helvetica',36,42);
elseif  strcmpi(type,'Helvet');
    sfig(fig,ax,filename,format,m,'Helvetica',24,30);
else
    error('Undefined Type');
end

set(fig,'PaperUnits','default');
set(fig,'PaperSize','default');
end

function resizefig(fig,m)
% Resize
w = m(1,1);
h = m(1,2);
set(fig,'PaperPositionMode','auto')
pos=get(fig,'Position');
pos(3)=w;
pos(4)=h;
set(fig,'Position',pos);
set(fig,'PaperSize',m);
end
function fsave(filename,format)
fname = sprintf('%s.%s',filename,format);
if strcmp(format,'eps')
    sub = '-depsc';
    print('-r300',sub,'-tiff',fname);
else
    sub = sprintf('-d%s',format);
    print('-r300',sub,fname);
%     savefig(gcf,filename);
end
end % end of resizefig

% save figure
function sfig(fig,ax,filename,format,m,type,ssize,lsize)
set(fig,'PaperUnits','points');
resizefig(fig,m);
% font
set( ax, 'FontName',type,'FontSize',ssize );
t = get(ax,'Title');
x = get(ax,'XLabel');
y = get(ax,'YLabel');
z = get(ax,'ZLabel');

lw = round(ssize/6);

% For subplot
if iscell(t)
    for i = 1:length(t)
        % Change large font size
        set(t{i},'FontSize',lsize);
        set(x{i},'FontSize',lsize);
        set(y{i},'FontSize',lsize);
        set(z{i},'FontSize',lsize);

        % Change small font size and linewidth
        n = size(ax(i).Children);
        for j = 1:n(1)
            if isa(ax(i).Children(j),'matlab.graphics.primitive.Text')
                set(ax(i).Children(j),'FontName',type,'FontSize',ssize);
            elseif isa(ax(i).Children(j),'matlab.graphics.chart.primitive.Line')
                set(ax(i).Children(j),'LineWidth',lw);
            elseif isa(ax(i).Children(j),'matlab.graphics.chart.primitive.ErrorBar')
                set(ax(i).Children(j),'LineWidth',lw);
            end
        end
    end
else
    set(t,'FontSize',lsize);
    set(x,'FontSize',lsize);
    set(y,'FontSize',lsize);
    set(z,'FontSize',lsize);
    n = size(ax.Children);
    for i = 1:n(1)
        if isa(ax.Children(i),'matlab.graphics.primitive.Text')
            set(ax.Children(i),'FontName',type,'FontSize',ssize);
        elseif isa(ax.Children(i),'matlab.graphics.chart.primitive.Line')
            set(ax.Children(i),'LineWidth',lw);
        elseif isa(ax.Children(i),'matlab.graphics.chart.primitive.ErrorBar')
            set(ax.Children(i),'LineWidth',lw);
        end
    end
end

% Save figure
fsave(filename,format);
end