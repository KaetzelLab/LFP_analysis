% ----------------------------------------------------------------------- %
% % Function polygonplot plots a kind of radial plot whose shape is a N-poly%
% gon, depending on the size of the data. It have also provide users the  %
% functionality of plotting a shaded error area if the dimension of the   %
% data is highly enough.                                                  %
%                                                                         %
%   Input parameters:                                                     %
%       - data:     Data matrix.                                          %
%           * Mean case: [N x nlines], with N corresponding to the number %
%             of points of each line (number of polygon sides) and nlines %
%             corresponding to the number of lines to plot.  
%           * Error bar case: [N x nlines x M], with M corresponding to   %
%             the number of observations.       
%             Connections/regions x groups x samples
%       - opt_axes: (Optional) Struct with the axes options.              %
%           * opt_axes.Ticks:       Vector that contains the axis ticks.  %
%           * opt_axes.Background:  Background color of the tick labels.  %
%           * opt_axes.Labels:      {1 x N} cell-vector with the axes lbls%
%       - opt_lines:(Optional) Struct with the lines options.             %
%           * opt_lines.LineWidth:  Width of the means.                   %
%           * opt_lines.LineStyle:  Line style of the means.              %
%           * opt_lines.Marker:     Marker of the means.                  %
%           * opt_lines.Color:      [nlines x 3] matrix with RGB colors.  %
%           * opt_lines.Labels:     Boolean. If true, text labels of each %
%                                   point are plotted.                    %
%           * opt_lines.Legend:     {1 x nlines} cell-matrix with the     %
%                                   legend of each line.                  %
%       - opt_area: (Optional) Struct with the shaded area options.       %
%           * opt_area.err:         Type of error to plot:                %
%                   if 'std',       one standard deviation;               %
%                   if 'sem',       standard error mean;                  %
%                   if 'var',       one variance;                         %
%                   if 'c95',       95% confidence interval.              %
%           * opt_area.FaceAlpha:   Alpha transparency constant.          %
%           * opt_area.Color:       [nlines x 3] matrix with RGB colors.  %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       d_ex = [2 3; 1 0; 0.1 3; -1 7; -0.2 0.9];                         %
%       data = cat(3,d_ex-0.5,d_ex,d_ex+0.7);                             %
%       polygonplot(data,~,~,9);                                                %
% ----------------------------------------------------------------------- %
%   Author:  Vctor Martnez Cagigal                                        %
%   Date:    21/03/2017                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %
function polygonplot(data,opt_axes,opt_lines,opt_area,minvalue,maxvalue)
    
    % Defaults
    if(nargin<4)
        opt_area.err = 'std';
        opt_area.FaceAlpha = 0.5;
    end
    if(nargin<3)
        opt_lines.LineWidth = 2;
        opt_lines.LineStyle = '-';
        opt_lines.Marker    = 'none';
    end
    if(nargin<2)
        opt_axes = [];
        opt_axes.Background = 'none';
    end
    if(nargin<1)
        error('Not enough parameters');
    end
    
    % Error detection
    if(length(size(data))>3)
        error('Data must be a MxNxnlines or Nxnlines matrix.');
    elseif(length(size(data))==3)       % Shaded
        [N,nlines,M] = size(data);
    elseif(length(size(data))==2)       % Only mean
        [N,nlines]   = size(data);
        M = 0;
    end  
    if(~isfield(opt_lines,'Color'))     % Color properties
        line_color = [];
    else
        if(size(opt_lines.Color,1)~=nlines)
            error('Number of colors must be equal to the number of lines to plot.');
        else
            if(~isfield(opt_area,'Color') && M~=0)
                error('Please, specify also the colors of the shaded areas.');
            else
                if(size(opt_area.Color,1)~=nlines)
                    error('Number of shaded areas must be equal to the number of lines to plot.');
                else
                    line_color = opt_lines.Color;
                    opt_lines = rmfield(opt_lines,'Color');
                end
            end
        end
    end
    if(~isfield(opt_lines,'Labels'))    % Text labels
        labels = false;
    else
        labels = opt_lines.Labels;
        opt_lines = rmfield(opt_lines,'Labels');
    end
    if(~isfield(opt_area,'Color'))
        opt_area.Color = 0.5.*ones(nlines,3);
    end
    if(~isfield(opt_area,'FaceAlpha'))
        opt_area.FaceAlpha = 0.5;
    end
    if(~isfield(opt_lines,'Legend'))
        leg = [];
    else
        leg = opt_lines.Legend;
        opt_lines = rmfield(opt_lines,'Legend');
    end
    if(isfield(opt_axes,'Labels'))
        if(length(opt_axes.Labels)~=N)
            error('You must provide N axis labels.');
        end
    end
    set(gcf,'Color',[1 1 1]);
    
    % Compute the isocurve Ticks
%     data_min = min(data(:));
%     data_max = max(data(:));
data_min = minvalue;
data_max = maxvalue;
    if(~isfield(opt_axes,'Ticks'))
        opt_axes.Ticks = linspace(data_min,data_max,6);
%         opt_axes.Ticks = logspace(data_min,data_max,6);
%         opt_axes.Ticks = [0 0.001 0.01 0.1 1];
    elseif(length(opt_axes.Ticks)<2)
        error('There should be more than 2 isocurves.');
    end
    r_iso = (opt_axes.Ticks(:)-data_min)/(data_max-data_min)*ones(1,N);
    th_iso = (2*pi/N)*(ones(length(opt_axes.Ticks),1)*(N:-1:1));
    [x,y] = pol2cart(th_iso, r_iso);
    h_iso = line([x,x(:,1)]',[y,y(:,1)]','LineWidth',0.5,'Color',0.85.*ones(1,3));
    for iso_id = 1:1:length(opt_axes.Ticks)   % Exclude axes from legend
        set(get(get(h_iso(iso_id),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    hold on;
    
    % Compute and plot the axes depending on N
    th_jump = (2*pi/N)*(ones(2,1)*(N:-1:1));
    radii   = [zeros(1,N); ones(1,N)];
    [x,y]   = pol2cart(th_jump, radii);
    h_axes  = line(x,y,'LineWidth',1,'Color','k');
    for ax_id = 1:1:N   % Exclude axes from legend
        set(get(get(h_axes(ax_id),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    hold on;
    
    % Display axis Ticks
    if(~isfield(opt_axes,'Background')), opt_axes.Background = 'none'; end
    loc = (opt_axes.Ticks(:)-data_min)/(data_max-data_min);
    odd_v = 0.03.*ones(size(opt_axes.Ticks));
    odd_v(2:2:end) = -0.7.*odd_v(2:2:end);
    for t = 1:1:length(opt_axes.Ticks)
        text(loc(t),odd_v(t),num2str(opt_axes.Ticks(t)),'Background',opt_axes.Background, 'fontweight', 'bold', 'FontSize',12);
    end
    
    % Display axis labels
    if(isfield(opt_axes,'Labels'))
        th_lbl = (2*pi/N)*(N:-1:1);
        r_lbl  = 1.1.*ones(1,N);
        [xlbl,ylbl] = pol2cart(th_lbl,r_lbl);
        for lid = 1:1:N
            text(xlbl(lid)-0.05,ylbl(lid),opt_axes.Labels{lid}, 'fontweight','bold','fontsize',14);
        end
    end
    
    % Plot the data
    opt_lines_str = adapt_options(opt_lines);
    if(M == 0)    % Only mean
        R  = ([data; data(1,:)]-data_min)/(data_max-data_min);
        TH = (2*pi/N)*((N:-1:0)'*ones(1,nlines));
        [X,Y] = pol2cart(TH, R);
        if(isempty(line_color)), plot(X,Y,opt_lines_str{:});
        else
            for n = 1:1:nlines,  plot(X(:,n),Y(:,n),opt_lines_str{:},'Color',line_color(n,:)); end
        end
        if(labels)
            for n = 1:1:nlines
                for t = 1:1:N
                    text(X(t,n),Y(t,n),num2str(data(t,n)),'FontSize',8);
                end
            end
        end
        axis square;
        axis off;
    else          % Shaded area
        % Computing the mean and standard deviation of the data matrix
        data_mean = nanmean(data,3);
        data_std  = nanstd(data,0,3);

        % Type of error plot
        switch(opt_area.err)
            case 'std', err = data_std;
            case 'sem', err = (data_std./sqrt(size(data,1)));
            case 'var', err = (data_std.^2);
            case 'c95', err = (data_std./sqrt(size(data,1))).*1.96;
        end
        
        % Plots
        m_down = data_mean-err./2; m_down = [m_down; m_down(1,:)];
        m_down = (m_down-data_min)/(data_max-data_min);
        m_up   = data_mean+err./2; m_up   = [m_up; m_up(1,:)];
        m_up   = (m_up-data_min)/(data_max-data_min);
        
        R  = ([data_mean; data_mean(1,:)]-data_min)/(data_max-data_min);
        TH = (2*pi/N)*((N:-1:0)'*ones(1,nlines));
        [X,Y] = pol2cart(TH, R);
        
        for n = 1:1:nlines
            % If error is not null, fill the area
            if nansum(nansum(err(:))) ~= 0
                [xa,ya] = pol2cart([TH(:,n); fliplr(TH(:,n))],[m_down(:,n); fliplr(m_up(:,n))]);
                pat = fill(xa,ya,opt_area.Color(n,:));
                hold on;
                set(pat, 'EdgeColor', 'none');
                set(pat, 'FaceAlpha', opt_area.FaceAlpha);
            else
                warning(['The ' opt_area.err 'in this dataset is null!']);
            end
        end
        if(isempty(line_color))
            h_leg = plot(X,Y,opt_lines_str{:}); hold on;
        else
            for n = 1:1:nlines
                h_leg(n) = plot(X(:,n),Y(:,n),opt_lines_str{:},'Color',line_color(n,:)); 
            end
        end
        if(labels)
            for n = 1:1:nlines
                for t = 1:1:N
                    text(X(t,n),Y(t,n),num2str(data_mean(t,n)),'FontSize',8);
                end
            end
        end
        axis equal;
        axis off;
        if(~isempty(leg))
            legend(h_leg([1:n]),leg, 'fontweight', 'bold', 'fontsize', 18);
        end
    end
    hold off;
end

% Function adapt_options adapts the input options struct in a dynamic
% cell-object that can directly pass Name-Value pair arguments through a function
function optList = adapt_options(optStruct)
    optList = {};
    for optField = fieldnames(optStruct)'
        optList{end+1} = char(optField);                % Name parameter
        optList{end+1} = optStruct.(char(optField));    % Value parameter
    end
end