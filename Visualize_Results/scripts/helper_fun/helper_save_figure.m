function helper_save_figure(cimec)
%% 'target right hemisphere'

conn_metric = cimec.conn_metric;
band_name = cimec.band_name;
seed = cimec.seed;
seed_target = cimec.seed_target;
analysis_type = cimec.statistics.analysis_type;
alpha = char(string(cimec.statistics.alpha));
correction = cimec.statistics.correction;
flag = cimec.flag;
duration = flag.duration;
dir_fig = flag.dir_fig;

if or(strcmp(seed_target, 'right-right'), strcmp(seed_target, 'left-right'))
    
    try
        if strcmp(flag.script, 'collapseTargets')
            seed_target = char(string(join(flag.seed_targets)));
        end
    catch
        fprintf('Almost ready...\n', i);
    end
    
    %% 0
    view(270,360)
    xlim([-100 100])
    ylim([-100 0])
    zlim([-20 150])
    set(gcf, 'Position', get(0, 'Screensize'));
    
    %% 1
    % save the figure
    lighting none
    material dull
    lighting phong
    cam1 = camlight('left');
    
    % for light
    view(360,270)
    cam2 = camlight('headlight');
    view(360,360)
    cam3 = camlight('left');
    
    if strcmp(analysis_type, 'exploratory')
        cam4 = camlight('right');
        
        %% 1
        number = '1';
        name = [number '_' conn_metric '_' band_name '_' seed '_' seed_target '_' analysis_type '_' alpha '_' correction '_' duration];
        print(gcf,[dir_fig name '.png'],'-dpng','-r250');
        
        %% 3
        view(180,360)
        cam3 =camlight('headlight');
        number = '3';
        name = [number '_' conn_metric '_' band_name '_' seed '_' seed_target '_' analysis_type '_' alpha '_' correction '_' duration];
        print(gcf,[dir_fig name '.png'],'-dpng','-r250');
    else
        %% 1
        number = '1';
        name = [number '_' conn_metric '_' band_name '_' seed '_' seed_target '_' analysis_type '_' alpha '_' correction '_' duration];
        print(gcf,[dir_fig name '.png'],'-dpng','-r250');

        %% 2
        view(315,5)
        number = '2';
        name = [number '_' conn_metric '_' band_name '_' seed '_' seed_target '_' analysis_type '_' alpha '_' correction '_' duration];
        print(gcf,[dir_fig name '.png'],'-dpng','-r250');

        %% 3
        view(180,360)
        cam3 =camlight('headlight');
        number = '3';
        name = [number '_' conn_metric '_' band_name '_' seed '_' seed_target '_' analysis_type '_' alpha '_' correction '_' duration];
        print(gcf,[dir_fig name '.png'],'-dpng','-r250');

        %% 4
        if strcmp(conn_metric, 'pdc')
            view(360,270)
            number = '4';
            name = [number '_' conn_metric '_' band_name '_' seed '_' seed_target '_' analysis_type '_' alpha '_' correction '_' duration];
            print(gcf,[dir_fig name '.png'],'-dpng','-r250');
        end
    end

elseif or(strcmp(seed_target, 'right-left'), strcmp(seed_target, 'left-left'))
    
    try
        if strcmp(flag.script, 'collapseTargets')
            seed_target = char(string(join(flag.seed_targets)));
        end
    catch
        fprintf('Almost ready...\n', i);
    end

    
    
    %% RUN
    %% 0
    view(270,360)
    xlim([-100 100])
    ylim([0 100])
    zlim([-20 150])
    set(gcf, 'Position', get(0, 'Screensize'));
    
    %% 1
    % save the figure
    lighting none
    material dull
    lighting phong
    cam1 = camlight('right');
    
    % for light
    view(360,270)
    cam2 = camlight('headlight');
    view(180,360)
    cam3 = camlight('right');

    if strcmp(analysis_type, 'exploratory')
        %% 1
        cam4 = camlight('left');
        
        number = '1';
        name = [number '_' conn_metric '_' band_name '_' seed '_' seed_target '_' analysis_type '_' alpha '_' correction '_' duration];
        print(gcf,[dir_fig name '.png'],'-dpng','-r250');
    
        %% 3
        view(360,360)
        cam3 =camlight('headlight');
        number = '3';
        name = [number '_' conn_metric '_' band_name '_' seed '_' seed_target '_' analysis_type '_' alpha '_' correction '_' duration];
        print(gcf,[dir_fig name '.png'],'-dpng','-r250');

    else
        %% 1
        number = '1';
        name = [number '_' conn_metric '_' band_name '_' seed '_' seed_target '_' analysis_type '_' alpha '_' correction '_' duration];
        print(gcf,[dir_fig name '.png'],'-dpng','-r250');

        %% 2
        view(225,5)
        number = '2';
        name = [number '_' conn_metric '_' band_name '_' seed '_' seed_target '_' analysis_type '_' alpha '_' correction '_' duration];
        print(gcf,[dir_fig name '.png'],'-dpng','-r250');

        %% 3
        view(360,360)
        cam3 =camlight('headlight');
        number = '3';
        name = [number '_' conn_metric '_' band_name '_' seed '_' seed_target '_' analysis_type '_' alpha '_' correction '_' duration];
        print(gcf,[dir_fig name '.png'],'-dpng','-r250');

        %% 4
        if strcmp(conn_metric, 'pdc')
            view(360,270)
            number = '4';
            name = [number '_' conn_metric '_' band_name '_' seed '_' seed_target '_' analysis_type '_' alpha '_' correction '_' duration];
            print(gcf,[dir_fig name '.png'],'-dpng','-r250');
        end
    end
end

end
%%