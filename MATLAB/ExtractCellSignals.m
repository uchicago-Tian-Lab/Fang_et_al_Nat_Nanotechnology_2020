close all; clear all;clearvars; close all; clc;

folder = 'J:\WeiLabHD\recordings\NewRig\20190305';

moviefiles = dir([folder '\*.zip']);
cell_area = 0;

pixel_xcorr = input('Pixel Cross-Correlation Analysis (For XCorr Bkg must be substracted)?\n');

smoothing = input('Smoothing?\n');

unmixing = input('Unmixing?\n');

substract_bgd = input('Substract Background?\n');

%retinaInfo = dlmread([folder '\RetinaInfo.txt']);


for i = 1:size(moviefiles,1)
    close all
    clearvars -except paras folder moviefiles i retinaInfo cell_area num_stimuli compare pixel_xcorr smoothing unmixing substract_bgd
    movieNum = str2num(moviefiles(i,1).name(7:10));
    moviestring = moviefiles(i,1).name(7:10);
    pname = dir([folder '\TSeries*' moviestring]);
    tseriesxml = [pname.name '.xml'];
    pname = [folder '\' pname.name];

%     movieInf = find(retinaInfo(:,1) == movieNum);
%     ch_num = retinaInfo(movieInf, 7);

    %% IMPORT TIFF FILES%%%%


    [~, tseries] = fileparts(pname);
    xmlname = dir([pname '\TSeries*' '.xml']);
    xmlfile = fileread([pname '\' xmlname(1).name]);
    chImgNames =  regexp(xmlfile, ...
        'File channel="(\d*)".*?filename="(.*?)"', 'tokens');
    %framerate = regexp(xmlfile, 'key="framePeriod" permissions="ReadWriteSave" value="(\d*)"' ;
    ch = cellfun(@(x) str2double(x{1}),chImgNames);
    filtChImageNames = chImgNames(ch == 2);
    imgNames = cellfun(@(x) x{2},filtChImageNames, 'UniformOutput', 0);



    %%TURN TIFF FILES INTO MATRICES%%%%
    % Read individual tif files  img
    for n = 1:numel(imgNames)
        I=imread([pname '\' imgNames{n}]);
        movieM{n} = I;
    end

    clear chImgNames ch filtChImageNames imgNames

    if unmixing == 1
        chImgNames =  regexp(xmlfile, ...
        'File channel="(\d*)".*?filename="(.*?)"', 'tokens');
        %framerate = regexp(xmlfile, 'key="framePeriod" permissions="ReadWriteSave" value="(\d*)"' ;
        ch = cellfun(@(x) str2double(x{1}),chImgNames);
        filtChImageNames = chImgNames(ch == 1);
        imgNames = cellfun(@(x) x{2},filtChImageNames, 'UniformOutput', 0);

        for n = 1:numel(imgNames)
            I_Ch1=imread([pname '\' imgNames{n}]);
            movieM_Ch1{n} = I_Ch1*0.08;

            movieM{n} = movieM{n} - movieM_Ch1{n};

        end
    end



    clear xmlfile;
    [rows, cols] = size(movieM{1,1});


    %%%% IMPORT IMAGEJ ROIS%%%%
    [rois] = ReadImageJROI([folder '\RoiSet' num2str(movieNum) '.zip']);

    cellSignals = zeros(length(movieM), length(rois));
    coords = zeros(length(rois),4);

    for i = 1:length(rois)
        coords(i,:) = rois{1,i}.vnRectBounds;
        if coords(i, 3) > rows
            coords(i,3) = rows;
        elseif coords(i,4) > cols
            coords(i,4) = cols;
        end
    end

    coords(coords < 1) = 1;
    coords(find(coords(:,3) > rows),3) = rows;
    coords(find(coords(:,4) > cols),4) = cols;

    %%%%%EXTRACT MEAN VALUE FOR EACH ROI IN EVERY FRAME%%%%%
    for i=1:length(movieM)
        for k= 1:length(rois)
            regions{i,k} = movieM{1,i}(coords(k,1):coords(k,3), coords(k,2): coords(k,4));
            cellSignals(i,k) = mean(mean(regions{i,k}));
        end
    end

    %%PIXEL CROSS CORRELATION ANALYSIS FOR EACH ROI%%%%%

    if substract_bgd == 1

        for n = 1:(size(cellSignals,2) - 1)
            cellSignals_bgS(:,n) = cellSignals(:,n) - cellSignals(:,end);
        end

        if pixel_xcorr == 0
            cellSignals2 = cellSignals;
            cellSignals = [cellSignals_bgS, cellSignals2(:,end)];
        end
    end


    if pixel_xcorr == 1
        for q = 1:(length(rois)-1)
            num_pixels = numel(regions{1,q});
            [nrows ncol] = size(regions{1,q});
            for z = 1:length(regions)
                for y = 1:num_pixels
                    row = floor(y/ncol) + 1;
                    if row == (nrows + 1)
                        row = nrows;
                    end
                    col = mod(y, ncol);
                    if col == 0
                        col = ncol;
                    end
                    region_deconv(z,y) = double(regions{z,q}(row,col));
                    pixel_coords{1,q}(y,:) = [row, col];
                end
            end
            for d = 1:size(region_deconv,2)
                region_deconv(:,d) = region_deconv(:,d) - cellSignals(:,end);
            end

            std_pixels = std(region_deconv);
            region_deconv = region_deconv(:,find(std_pixels > (mean(std_pixels) - std(std_pixels))));
            pixel_coords{1,q} = pixel_coords{1,q}(find(std_pixels > (mean(std_pixels) - std(std_pixels))),:);

            rho_corr = zeros(size(region_deconv,2));
            distance_matrix = zeros(size(region_deconv,2));

            for b = 1:size(region_deconv,2)
                for a = 1:size(region_deconv,2)
                    [crosscorr , l] = xcorr(region_deconv(:,b), region_deconv(:,a), 'coeff');
                    rho_corr(b,a) = crosscorr(find(l == 0));
                    if b == a
                        rho_corr(b,a) = 0;
                    end
                    distance_matrix(b,a) = pdist([pixel_coords{1,q}(b,:);pixel_coords{1,q}(a,:)],'euclidean');

                end
            end

            elements = intersect(find(distance_matrix > 3 & distance_matrix < 11 ), find(rho_corr > (mean(rho_corr(:)) + std(rho_corr(:)))));
            [m, j] = ind2sub([size(region_deconv,2),size(region_deconv,2)], elements);
            elements = unique(m);
            region_deconv = region_deconv(:,elements);
            cellSignals2(:,q) = cellSignals(:,q);
            cellSignals(:,q) = mean(region_deconv,2);
            if unique(isnan(cellSignals(:,q))) == 1
                cellSignals(:,q) = cellSignals_bgS(:,q);
            end

            clear region_deconv c l
        end
    end

    if smoothing == 1
        %cellSignals = DataSmooth(cellSignals, 3, 'average');


        for q = 1:(length(rois)-1)
            plot(cellSignals(:,q))
            title(['Movie ' num2str(movieNum) ' Cell ' num2str(q)])
            ylim([0 1500])
            fig = gcf;
            saveas(fig,strcat(folder,'\CellSignals\','Movie',num2str(movieNum),'Cell',num2str(q)),'fig')
            saveas(fig,strcat(folder,'\CellSignals\','Movie',num2str(movieNum),'Cell',num2str(q)),'eps')
            saveas(fig,strcat(folder,'\CellSignals\','Movie',num2str(movieNum),'Cell',num2str(q)),'png')
            close all
        end
    end



    output_fn = strcat('CellSignals', num2str(movieNum),'.txt');
    output_path = [folder '\CellSignals\' output_fn];
    dlmwrite(output_path, cellSignals, 'delimiter', '\t','newline','pc')
end
