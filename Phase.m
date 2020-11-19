function [FileList,Result] = Phase(InputDir)

    clc;
    close all;
    clear all;

    %InputDir='D:\Temp\Maya\Phase\Images\22.2.19';
    ShowFigure = false;
    SaveImages = false;

    FileList=GetFileNames(InputDir,'tif');

    % Initialize results and template
    Results(1, size(FileList,1)) = struct('Directory',string([]),...
                                          'FileName',string([]),...
                                          'Well',string([]),...
                                          'WellCoordinate',string([]),...
                                          'Date',string([]),...
                                          'Time',string([]),...
                                          'NumberOfCells',[]);

    % Loop and calculate
    WaitBarHandle = waitbar(0,'Please wait...');
    for n=1:size(FileList,1)

        % Parse name
        FullFileName = string(FileList(n,:));
        PartialFileName = FullFileName.split('\');
        Results(n).Directory = join(PartialFileName(1:(end-1)),'/');
        PartialFileName = PartialFileName(end);
        Results(n).FileName = PartialFileName;    

        % Find parameters
        Parameters = PartialFileName.split('.');
        Parameters = Parameters(1);
        Parameters = Parameters.split('_');
        Results(n).Well = Parameters(2);
        Results(n).WellCoordinate = Parameters(3);
        Results(n).Date = Parameters(4);
        Results(n).Time = Parameters(5);
        clear Parameters PartialFileName;  
        
        % Set output file name
        OutputFileName = FullFileName.split('.');
        OutputFileName = join(OutputFileName(1:(end-1)),'.') + string('.pgm');
        
        if(~exist(OutputFileName.char,'file'))

            % Load image
            Im=imread(FullFileName.char());

            % Convert to gray
            Im = rgb2gray(Im);        

            % Find cells
            ImCells = zeros(size(Im), 'uint8');
            ImCells(Im < 20) = 255;
            ImCells(Im > 235) = 255;
            SE = strel('disk',3);
            ImCells = imdilate(ImCells,SE);

            % Show image with cells
            if(ShowFigure)
                figure;
                H(1) = subplot(1,2,1);
                image(Im);
                colormap(gray(256));
                axis image; 
                H(2) = subplot(1,2,2);
                image(ImCells);
                colormap(gray(256));
                axis image;
                linkaxes(H,'xy');
            end;

            % Find the lines
            ImResult = FindLines(Im, ImCells, ShowFigure);

            % Show image with cells
            if(ShowFigure)
                figure;
                H(1) = subplot(1,2,1);
                image(Im);
                colormap(gray(256));
                axis image;
                H(2) = subplot(1,2,2);
                imagesc(ImResult);
                colormap(gray(256));
                axis image;
                linkaxes(H,'xy');
            end;

            % Save output image
            if(SaveImages)
                OutputFileName = FullFileName.split('.');
                OutputFileName = join(OutputFileName(1:(end-1)),'.') + string('.pgm');
                imwrite(ImResult,OutputFileName.char(),'pgm');
            end;
        else
            ImResult=imread(OutputFileName.char());
        end;            

        % Update results
        Results(n).NumberOfCells=sum(ImResult(:)>0);
        clear ImResult;

        % Update progress bar
        waitbar(n/size(FileList,1),WaitBarHandle);        
    end;
    close(WaitBarHandle);
    
    % Save results in matlab format
    save([InputDir '\MayaResults.mat'],'Results');
    
    % Find wells
    Wells = Results(1).Well;
    for m=2:numel(Results) 
        if(~Wells.contains(Results(m).Well))
            Wells(end+1) = Results(m).Well;
        end;
    end;
    
    % Loop over all wells
    for m=1:numel(Wells)
        
        % Find well coordinates number of times
        Indices = find([Results(:).Well] == Wells(m));
        WellCoordinates = Results(Indices(1)).WellCoordinate;
        NumberOfTimes = sum([Results(Indices).WellCoordinate] == WellCoordinates(1));
        for n=1:numel(Indices)
            if(~WellCoordinates.contains(Results(Indices(n)).WellCoordinate))
                WellCoordinates(end+1) = Results(Indices(n)).WellCoordinate;
                NumberOfTimes = max(NumberOfTimes, sum([Results(Indices).WellCoordinate] == WellCoordinates(end)));
            end;            
        end;
       
        % Initialize output
        OutputData = zeros(NumberOfTimes,numel(WellCoordinates),'single');
        
        % Loop over all well coordinates
        DateTime = strings(NumberOfTimes,1);
        for n=1:numel(WellCoordinates)
            
            % Find all indices matching the well and coordinates
            Indices = find(([Results(:).Well] == Wells(m)) & ([Results(:).WellCoordinate] == WellCoordinates(n)));
            
            % Loop and fill data
            for k=1:numel(Indices)
                OutputData(k,n) = Results(Indices(k)).NumberOfCells;
                DateTime(k) = Results(Indices(k)).Date + Results(Indices(k)).Time;
            end;
        end;
        
        figure;
        plot(1:size(OutputData,1),mean(OutputData,2),1:size(OutputData,1),median(OutputData,2));
        title(Wells(m).char);
        
        % Write results
        fid = fopen([InputDir '/MayaResults_' Wells(m).char '.csv'],'wt');
        fprintf(fid,'DateTime');
        for n=1:numel(WellCoordinates)
            fprintf(fid,',%s',WellCoordinates(n).char);
        end;
        for n=1:numel(WellCoordinates)
            fprintf(fid,',%s',WellCoordinates(n).char);
        end;
        fprintf(fid,'\n');
        for n=1:numel(DateTime)
            fprintf(fid,'%s',DateTime(n).char);
            for k=1:size(OutputData,2)
                fprintf(fid,',%f',OutputData(n,k));
            end;
            for k=1:size(OutputData,2)
                fprintf(fid,',%f',OutputData(n,k)/(OutputData(1,k)+eps));
            end;
            fprintf(fid,'\n');
        end;
        fclose(fid);            
    end;
end

function OutputImage = FindLines(InputImage, CellImage, ShowFigure)

    % Initialize show figure
    if(nargin < 2)
        ShowFigure = false;
    end;
    
    % Estimate the MED and median
    MED = single(median(InputImage(:)));
    MAD = median(abs(single(InputImage(:))-MED));
    
    % Remove everything that is not backgroud
    OutputImage = zeros(size(InputImage),'uint8');
    OutputImage(InputImage < (MED-1*1.4826*MAD)) = 255;
    
    % Apply morphological operations to clean
    SE = strel('disk',1,0);
    OutputImage = imclose(OutputImage,SE);
    
    % Remove cells from results
    OutputImage = bitand(OutputImage,uint8(255)-CellImage);    
    
    % Perform connected components
    CC = bwconncomp(OutputImage);
    StatsCC = regionprops(CC,'BoundingBox','Area','Extent','Eccentricity','MajorAxisLength','MinorAxisLength','Orientation');
    
    % Remove box shaped 
    MinSizeX = 5;
    MinSizeY = 5;
    BoundingBoxes = reshape([StatsCC(:).BoundingBox],4,numel(StatsCC));
    Indices = find((BoundingBoxes(3,:) <= MinSizeX) & (BoundingBoxes(4,:) <= MinSizeY));
    for m=1:numel(Indices)
        OutputImage(CC.PixelIdxList{Indices(m)}) = 0;
    end;
    
    % Remove by length
    MinLength = 15;
    EllipseAxis = max([[StatsCC(:).MinorAxisLength];[StatsCC(:).MajorAxisLength]]);
    Indices = find(EllipseAxis < MinLength);
    for m=1:numel(Indices)
        OutputImage(CC.PixelIdxList{Indices(m)}) = 0;
    end;
    
    
    
    
    
    
    
    
%     TempImage = conv2(single(OutputImage),ones(5,5),'same');
%     OutputImage(TempImage < (9*255+1)) = 0;
%     
%     TempImage = conv2(single(OutputImage),ones(3,3),'same');
%     OutputImage(TempImage < (2*255+1)) = 0;
       
    % Display images
    if(ShowFigure)
        figure;
        H(1)=subplot(1,2,1);
        image(InputImage);
        colormap(gray(256));
        axis image;
        H(2)=subplot(1,2,2);    
        image(OutputImage);
        colormap(gray(256));
        axis image;
        linkaxes(H,'xy');
    end;  
end

% ************** Count cells function
function [NumberOfCells,CellsPerArea,Polygon]=CountCellsMaya(ImageName)

    % Read image
    ImOrig=imread(ImageName,'tif');

    % Show image
    H=figure('units','normalized','outerposition',[0 0 1 1]);
    image(ImOrig);
    hold on;
    axis equal;

    % Get Polygon of interest
    X=[];
    Y=[];
    while(1)
        [X1,Y1,Button]=ginput(1);
        if(Button == 3)
            X=[X X(1)];
            Y=[Y Y(1)];
            plot(X((end-1):end),Y((end-1):end),'-gx');
            break;
        end;
        X=[X X1];
        Y=[Y Y1];
        if(size(X,2) >= 2)
            plot(X((end-1):end),Y((end-1):end),'-gx');
        else
            plot(X,Y,'gx');        
        end;
    end;
    Polygon=zeros(2,size(X,2));
    Polygon(1,:)=X;
    Polygon(2,:)=Y;
    clear X Y X1 Y1;

    % Close figure
    pause(0.1);
    close(H);
    clear H;

    % Transform to gray
    Im=double(ImOrig(:,:,2));

    % Create template
    Std=5;
    Template=zeros(3*Std*2+1,3*Std*2+1);
    for m=(-3*Std):1:3*Std
        for n=(-3*Std):1:3*Std
            Template(m+3*Std+1,n+3*Std+1)=exp(-(m^2+n^2)/2/(Std^2));
        end;
    end;

    % % Show template
    % figure;
    % imagesc(Template);

    % Perform normalized cross correlation
    ImCorr=normxcorr2(Template,Im);
    ImCorr=ImCorr(((size(Template,1)+1)/2):((size(Template,1)+1)/2+size(Im,1)-1),((size(Template,2)+1)/2):((size(Template,2)+1)/2+size(Im,2)-1));

    % % Show results of correlation
    % figure;
    % H(1)=subplot(2,1,1);
    % imagesc(Im);
    % colormap(gray(256));
    % axis equal;
    % H(2)=subplot(2,1,2);
    % imagesc(ImCorr);
    % colormap(gray(256));
    % axis equal;
    % linkaxes(H,'xy');

    % Count number of cells
    NumberOfCells=0;
    Cells=[];
    for m=1:(2*Std):size(ImCorr,1)
        for n=1:(2*Std):size(ImCorr,2) 

            Index1=min(m+2*Std+2,size(ImCorr,1));
            Index2=min(n+2*Std+2,size(ImCorr,2));      

            ImBin=ImCorr(m:Index1,n:Index2);
            MaxBin=max(ImBin(:));
            [I,J]=find(ImBin==MaxBin);
            if((MaxBin >= 0.4) && (I ~= 1) && (J ~= 1) && (I ~= size(ImBin,1)) && (J ~= size(ImBin,2)))                      

                % Check if inside polygon
                IN=inpolygon(n+J-1,m+I-1,Polygon(1,:),Polygon(2,:));
                if(IN)
                    NumberOfCells=NumberOfCells+1;
                    Cells(NumberOfCells,:)=[m+I-1,n+J-1];
                end;
            end;

        end;
    end;

    % Calculate area of polygon
    Area=polyarea(Polygon(1,:),Polygon(2,:));

    % Calculate cells per area
    CellsPerArea=NumberOfCells/Area;

    % Display results
    % hold on;
    % plot(Cells(:,2),Cells(:,1),'bx');
    H=figure('units','normalized','outerposition',[0 0 1 1]);
    image(ImOrig);
    axis equal;
    hold on;
    plot(Cells(:,2),Cells(:,1),'gx');
    plot(Polygon(1,:),Polygon(2,:),'b');
    title(['Number of cells: ' num2str(NumberOfCells) ', Density of cells: ' num2str(CellsPerArea) ' [Cells/Pixels^2]']);
    KeyDown=waitforbuttonpress;
    close(H);

    % Save results
    ResultsFileName=[ImageName(1:strfind(ImageName,'.tif')) 'mat'];
    save(ResultsFileName,'ImageName','Polygon','Cells','NumberOfCells','CellsPerArea');
end


% ************** Get file names function
function FoundFileNames=GetFileNames(SearchDir,FileType)

    FoundFileNames=char([]);
    if(isempty(SearchDir)||isempty(FileType))
        return;
    end;

    FileNames=dir(SearchDir);
    for n=1:size(FileNames,1)    
        FullFileName=[SearchDir '\' FileNames(n).name];    
        if(FileNames(n).isdir)
            if(~strcmp(FileNames(n).name,'.') && ~strcmp(FileNames(n).name,'..'))
                NextLevelFileNames=GetFileNames(FullFileName,FileType);
                if(~isempty(NextLevelFileNames))
                    FoundFileNames((size(FoundFileNames,1)+1):(size(FoundFileNames,1)+size(NextLevelFileNames,1)),1:size(NextLevelFileNames,2))=NextLevelFileNames;            
                end;
            end;
        elseif(~isempty(strfind(FileNames(n).name,['.' FileType])))        
            FoundFileNames(size(FoundFileNames,1)+1,1:size(FullFileName,2))=FullFileName;        
        end;
    end;
end







