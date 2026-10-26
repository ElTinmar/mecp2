clear all; close all; clc;

DevelopperMode=0; %Advanced parameters.
ReshapeDataDevelopperMode=0; %Advanced parameters.
MaskDevelopperMode=0; %Advanced parameters.
SegmentationMode=1; %Advanced parameters.
NewAxisLineDevelopperMode=0;
Normalization_Parameter=1;

% % % % % % %fis={'Fish 1-MeCP2-6dpf_060318' 'Fish 2-MeCP2-6dpf_120318' 'Fish 3-MeCP2-6dpf_120318' 'Fish 4-MeCP2-6dpf_120318' 'Fish 5-MeCP2-6dpf_120318' 'Fish 6-MeCP2-6dpf_120318' 'Fish 7-MeCP2-6dpf_120318' 'Fish 9-MeCP2-6dpf_120318'}; %
fis={'Fish 1-wt-6dpf_050318' 'Fish 2-wt-6dpf_050318' 'Fish 3-wt-6dpf_050318' 'Fish 4-wt-6dpf_060318' 'Fish 5-wt-6dpf_060318' 'Fish 6-wt-6dpf_060318' 'Fish 7-wt-6dpf_060318' 'Fish 9-wt-6dpf_060318' 'Fish 1-MeCP2-6dpf_060318' 'Fish 2-MeCP2-6dpf_120318' 'Fish 3-MeCP2-6dpf_120318' 'Fish 4-MeCP2-6dpf_120318' 'Fish 5-MeCP2-6dpf_120318' 'Fish 6-MeCP2-6dpf_120318' 'Fish 7-MeCP2-6dpf_120318' 'Fish 9-MeCP2-6dpf_120318'}; %'Fish 1-wt-6dpf_050318' 'Fish 2-wt-6dpf_050318' 'Fish 3-wt-6dpf_050318' 'Fish 4-wt-6dpf_060318' 'Fish 5-wt-6dpf_060318' 'Fish 6-wt-6dpf_060318' 'Fish 7-wt-6dpf_060318' 'Fish 8-wt-6dpf_060318' 'Fish 9-wt-6dpf_060318' 'Fish 8-wt-6dpf_060318'
side={'Left', 'Rigth'};
Tectum_Begin={24 15 15 1 1 1 1 1 1 1 1 1 21 1 3 1}; %{80 80};%
Tectum_Stop={110 101 101 87 87 87 87 87 87 87 87 87 106 87 89 87}; %{84 84};% 110
Tectum_Start={77 56 61 27 43 27 33 20 29 37 31 27 58 38 23 28}; %{84 84};% 77
XPoint.Left={900 400 500 400 700 700 300 800 350 450 250 700 800 550 450 700};
YPoint.Left={500 900 1000 50 1000 900 150 800 100 0 750 50 100 1000 0 100};
XPoint.Rigth={250 550 400 800 150 100 850 150 800 800 600 750 700 200 800 700};
YPoint.Rigth={900 100 100 600 500 400 550 450 650 800 0 950 1000 400 800 1000};

for overlap=[0.45]     %0.45 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
for ANI=1:length(fis)%[6]%[2]%[16]%
for si=1:length(side)
clearvars -except ANI fis overlap side Tectum_Start Tectum_Stop Tectum_Begin DevelopperMode ReshapeDataDevelopperMode MaskDevelopperMode SegmentationMode NewAxisLineDevelopperMode si XPoint YPoint Normalization_Parameter
%%%%%%%%%%%%%%%%%%%%%%
    %%    KEY SCRIPT    %%
%%%%%%%%%%%%%%%%%%%%%%

%Version 07102016. Analysis cell type density in WT and MeCP2 mutant fish.
%Parts i. Preamble: Constants and display parameters ii. Load data, single
%or multiple, WT and MeCP2 folders. iii. Analyses features: gaussian
%filtering and masks labels iv. Separate/Process layers: files display,
%color separation and DAPI contrast adaptation v. Get colors for each patch
%vi. Compute cell types proportions/display/save data

%% ----PREAMBLE -----%

set(0,'DefaultFigureWindowStyle','docked');
set(0,'DefaultAxesFontName', 'FreeSans'); set(0,'DefaultAxesFontSize', 14); %Axes
set(0,'DefaultTextFontname', 'FreeSans'); set(0,'DefaultTextFontSize', 14); %Text
load('/media/emiliano/CODES/Segmentation/cmapBlackGreen.mat');
load('/media/emiliano/CODES/Segmentation/cmapBlackRed.mat');
load('/media/emiliano/CODES/Segmentation/cmapBlackCyan.mat');

% -- Constants
P = struct('vGlut2',1100,'Gad1b',1100,'DAPI',400); %IMPORTANT(biological), all the region with areas higher than this value will be removed, this is a necessary step when using 2d correlation
% P.vGlut2=1100;  %IMPORTANT(biological), all the region with areas higher than this value will be removed, this is a necessary step when using 2d correlation
% P.Gad1b=1100;  %IMPORTANT(biological), all the region with areas higher than this value will be removed, this is a necessary step when using 2d correlation
% P.DAPI=400;  %IMPORTANT(biological), all the region with areas higher than this value will be removed, this is a necessary step when using 2d correlation
% open=20;   %IMPORTANT(biological), all region smaller than this value will be removed
open=struct('vGlut2',15,'Gad1b',15,'DAPI',15);  %IMPORTANT(biological), all region smaller than this value will be removed
% open.vGlut2=80;   %IMPORTANT(biological), all region smaller than this value will be removed
% open.Gad1b=80;   %IMPORTANT(biological), all region smaller than this value will be removed
% open.DAPI=15;   %IMPORTANT(biological), all region smaller than this value will be removed
taille=7;   %IMPORTANT(biological), this is the more or less the diameter of the nuclei, it's necesseray to create a gaussian filter and to do a 2D correlation
FG=5;   %IMPORTANT, this parameter will be use to filter the 2d correlation map
% cut=1; %IMPORTANT, used to determine the threshold to create the mask
cut=struct('vGlut2',1.1,'Gad1b',0.7,'DAPI',1.1);  %IMPORTANT, used to determine the threshold to create the mask
% cut.vGlut2=1.2; %IMPORTANT, used to determine the threshold to create the mask
% cut.Gad1b=1.2; %IMPORTANT, used to determine the threshold to create the mask
% cut.DAPI=1.6; %IMPORTANT, used to determine the threshold to create the mask
% -- Options

Contrast_average = [133.8120  131.6593  135.1084];
Normalize_contrast_average = [133.8120  131.6593  135.1084];
Image_average = [29.9337   38.1002   71.5755];
Normalize_image_average =  [29.9337   38.1002   71.5755];

%% ----ANALYSES FEATS -----%        
    % -- Gaussian Filter
    mu = taille/2-0.5;sigma = taille/2-0.5;gauss=zeros(taille);
    for X=1:taille
        for Y=1:taille
            gauss(X,Y) = exp((-(X-mu).^2-(Y-mu).^2)./sigma);
        end
    end
    h = fspecial('gaussian', FG, FG);
%% ----LOAD DATA -----%
% -- Load several files
%%if DevelopperMode==1
    %%WTFolder='/media/jboulang/DATA2/DATA2/MECP2_PROJECT/WT/';
    %%MeCP2Folder='/media/jboulang/DATA2/DATA2/MECP2_PROJECT/MeCP2/';
    %%EnucFolder1='/home/timothee/Bureau/PoissonEnuclee1.tif/';
    %%EnucFolder2='/home/timothee/Bureau/PoissonEnuclee2.tif/';
    %%EnucFolder3='/home/timothee/Bureau/PoissonEnuclee3.tif/';
    %%[WT,MeCP2]=Load_Immunos(WTFolder,MeCP2Folder);
    %%[PoissonEnuclee1.tif,PoissonEnuclee2.tif,PoissonEnuclee3.tif]=Load_Immunos(EnucFolder1,EnucFolder2,EnucFolder3);
%%end

% -- Load single file
%Rq: Clear all files within the folder except the loaded file.
%%picPath='E:\Recherche_Virginie\MeCP2_Confocal\Image_Confocal\Immuno_WT\Experiment_050916_6j_WT';
%%resultsFolder='E:\Recherche_Virginie\MeCP2_Confocal\Segmentation';
FolderIn='/media/emiliano/DATA/Segmentation/ALL/';
FileIn=string(fis(ANI)); %'poiss-enucl-110817.lif - TileScan_002_Merging001'; % 'Fish 9-MeCP2-6dpf_120318';  string(fis(ANI))
picPath=sprintf('%s%s%s',FolderIn,FileIn,'.tif');
cutName=strfind(picPath,'.tif');
tiff_info = imfinfo(picPath); 
resultsFolder='/media/emiliano/CODES/Segmentation/Resultats'; %
% -- Labelings/ROIs
maskNames={'Tectum'};%,'Habenula'}; %For each brain part (eg. Tectum) a mask will be created, just add annother argument (eg. 'Habenula') to make another mask per stack fraome
labelingNames={'vGlut2','Gad1b','DAPI'}; %Names of the labeling, the order corresponses to the loaded colormaps, DAPI is necessary for neuronal segmentation.
colorsCount=length(labelingNames);
TotalFrames=numel(tiff_info)/colorsCount;
Masks=struct();
Masks.Tectum.Start=Tectum_Begin{ANI};%1; 
Masks.Tectum.Stop=Tectum_Stop{ANI}; %TotalFrames; %Start and Stop frames of the masks (eg. Masks.Habenula.Start=30; Masks.Habenula.Stop=31;)
Masks.Habenula.Start=1; Masks.Habenula.Stop=3; %Start and Stop frames of the masks (eg. Masks.Habenula.Start=30; Masks.Habenula.Stop=31;)
MaSt=Tectum_Start{ANI};%31; %stak to start to make mask
%% ----SEPARATE LAYERS-----%
if ReshapeDataDevelopperMode==1
    % -- Display frames one by one
    %Check whether each color corresponds to the expected labeling (eg. DAPI is blue).
   
     for j = 1:colorsCount:size(tiff_info, 1)
         figure(10);
         temp_tiff= imread(picPath, j);
         sub1=subplot(1,colorsCount,1); imagesc(imadjust(temp_tiff)); axis image; title(labelingNames{1})
         colormap(sub1,cmapBlackRed); freezeColors
         
         temp_tiff= imread(picPath, j+1);
         sub2=subplot(1,colorsCount,2);
         imagesc(imadjust(temp_tiff)); axis image; title(labelingNames{2})
         colormap(sub2,cmapBlackGreen); freezeColors
         if colorsCount>2
             temp_tiff= imread(picPath, j+2);
             subplot(1,colorsCount,3);
             imagesc(imadjust(temp_tiff)); axis image; title(labelingNames{3})
             colormap(cmapBlackCyan); freezeColors
         end
     end
     % -- Reshape Data
     %Creates a single file per color where color patches will be inspected.
     %Make sure that these files have not been created before other wise you
     %will end up with twice as much frames.
     
     countDAPI=0; countvGlut2=0; countGad1b=0;
     for i=1:length(labelingNames)
         outputFileName=[picPath(1:cutName-1) '_' labelingNames{i} '.tif'];
         outputFileNameCAdapted=[picPath(1:cutName-1) '_' labelingNames{i} '_ContrastAdapted.tif'];
         if strcmp(labelingNames{i},'DAPI')==1
             for j = 1:colorsCount:size(tiff_info, 1)
                 countDAPI=countDAPI+1;
                 % -- Save DAPI channel into a specific file (function imwrite)
                 figure(10);
                 temp_tiff= imread(picPath, j+colorsCount-1);
                 imwrite(imread(picPath, j+colorsCount-1), outputFileName, 'WriteMode', 'append');
                 % -- Adapt Contrast of the DAPI layer for easier segmentation
                 % (function adapthisteq) and save
                 subplot(1,2,1); imagesc(temp_tiff); axis image; colormap gray; title(labelingNames{i});
                 temp_contrast= adapthisteq(imread(picPath, j+colorsCount-1),'NumTiles',[100 100],'clipLimit',0.8,'Distribution','uniform');
                 subplot(1,2,2); imshow(temp_contrast); colormap gray; title([labelingNames{i} ' Contrast_Adapted'], 'Interpreter', 'none');
                 imwrite(temp_contrast, outputFileNameCAdapted, 'WriteMode', 'append');
                               %Layers.(labelingNames{i})(:,:,countDAPI)=temp_tiff; %temp_contrast;
             end
         else
             % -- Save vGlut2 channel into a specific file
             if strcmp(labelingNames{i},'vGlut2')==1
                 for j = 1:colorsCount:size(tiff_info, 1)
                     countvGlut2=countvGlut2+1;
                     temp_tiff= imread(picPath, j);
                     imwrite(imread(picPath, j), outputFileName, 'WriteMode', 'append');
                                %Layers.(labelingNames{i})(:,:,countvGlut2)=temp_tiff;
                     % -- Adapt Contrast of the DAPI layer for easier segmentation
                     % (function adapthisteq) and save
                     subplot(1,2,1); imagesc(temp_tiff); axis image; colormap gray; title(labelingNames{i});
                     temp_contrast= adapthisteq(imread(picPath, j+i-1),'NumTiles',[100 100],'clipLimit',0.8,'Distribution','uniform');
                     subplot(1,2,2); imshow(temp_contrast); colormap gray; title([labelingNames{i} ' Contrast_Adapted'], 'Interpreter', 'none');
                     imwrite(temp_contrast, outputFileNameCAdapted, 'WriteMode', 'append');
                     %Layers.(labelingNames{i})(:,:,countvGlut2)=temp_contrast;
                 end
                 % -- Save Gad1b channel into a specific file
             elseif strcmp(labelingNames{i},'Gad1b')==1
                 for j = 1:colorsCount:size(tiff_info, 1)
                     countGad1b=countGad1b+1;
                     temp_tiff= imread(picPath, j+1);
                     imwrite(imread(picPath, j+1), outputFileName, 'WriteMode', 'append');
                                 %Layers.(labelingNames{i})(:,:,countGad1b)=temp_tiff;
                     % -- Adapt Contrast of the DAPI layer for easier segmentation
                     % (function adapthisteq) and save
                     subplot(1,2,1); imagesc(temp_tiff); axis image; colormap gray; title(labelingNames{i});
                     temp_contrast= adapthisteq(imread(picPath, j+i-1),'NumTiles',[100 100],'clipLimit',0.8,'Distribution','uniform');
                     subplot(1,2,2); imshow(temp_contrast); colormap gray; title([labelingNames{i} ' Contrast_Adapted'], 'Interpreter', 'none');
                     imwrite(temp_contrast, outputFileNameCAdapted, 'WriteMode', 'append');
                     %Layers.(labelingNames{i})(:,:,countGad1b)=temp_contrast;
                 end
             end
         end
     end
end

%% ----MASK-----%
if MaskDevelopperMode==1
    %Masks/Save Masks
    Layers=struct();
    Border=struct();
    Data=struct();
    close all;
    outputFileNameCAdapted=[picPath(1:cutName-1) '_' labelingNames{3} '_ContrastAdapted.tif'];
    DAPIAdaptedPath=outputFileNameCAdapted; sizePic_1=size(imread(outputFileNameCAdapted),1);sizePic_2=size(imread(outputFileNameCAdapted),2);
    for maskCount=1:size(maskNames,1)              % if you want more than one mask you should change "1" to "2" and the mask will be the maskNames (also change down)
        Masks.(maskNames{maskCount}).tiff=zeros(sizePic_1,sizePic_2,Masks.(maskNames{maskCount}).Stop-Masks.(maskNames{maskCount}).Start);
        NumMas=[MaSt:Masks.(maskNames{maskCount}).Stop,MaSt:-1:Masks.(maskNames{maskCount}).Start];
        prevPos=cell(max(NumMas));
        
        for j = 1:length(NumMas)
            figure;
            imagesc(imread(DAPIAdaptedPath, NumMas(j))); axis image;colormap gray; title(maskNames{maskCount});
            if NumMas(j)==MaSt
                if isempty(prevPos{NumMas(j)})==1
                    polyhandle = impoly;
                else
                    clf
                end
            else
                polyhandle = impoly(gca,prevPos{NumMas(j-1)});
            end
            if isempty(prevPos{NumMas(j)})==1
                position = wait(polyhandle);
                clf
                prevPos{NumMas(j)}=position;
            end
        end
    end
    save(sprintf('%s%s%s%s',resultsFolder,'/',string(side{si}),'_poly_',FileIn,'.mat'),'prevPos');
end


if MaskDevelopperMode==2
    %Correct Masks/Save Masks
    close all;
    outputFileNameCAdapted=[picPath(1:cutName-1) '_' labelingNames{3} '_ContrastAdapted.tif'];
    DAPIAdaptedPath=outputFileNameCAdapted; sizePic_1=size(imread(outputFileNameCAdapted),1);sizePic_2=size(imread(outputFileNameCAdapted),2);
    load(sprintf('%s%s%s%s',resultsFolder,'/',string(side{si}),'_poly_',FileIn,'.mat'));
    for maskCount=1:size(maskNames,1)              % if you want more than one mask you should change "1" to "2" and the mask will be the maskNames (also change down)
        Masks.(maskNames{maskCount}).tiff=zeros(sizePic_1,sizePic_2,Masks.(maskNames{maskCount}).Stop-Masks.(maskNames{maskCount}).Start);
        for j = Masks.(maskNames{maskCount}).Start:Masks.(maskNames{maskCount}).Stop
            figure;
            imagesc(imread(DAPIAdaptedPath, j)); axis image;colormap gray; title(maskNames{maskCount});
            polyhandle = impoly(gca,prevPos{j});          
            %if isempty(prevPos{j})==1
                position = wait(polyhandle);
                clf
                prevPos{j}=position;
            %end
        end
    end
    save(sprintf('%s%s%s%s',resultsFolder,'/',string(side{si}),'_poly_',FileIn,'.mat'),'prevPos');
end

load(sprintf('%s%s%s%s',resultsFolder,'/',string(side{si}),'_poly_',FileIn,'.mat'));
%% ----NEW AXIS LINE-----%
if NewAxisLineDevelopperMode==1
    outputFileNameCAdapted=[picPath(1:cutName-1) '_' labelingNames{3} '_ContrastAdapted.tif'];
    DAPIAdaptedPath=outputFileNameCAdapted; sizePic_1=size(imread(outputFileNameCAdapted),1);sizePic_2=size(imread(outputFileNameCAdapted),2);
 
    MedianLine=cell(Masks.(maskNames{1}).Stop,1);
    for maskCount=1:size(maskNames,1)                  % if you want more than one mask you should change "1" to "2" and the mask will be the maskNames (also change top)
        for j = Masks.(maskNames{1}).Start:Masks.(maskNames{1}).Stop
            
        % -- Images
        Layers.(labelingNames{1})=imread(picPath, j*colorsCount-2);
        Layers.(labelingNames{2})=imread(picPath, j*colorsCount-1);
        Layers.(labelingNames{3})=imread(picPath, j*colorsCount);
        
        % -- Apply Mask
        Masks.(maskNames{maskCount}).tiff(:,:,j) = poly2mask(prevPos{j}(:,1),prevPos{j}(:,2),sizePic_1,sizePic_2);
        Masks.(maskNames{maskCount}).tiff(1,:,j)=0;
        Masks.(maskNames{maskCount}).tiff(sizePic_1,:,j)=0;
        Masks.(maskNames{maskCount}).tiff(:,1,j)=0;
        Masks.(maskNames{maskCount}).tiff(:,sizePic_2,j)=0;
        Border.(maskNames{maskCount}).tiff(:,:,j)= boundarymask(Masks.(maskNames{maskCount}).tiff(:,:,j));
        %Masks.(maskNames{maskCount}).MaskData=zeros(size(temp_tiff,1),size(temp_tiff,2),Masks.Tectum.Stop-Masks.Tectum.Start);
        temp_tiff= imread(DAPIAdaptedPath, j);
        Data.(maskNames{maskCount}).(labelingNames{3}).Tiff=zeros(sizePic_1,sizePic_2,Masks.(maskNames{maskCount}).Stop-Masks.(maskNames{maskCount}).Start);
        temp_tiff(~Masks.(maskNames{maskCount}).tiff(:,:,j))= 0; Data.(maskNames{maskCount}).(labelingNames{3}).Tiff(:,:,j)=temp_tiff;
        MedianLine{j}=MakeNewAxisLine(temp_tiff);
%         save(sprintf('%s%s%s%s',resultsFolder,'/',string(side{si}),'_line_',FileIn,'.mat'),'MedianLine');
        end
    end
        
    save(sprintf('%s%s%s%s',resultsFolder,'/MeanLine/',string(side{si}),'_line_',FileIn,'.mat'),'MedianLine');
end


load(sprintf('%s%s%s%s',resultsFolder,'/MeanLine/',string(side{si}),'_line_',FileIn,'.mat'));

%% ---- NORMALISATION ---- %%
if Normalization_Parameter==1
    for i=1:length(labelingNames)   
        
        outputFileNameCAdapted=[picPath(1:cutName-1) '_' labelingNames{i} '_ContrastAdapted.tif'];
        DAPIAdaptedPath=outputFileNameCAdapted; sizePic_1=size(imread(outputFileNameCAdapted),1);sizePic_2=size(imread(outputFileNameCAdapted),2);
        
        for maskCount=1:size(maskNames,1)                  % if you want more than one mask you should change "1" to "2" and the mask will be the maskNames (also change top)
            for j = Masks.(maskNames{maskCount}).Start:Masks.(maskNames{maskCount}).Stop
                
                % -- Apply Mask
                Masks.(maskNames{maskCount}).tiff(:,:,j) = poly2mask(prevPos{j}(:,1),prevPos{j}(:,2),sizePic_1,sizePic_2);
                Masks.(maskNames{maskCount}).tiff(1,:,j)=0;
                Masks.(maskNames{maskCount}).tiff(sizePic_1,:,j)=0;
                Masks.(maskNames{maskCount}).tiff(:,1,j)=0;
                Masks.(maskNames{maskCount}).tiff(:,sizePic_2,j)=0;
                Border.(maskNames{maskCount}).tiff(:,:,j)= boundarymask(Masks.(maskNames{maskCount}).tiff(:,:,j));
                
                temp_tiff= imread(DAPIAdaptedPath, j);
                Data.(maskNames{maskCount}).(labelingNames{i}).Tiff=zeros(sizePic_1,sizePic_2,Masks.(maskNames{maskCount}).Stop-Masks.(maskNames{maskCount}).Start);
                temp_tiff(~Masks.(maskNames{maskCount}).tiff(:,:,j))= 0; contrast_tiff_matrix(:,:,j)=temp_tiff;
                
                
                temp_tiff=  imread(picPath, j*colorsCount-3+i);
                Data.(maskNames{maskCount}).(labelingNames{i}).Tiff=zeros(sizePic_1,sizePic_2,Masks.(maskNames{maskCount}).Stop-Masks.(maskNames{maskCount}).Start);
                temp_tiff(~Masks.(maskNames{maskCount}).tiff(:,:,j))= 0; image_tiff_matrix{i}(:,:,j)=temp_tiff;
                
                
                
            end
        end
        
        Contrast_average(i)=mean(contrast_tiff_matrix(find(contrast_tiff_matrix)));
        Image_average(i)=mean(image_tiff_matrix{i}(find(image_tiff_matrix{i})));
    end
              
end


%% ----SEGMENTATION-----%
if SegmentationMode==1


%Adata Contrast/Perform Segmentation/Check results
for i=1:length(labelingNames)   
    outputFileNameCAdapted=[picPath(1:cutName-1) '_' labelingNames{i} '_ContrastAdapted.tif'];
    DAPIAdaptedPath=outputFileNameCAdapted; sizePic_1=size(imread(outputFileNameCAdapted),1);sizePic_2=size(imread(outputFileNameCAdapted),2);
 
    for maskCount=1:size(maskNames,1)                  % if you want more than one mask you should change "1" to "2" and the mask will be the maskNames (also change top)
        for j = Masks.(maskNames{maskCount}).Start:Masks.(maskNames{maskCount}).Stop
            
            % -- Images
            Layers.(labelingNames{1})=imread(picPath, j*colorsCount-2);
            Layers.(labelingNames{2})=imread(picPath, j*colorsCount-1);
            Layers.(labelingNames{3})=imread(picPath, j*colorsCount);
        
            % -- Apply Mask
            Masks.(maskNames{maskCount}).tiff(:,:,j) = poly2mask(prevPos{j}(:,1),prevPos{j}(:,2),sizePic_1,sizePic_2);
            Masks.(maskNames{maskCount}).tiff(1,:,j)=0;
            Masks.(maskNames{maskCount}).tiff(sizePic_1,:,j)=0;
            Masks.(maskNames{maskCount}).tiff(:,1,j)=0;
            Masks.(maskNames{maskCount}).tiff(:,sizePic_2,j)=0;
            Border.(maskNames{maskCount}).tiff(:,:,j)= boundarymask(Masks.(maskNames{maskCount}).tiff(:,:,j));
            %Masks.(maskNames{maskCount}).MaskData=zeros(size(temp_tiff,1),size(temp_tiff,2),Masks.Tectum.Stop-Masks.Tectum.Start);
            Mask_Area(j,:)=[j-Tectum_Begin{ANI}, sum(sum(Masks.(maskNames{maskCount}).tiff(:,:,j)))];
            temp_tiff= imread(DAPIAdaptedPath, j);
            Data.(maskNames{maskCount}).(labelingNames{i}).Tiff=zeros(sizePic_1,sizePic_2,Masks.(maskNames{maskCount}).Stop-Masks.(maskNames{maskCount}).Start);
            temp_tiff(~Masks.(maskNames{maskCount}).tiff(:,:,j))= 0; Data.(maskNames{maskCount}).(labelingNames{i}).Tiff(:,:,j)=temp_tiff;
        
            % -- Segmentation
            Cmatrix =  normxcorr2(double(gauss),double(temp_tiff));                                   % Cmatrix =  normxcorr2(double(gauss),double(Masks.(maskNames{maskCount}).Img(:,:,j)));        % Perform correlation between the image and the previously defined gaussian with size "taille"
            f = imfilter(Cmatrix, h);                                                                 % Filtered this 2D correlation images with the gaussian beam "h"
            bw= im2bw(f, cut.(labelingNames{i})*graythresh(f)*(Normalize_contrast_average(i)/Contrast_average(i)));                                                          % Build a binary image taking into account the noise of f and a given pourcentage of it
            CC = bwconncomp(bw, 8);
            L = labelmatrix(CC);
            S = regionprops(CC, 'Area');
            bw2 = ismember(L, find([S.Area] <= P.(labelingNames{i})));                     % Keep all the regions with an area smaller than P
            bw3 = imfill(bw2,'holes');                                  % Fill the holes in the regions
            bw4 = bwareaopen(bw3, open.(labelingNames{i}));                                % Remove region with area smaller than parameter "open"
        
            % -- Check segmentation efficiency
% % % % % %             if DevelopperMode==1
% % % % % %                 figure(20); subplot(221); imshow(Data.(maskNames{maskCount}).(labelingNames{i}).Tiff(:,:,j),[])
% % % % % %                 subplot(222); imshow(Cmatrix,[]);
% % % % % %                 subplot(223); imshow(f,[]);
% % % % % %                 subplot(224); imshow(bw4,[]);
% % % % % %             end
% % % % % %         
            % -- Check centroids
            CC = bwconncomp(bw4,8); L = labelmatrix(CC);
             Data.(maskNames{maskCount}).(labelingNames{i}).Seg{j}=L;
            A= regionprops(L,'Centroid','Area'); dd=(taille-1)/2; % determine the number, the centroid and the area of the regions
% % % % % %             clf; clear ax; figure(32);
% % % % % %             ax=struct();
% % % % % %             ax(1)=subplot(2,3,1); imshow(Layers.(labelingNames{1})), hold on; title(labelingNames{1})
% % % % % %             ax(2)=subplot(2,3,2); imshow(Layers.(labelingNames{2})), hold on; title(labelingNames{2})
            %Merge Green and Red channels
            C = imfuse(Layers.(labelingNames{1}),Layers.(labelingNames{2}),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
% % % % % %             ax(3)=subplot(2,3,3); imshow(imread(DAPIAdaptedPath, j)), hold on; %title(labelingNames{3})
        
            % -- Get/plot patches where fluorescence will be inspected
            squarePatch=zeros(numel(A),4);
            Masks.projections{j}=zeros(numel(A),7);
      
            for k=1:numel(A)
% % % % % %                 for plotIdx=1:length(ax)
% % % % % %                     plot(ax(plotIdx),A(k).Centroid(1)-dd,A(k).Centroid(2)-dd, 'r*');
% % % % % %                 end
                centerX=A(k).Centroid(1)-dd; centerY=A(k).Centroid(2)-dd;
                squarePatch(k,:)=[centerX-2 centerY-2 4 4];  %%% [centerX-2 centerY-2.5 4 4];
                idx=knnsearch(MedianLine{j}{:},[centerY centerX]); %%%  squarePatch(k,1:2)
                Masks.projections{j}(k,1)=sum(sqrt(diff(MedianLine{j}{:}(1:idx,1)).^2+diff(MedianLine{j}{:}(1:idx,2)).^2))/sum(sqrt(diff(MedianLine{j}{:}(:,1)).^2+diff(MedianLine{j}{:}(:,2)).^2));
                Masks.projections{j}(k,2)=j-Tectum_Begin{ANI};
                Masks.projections{j}(k,3:4)=[centerX centerY];
                Masks.projections{j}(k,5:6)=MedianLine{j}{:}(idx,:);
                Masks.projections{j}(k,7)=sign(dot([centerX centerY]-MedianLine{j}{:}(idx,[2,1]),MedianLine{j}{:}(idx,[2 1])-[XPoint.(side{si}){ANI} YPoint.(side{si}){ANI}]))*sqrt([centerX-MedianLine{j}{:}(idx,2)].^2 + (centerY-MedianLine{j}{:}(idx,1)).^2);
 
            
% % %                 rectangle('Position',squarePatch(k,:),'EdgeColor','b');
                for colorLayers=1:length(labelingNames)
                    valsTemp = mean2(imcrop(Layers.(labelingNames{colorLayers}),squarePatch(k,:)))*(Normalize_image_average(colorLayers)/Image_average(colorLayers));
                    SDvalsTemp = std2(imcrop(Layers.(labelingNames{colorLayers}),squarePatch(k,:)));
                    %labelingNames={'vGlut2','Gad1b','DAPI'};
                    Data.(maskNames{maskCount}).(labelingNames{i}).Layer{j}(k,colorLayers)=valsTemp;                           % Mean fluorescence value for layer for marker
                    Data.(maskNames{maskCount}).(labelingNames{i}).SDLayer{j}(k,colorLayers)=SDvalsTemp;
                end
            end
            
            
            % -- Plot color patches on DAPI and Merged corresponding to the detected fluorescence
% % % % % %             ax(5)=subplot(2,3,5); imshow(C), hold on
% % % % % %             plotPatchesE(A,Data,j,dd,maskNames{maskCount},labelingNames{i}); 
% % % % % %             ax(6)=subplot(2,3,6); imshow(imread(DAPIAdaptedPath, j)), hold on; %title(labelingNames{3})
% % % % % %             plotPatchesE(A,Data,j,dd,maskNames{maskCount},labelingNames{i})
% % % % % %             linkaxes([ax(1),ax(2),ax(3),ax(4),ax(5),ax(6)],'xy'); %Enables to zoom on all subplots simultaneously. 
        end
    end
end

%Export Figure
% % % % % % export_fig(fullfile(resultsFolder,['ColorPatch '  date]),'-eps','-painters');
%% ---- COLOR DETECTION PLOT-----%
%%%  Remove elements in the border  %%%
ifB=struct();
% aDapi=struct();
% aGad=struct();
% aGlu=struct;
ifG=struct();
ifGlu=struct();
% aGl=struct();
ifGad=struct();
% aGa=struct();
for j = Masks.(maskNames{maskCount}).Start:Masks.(maskNames{maskCount}).Stop
dapi=Data.(maskNames{maskCount}).(labelingNames{3}).Seg{j};
glu=Data.(maskNames{maskCount}).(labelingNames{1}).Seg{j};
gad=Data.(maskNames{maskCount}).(labelingNames{2}).Seg{j};
bord=Border.(maskNames{maskCount}).tiff(:,:,j);
pixelDapi=regionprops(dapi,'PixelList');
pixelglu=regionprops(glu,'PixelList');
pixelgad=regionprops(gad,'PixelList');
gad(gad>0)=1;
glu(glu>0)=1;
dapi(dapi>0)=1;

ifB.dapi{j}=zeros(length(pixelDapi),1);
 for k=1:length(pixelDapi)
     pix_DAPI=[pixelDapi(k).PixelList(:,2), pixelDapi(k).PixelList(:,1)];
     pix_DAPI(any(pix_DAPI>=1025,2),:) = [];
     aDapi=mean(diag(bord(pix_DAPI)));
    if(aDapi>0)
       ifB.dapi{j}(k)=1;
    end  
 end
 
 ifB.gad{j}=zeros(length(pixelgad),1);
  for k=1:length(pixelgad)
      pix_Gad=[pixelgad(k).PixelList(:,2),pixelgad(k).PixelList(:,1)];
      pix_Gad(any(pix_Gad>=1025,2),:) = [];  
      aGad=mean(diag(bord(pix_Gad)));
    if(aGad>0)
         ifB.gad{j}(k)=1;
    end  
  end
  
  ifB.glu{j}=zeros(length(pixelglu),1);
  for k=1:length(pixelglu)
      pix_Glut=[pixelglu(k).PixelList(:,2),pixelglu(k).PixelList(:,1)];
      pix_Glut(any(pix_Glut>=1025,2),:) = [];
      aGlu=mean(diag(bord(pix_Glut)));
    if(aGlu>0)
        ifB.glu{j}(k)=1;
    end  
  end
 
end


for j = Masks.(maskNames{maskCount}).Start:Masks.(maskNames{maskCount}).Stop
dapi=Data.(maskNames{maskCount}).(labelingNames{3}).Seg{j};
glu=Data.(maskNames{maskCount}).(labelingNames{1}).Seg{j};
gad=Data.(maskNames{maskCount}).(labelingNames{2}).Seg{j};
pixelDapi=regionprops(dapi,'PixelList');
gad(gad>0)=1;
glu(glu>0)=1;
dga=0;
dgl=0;
ifG.Leyer{j}=zeros(length(pixelDapi),2);
 for k=1:length(pixelDapi)
    aGad(k)=mean(diag(gad(pixelDapi(k).PixelList(:,2),pixelDapi(k).PixelList(:,1))));
    if(aGad(k)>overlap) %0.55
        ifG.Leyer{j}(k,2)=1;
        dga=dga+1;
    end  
    aGlu=mean(diag(glu(pixelDapi(k).PixelList(:,2),pixelDapi(k).PixelList(:,1))));
    if(aGlu>overlap) %0.55
        ifG.Leyer{j}(k,1)=1;
        dgl=dgl+1;
    end  
 end
 
% ifG.Leyer{j}([ifB.dapi{j}==1],:)=NaN;
dapi(dapi>0)=1;
X = cat(3,glu,gad,dapi);
% % % % % % % figure
% % % % % % % imshow(double(X))
end

ver=[ifG.Leyer{1},Data.(maskNames{maskCount}).(labelingNames{3}).Layer{1}];

glu(glu>0)=1;
gad(gad>0)=1;
dapi(dapi>0)=1;
X = cat(3,glu,gad,dapi);
% figure
% imshow(double(X))


%%% count glu gad without dapi  %%%%

for j = Masks.(maskNames{maskCount}).Start:Masks.(maskNames{maskCount}).Stop
dapi=Data.(maskNames{maskCount}).(labelingNames{3}).Seg{j};
glu=Data.(maskNames{maskCount}).(labelingNames{1}).Seg{j};
gad=Data.(maskNames{maskCount}).(labelingNames{2}).Seg{j};
pixelGlu=regionprops(glu,'PixelList');
pixelGad=regionprops(gad,'PixelList');
dapi(dapi>0)=1;

gl=0;
ifGlu.Leyer{j}=zeros(length(pixelGlu),1);
 for k=1:length(pixelGlu)
    aGl=mean(diag(dapi(pixelGlu(k).PixelList(:,2),pixelGlu(k).PixelList(:,1))));
    if(aGl>0)
        ifGlu.Leyer{j}(k,1)=1;
        gl=gl+1;
    end  
 end
%  ifGlu.Leyer{j}([ifB.glu{j}==1],:)=NaN;
 
 ga=0;
 ifGad.Leyer{j}=zeros(length(pixelGad),1);
 for k=1:length(pixelGad)
    aGa=mean(diag(dapi(pixelGad(k).PixelList(:,2),pixelGad(k).PixelList(:,1))));
    if(aGa>0)
        ifGad.Leyer{j}(k,1)=1;
        ga=ga+1;
    end  
 end
%  ifGad.Leyer{j}([ifB.gad{j}==1],:)=NaN;
end

IsGad=vertcat(ifGad.Leyer{:});
IsGlu=vertcat(ifGlu.Leyer{:});

%How was the detection of Green/Red performed
% % % % % % figure; set(gcf,'color','w');
tempvGlut2=Data.(maskNames{maskCount}).(labelingNames{i}).Layer{j}(:,1); %vGlut2 fluorescence values for each cell
tempGad1b=Data.(maskNames{maskCount}).(labelingNames{i}).Layer{j}(:,2); %Gad1b fluo. values 
% % % % % % plot(tempvGlut2,tempGad1b,'.'); axis square; hold on;

%Find when fluo of one channel is twice higher than the other.  
vGlutPlus=find(tempvGlut2./tempGad1b>2); Gad1bPlus=find(tempGad1b./tempvGlut2>2); 

%Plot with associated color
% % % % % % plot(tempvGlut2(vGlutPlus),tempGad1b(vGlutPlus),'r.'); axis square; hold on;
% % % % % % plot(tempvGlut2(Gad1bPlus),tempGad1b(Gad1bPlus),'g.'); axis square; hold on;

% %Histogram
% for i = 1:3
% for j = Masks.(maskNames{maskCount}).Start:Masks.(maskNames{maskCount}).Stop
% 
%     figure; 
%     h=histogram(Data.(maskNames{maskCount}).(labelingNames{i}).Layer{j}(:,1),200); hold on;
%     h.FaceColor = 'r';
% %     set(gca,'xscale','log')
%     figure; 
%     h=histogram(Data.(maskNames{maskCount}).(labelingNames{i}).Layer{j}(:,2),200); hold on;
%     h.FaceColor = 'g';
% %     set(gca,'xscale','log')
%     figure; 
%     h=histogram(Data.(maskNames{maskCount}).(labelingNames{i}).Layer{j}(:,3),200); hold on;
%     h.FaceColor = 'b';
% %     set(gca,'xscale','log')
% end
% end


% %Histogram
% for i = 1:3
%     for j = Masks.(maskNames{maskCount}).Start:Masks.(maskNames{maskCount}).Stop
%           vGlut2=Data.(maskNames{maskCount}).(labelingNames{i}).Layer{j}(:,1);
%           Gad1b=Data.(maskNames{maskCount}).(labelingNames{i}).Layer{j}(:,2);
%           DAPI=Data.(maskNames{maskCount}).(labelingNames{i}).Layer{j}(:,3);
%           figure;
%           scatterhist(vGlut2, Gad1b,'NBins',[300,300],'Style','stairs','Marker','.','MarkerSize',7); 
%           title(sprintf('%s%d%s%s','frame=',j,'  segmentation=',labelingNames{i}));
%           hold on;
% 
%           figure; 
%           scatterhist(vGlut2, DAPI,'NBins',[300,300],'Style','stairs','Marker','.','MarkerSize',7); 
%           title(sprintf('%s%d%s%s','frame=',j,'  segmentation=',labelingNames{i}));
%           hold on;
% 
%           figure; 
%          scatterhist(Gad1b, DAPI,'NBins',[300,300],'Style','stairs','Marker','.','MarkerSize',7); 
%          title(sprintf('%s%d%s%s','frame=',j,'  segmentation=',labelingNames{i}));
%          hold on;
% 
%     end
%     TODO=vertcat(Data.(maskNames{maskCount}).(labelingNames{i}).Layer{:});
%     vGlut2=TODO(:,1);
%     Gad1b=TODO(:,2);
%     DAPI=TODO(:,3);
%     figure;
%           scatterhist(vGlut2, Gad1b,'NBins',[300,300],'Style','stairs','Marker','.','MarkerSize',7); 
%           title(sprintf('%s%s%s%s','frame=','all','  segmentation=',labelingNames{i}));
%           hold on;
% 
%           figure; 
%           scatterhist(vGlut2, DAPI,'NBins',[300,300],'Style','stairs','Marker','.','MarkerSize',7); 
%           title(sprintf('%s%s%s%s','frame=','all','  segmentation=',labelingNames{i}));
%           hold on;
% 
%           figure; 
%          scatterhist(Gad1b, DAPI,'NBins',[300,300],'Style','stairs','Marker','.','MarkerSize',7); 
%          title(sprintf('%s%s%s%s','frame=','all','  segmentation=',labelingNames{i}));
%          hold on;
% end

VerTodo=[vertcat(ifB.dapi{:}),vertcat(ifG.Leyer{:}),vertcat(Data.(maskNames{maskCount}).(labelingNames{3}).Layer{:}),vertcat(Masks.projections{:})]; 
VerTodo=VerTodo(~(VerTodo(:,1)==1),:);
    vGlut2=VerTodo(:,4);
    Gad1b=VerTodo(:,5);
    DAPI=VerTodo(:,6);
%     figure;
%           scatterhist(vGlut2(VerTodo(:,2)==1), Gad1b(VerTodo(:,2)==1),'NBins',[300,300],'Style','stairs','Marker','.','MarkerSize',7); 
%           title(sprintf('%s%s%s%s','frame=','all','  segmentation=',labelingNames{i}));
%           hold on;
% 
%           figure; 
%           scatterhist(vGlut2(VerTodo(:,2)==1), DAPI(VerTodo(:,2)==1),'NBins',[300,300],'Style','stairs','Marker','.','MarkerSize',7); 
%           title(sprintf('%s%s%s%s','frame=','all','  segmentation=',labelingNames{i}));
%           hold on;
% 
%           figure; 
%          scatterhist(Gad1b(VerTodo(:,2)==1), DAPI(VerTodo(:,2)==1),'NBins',[300,300],'Style','stairs','Marker','.','MarkerSize',7); 
%          title(sprintf('%s%s%s%s','frame=','all','  segmentation=',labelingNames{i}));
%          hold on;
%     figure;
%           scatterhist(vGlut2(VerTodo(:,3)==1), Gad1b(VerTodo(:,3)==1),'NBins',[300,300],'Style','stairs','Marker','.','MarkerSize',7); 
%           title(sprintf('%s%s%s%s','frame=','all','  segmentation=',labelingNames{i}));
%           hold on;
% 
%           figure; 
%           scatterhist(vGlut2(VerTodo(:,3)==1), DAPI(VerTodo(:,3)==1),'NBins',[300,300],'Style','stairs','Marker','.','MarkerSize',7); 
%           title(sprintf('%s%s%s%s','frame=','all','  segmentation=',labelingNames{i}));
%           hold on;
% 
%           figure; 
%          scatterhist(Gad1b(VerTodo(:,3)==1), DAPI(VerTodo(:,3)==1),'NBins',[300,300],'Style','stairs','Marker','.','MarkerSize',7); 
%          title(sprintf('%s%s%s%s','frame=','all','  segmentation=',labelingNames{i}));
%          hold on;

         
VerTodo=VerTodo(~(VerTodo(:,6)<70 & VerTodo(:,5)<70 & VerTodo(:,4)<70),:);         
gru=zeros(size(VerTodo,1),1);
gru(VerTodo(:,2)==1 & VerTodo(:,3)==0)=1;
gru(VerTodo(:,2)==0 & VerTodo(:,3)==1)=2;
gru(VerTodo(:,2)==0 & VerTodo(:,3)==0)=3;
gru(VerTodo(:,2)==1 & VerTodo(:,3)==1 & VerTodo(:,4)>1.5*VerTodo(:,5))=1;
gru(VerTodo(:,2)==1 & VerTodo(:,3)==1 & 1.5*VerTodo(:,4)<VerTodo(:,5))=2;
gru(gru==0)=3;

gru(gru==3 & VerTodo(:,6)<70 & VerTodo(:,4)>1.5*VerTodo(:,5))=1;
gru(gru==3 & VerTodo(:,6)<70 & 1.5*VerTodo(:,4)<VerTodo(:,5))=2;

gru(gru==1 & VerTodo(:,4)<80)=3;
gru(gru==2 & VerTodo(:,5)<70)=3;


VerTodo=[VerTodo,gru];
VerTodo=VerTodo(~(gru==3 & VerTodo(:,6)<70),:);
color=VerTodo(:,14);
color(VerTodo(:,14)==1)=3;
color(VerTodo(:,14)==2)=2;
color(VerTodo(:,14)==3)=1;

% % % % % % figure;
% % % % % % scatterhist(VerTodo(:,4), VerTodo(:,5),'Group',color,'Color', 'bgr','NBins',[400,700],'Style','stairs','Marker','.','MarkerSize',7); %VerTodo(:,7)
% % % % % % xlabel('vGlut2 actuviy');
% % % % % % ylabel('Gad1b activity');
% % % % % % title(sprintf('%s%s','color cell distribution',' (Glut vs Gad)'));
% % % % % % hold on;


color2=VerTodo(:,14);
color2(:)=1;
color2(VerTodo(:,5)>2*VerTodo(:,4))=3;
color2(VerTodo(:,5)<VerTodo(:,4)/2)=2;

% % % % % % figure;
% % % % % % scatterhist(VerTodo(:,4), VerTodo(:,5),'Group',color2,'Color', 'brg','NBins',[400,700],'Style','stairs','Marker','.','MarkerSize',7); %VerTodo(:,7)
% % % % % % xlabel('vGlut2 actuviy');
% % % % % % ylabel('Gad1b activity');
% % % % % % title(sprintf('%s%s','color cell distribution',' (Glut vs Gad)'));
% % % % % % hold on;


Mask_Area=Mask_Area(find(Mask_Area(:,2)),:);

IsGlu=IsGlu(~(vertcat(ifB.glu{:})==1));
IsGad=IsGad(~(vertcat(ifB.gad{:})==1));

Numbers = {numel(VerTodo(VerTodo(:,14)==1));numel(VerTodo(VerTodo(:,14)==2));numel(VerTodo(VerTodo(:,14)==3));numel(VerTodo(:,14));numel(IsGlu(IsGlu==0));numel(IsGad(IsGad==0))};
Rows = {'vGlut2';'Gad1b';'DAPI';'Total Cells';'vGlut2 no DAPI';'Gad1b no DAPI'};
T = table(Numbers,'Rownames',Rows)
% writetable(T,sprintf('%s%s%s%s',resultsFolder,'/',string(side{si}),'_Table_',FileIn,'.txt'),'WriteRowNames',true)
% save(sprintf('%s%s%s%s',resultsFolder,'/',string(side{si}),'_VerTodo_',FileIn,'.mat'),'VerTodo','-v7.3')

writetable(T,sprintf('%s%s%s%s',resultsFolder,'/Comparation/Overlap/',string(overlap),'/',string(side{si}),'_Table_',FileIn,'.txt'),'WriteRowNames',true);
save(sprintf('%s%s%s%s',resultsFolder,'/Comparation/Overlap/',string(overlap),'/',string(side{si}),'_VerTodo_',FileIn,'.mat'),'VerTodo','-v7.3');
save(sprintf('%s%s%s%s',resultsFolder,'/Comparation/Overlap/',string(overlap),'/',string(side{si}),'_Mask_Area_',FileIn,'.mat'),'Mask_Area','-v7.3');

          
% for j = Masks.(maskNames{maskCount}).Start:Masks.(maskNames{maskCount}).Stop
% 
%     figure; 
%     h=histogram(Data.(maskNames{maskCount}).(labelingNames{i}).SDLayer{j}(:,1)); hold on;
%     h.FaceColor = 'r';
%     h=histogram(Data.(maskNames{maskCount}).(labelingNames{i}).SDLayer{j}(:,2)); hold on;
%     h.FaceColor = 'g';
%     h=histogram(Data.(maskNames{maskCount}).(labelingNames{i}).SDLayer{j}(:,3)); hold on;
%     h.FaceColor = 'b';
% end
% TODO=vertcat(Data.(maskNames{maskCount}).(labelingNames{i}).Layer{:});
%  figure; 
%     h=histogram(TODO(:,1)); hold on;
%     h.FaceColor = 'r';
% %     set(gca,'xscale','log')
%  figure;    
%     h=histogram(TODO(:,2)); hold on;
%     h.FaceColor = 'g';
% %     set(gca,'xscale','log')
%  figure; 
%     h=histogram(TODO(:,3)); hold on;
%     h.FaceColor = 'b';
% %     set(gca,'xscale','log')
% 
% TODOSD=vertcat(Data.(maskNames{maskCount}).(labelingNames{i}).SDLayer{:});
%  figure; 
%     h=histogram(TODOSD(:,1)); hold on;
%     h.FaceColor = 'r';
% %     set(gca,'xscale','log')
%  figure;    
%     h=histogram(TODOSD(:,2)); hold on;
%     h.FaceColor = 'g';
% %     set(gca,'xscale','log')
%  figure; 
%     h=histogram(TODOSD(:,3)); hold on;
%     h.FaceColor = 'b';
% %     set(gca,'xscale','log')
% 
% %How was the detection of Green/Red performed
% figure; set(gcf,'color','w');
% tempvGlut2=TODO(:,1); %vGlut2 fluorescence values for each cell
% tempGad1b=TODO(:,2); %Gad1b fluo. values 
% plot(tempvGlut2,tempGad1b,'.'); axis square; hold on;
% % set(gca,'xscale','log')
% % set(gca,'yscale','log')
% %Find when fluo of one channel is twice higher than the other.  
% vGlutPlus=find(tempvGlut2./tempGad1b>2); Gad1bPlus=find(tempGad1b./tempvGlut2>2); 
% 
% %Plot with associated color
% plot(tempvGlut2(vGlutPlus),tempGad1b(vGlutPlus),'r.'); axis square; hold on;
% plot(tempvGlut2(Gad1bPlus),tempGad1b(Gad1bPlus),'g.'); axis square; hold on;    
% 
% 
%  figure; 
%     h=histogram((TODO(:,1)-(TODO(:,2))./(TODO(:,1))+(TODO(:,2)))); hold on;
%     h.FaceColor = 'r';
%     h=histogram(TODO(:,2)./TODO(:,1)); hold on;
%     h.FaceColor = 'g';
%     
%     
% 
% TODO=vertcat(Data.(maskNames{maskCount}).(labelingNames{1}).Layer{:});
% figure; 
%     h=histogram(TODO(:,1)); hold on;
%     h.FaceColor = 'r';
% figure; 
%     h=histogram(TODO(:,2)); hold on;
%     h.FaceColor = 'g';
% figure; 
%     h=histogram(TODO(:,3)); hold on;
%     h.FaceColor = 'b';
% 
% 
% indx=TODO(:,3)>65;
% TODO=TODO(indx,:);
% 
% %How was the detection of Green/Red performed
% figure; set(gcf,'color','w');
% tempvGlut2=TODO(:,1); %vGlut2 fluorescence values for each cell
% tempGad1b=TODO(:,2); %Gad1b fluo. values 
% tempDAPI=TODO(:,3); %DAPI fluo. values 
% plot(tempvGlut2,tempGad1b,'.'); axis square; hold on;
% % set(gca,'xscale','log')
% % set(gca,'yscale','log')
% %Find when fluo of one channel is twice higher than the other.  
% vGlutPlus=find(tempvGlut2./tempGad1b>2); Gad1bPlus=find(tempGad1b./tempvGlut2>2); 
% 
% %Plot with associated color
% plot(tempvGlut2(vGlutPlus),tempGad1b(vGlutPlus),'r.'); axis square; hold on;
% plot(tempvGlut2(Gad1bPlus),tempGad1b(Gad1bPlus),'g.'); axis square; hold on;    
%     
% 
% 
% figure; set(gcf,'color','w');
% plot(tempvGlut2,tempDAPI,'.'); axis square; hold on;
% 
% %Plot with associated color
% plot(tempvGlut2(vGlutPlus),tempDAPI(vGlutPlus),'r.'); axis square; hold on;
% plot(tempvGlut2(Gad1bPlus),tempDAPI(Gad1bPlus),'g.'); axis square; hold on;   
% 
% 
% figure; set(gcf,'color','w');
% plot(tempGad1b,tempDAPI,'.'); axis square; hold on;
% 
% %Plot with associated color
% plot(tempGad1b(vGlutPlus),tempDAPI(vGlutPlus),'r.'); axis square; hold on;
% plot(tempGad1b(Gad1bPlus),tempDAPI(Gad1bPlus),'g.'); axis square; hold on;   
% 
%     
% 
% 
% TODO=vertcat(Data.(maskNames{maskCount}).(labelingNames{2}).Layer{:});
% figure; 
%     h=histogram(TODO(:,1)); hold on;
%     h.FaceColor = 'r';
% figure; 
%     h=histogram(TODO(:,2)); hold on;
%     h.FaceColor = 'g';
% figure; 
%     h=histogram(TODO(:,3)); hold on;
%     h.FaceColor = 'b';
%     
%     
% indx=TODO(:,3)>65;
% TODO=TODO(indx,:);
% 
% %How was the detection of Green/Red performed
% figure; set(gcf,'color','w');
% tempvGlut2=TODO(:,1); %vGlut2 fluorescence values for each cell
% tempGad1b=TODO(:,2); %Gad1b fluo. values 
% tempDAPI=TODO(:,3); %DAPI fluo. values 
% plot(tempvGlut2,tempGad1b,'.'); axis square; hold on;
% % set(gca,'xscale','log')
% % set(gca,'yscale','log')
% %Find when fluo of one channel is twice higher than the other.  
% vGlutPlus=find(tempvGlut2./tempGad1b>2); Gad1bPlus=find(tempGad1b./tempvGlut2>2); 
% 
% %Plot with associated color
% plot(tempvGlut2(vGlutPlus),tempGad1b(vGlutPlus),'r.'); axis square; hold on;
% plot(tempvGlut2(Gad1bPlus),tempGad1b(Gad1bPlus),'g.'); axis square; hold on;    
% 
% figure; set(gcf,'color','w');
% plot(tempvGlut2,tempDAPI,'.'); axis square; hold on;
% 
% %Plot with associated color
% plot(tempvGlut2(vGlutPlus),tempDAPI(vGlutPlus),'r.'); axis square; hold on;
% plot(tempvGlut2(Gad1bPlus),tempDAPI(Gad1bPlus),'g.'); axis square; hold on;   
% 
% 
% figure; set(gcf,'color','w');
% plot(tempGad1b,tempDAPI,'.'); axis square; hold on;
% 
% %Plot with associated color
% plot(tempGad1b(vGlutPlus),tempDAPI(vGlutPlus),'r.'); axis square; hold on;
% plot(tempGad1b(Gad1bPlus),tempDAPI(Gad1bPlus),'g.'); axis square; hold on;   
% 
% 
% 
% TODO=vertcat(Data.(maskNames{maskCount}).(labelingNames{3}).Layer{:});
% figure; 
%     h=histogram(TODO(:,1)); hold on;
%     h.FaceColor = 'r';
% figure; 
%     h=histogram(TODO(:,2)); hold on;
%     h.FaceColor = 'g';
% figure; 
%     h=histogram(TODO(:,3)); hold on;
%     h.FaceColor = 'b';
% 
% 
% indx=TODO(:,3)>65;
% TODO=TODO(indx,:);
% 
% %How was the detection of Green/Red performed
% figure; set(gcf,'color','w');
% tempvGlut2=TODO(:,1); %vGlut2 fluorescence values for each cell
% tempGad1b=TODO(:,2); %Gad1b fluo. values 
% tempDAPI=TODO(:,3); %DAPI fluo. values 
% plot(tempvGlut2,tempGad1b,'.'); axis square; hold on;
% % set(gca,'xscale','log')
% % set(gca,'yscale','log')
% %Find when fluo of one channel is twice higher than the other.  
% vGlutPlus=find(tempvGlut2./tempGad1b>2); Gad1bPlus=find(tempGad1b./tempvGlut2>2); 
% 
% %Plot with associated color
% plot(tempvGlut2(vGlutPlus),tempGad1b(vGlutPlus),'r.'); axis square; hold on;
% plot(tempvGlut2(Gad1bPlus),tempGad1b(Gad1bPlus),'g.'); axis square; hold on;    
% 
% figure; set(gcf,'color','w');
% plot(tempvGlut2,tempDAPI,'.'); axis square; hold on;
% 
% %Plot with associated color
% plot(tempvGlut2(vGlutPlus),tempDAPI(vGlutPlus),'r.'); axis square; hold on;
% plot(tempvGlut2(Gad1bPlus),tempDAPI(Gad1bPlus),'g.'); axis square; hold on;   
% 
% 
% figure; set(gcf,'color','w');
% plot(tempGad1b,tempDAPI,'.'); axis square; hold on;
% 
% %Plot with associated color
% plot(tempGad1b(vGlutPlus),tempDAPI(vGlutPlus),'r.'); axis square; hold on;
% plot(tempGad1b(Gad1bPlus),tempDAPI(Gad1bPlus),'g.'); axis square; hold on;   
% 
% 
% plot3(TODO(:,1),TODO(:,2),TODO(:,3))
    
%Table with the number of cells
% Numbers = {numel(vGlutPlus);numel(Gad1bPlus);numel(tempGad1b)};
% Rows = {'vGlut2';'Gad1b';'Total Cells'};
% T = table(Numbers,'Rownames',Rows)
% writetable(T,'TableauResultatsCD1Part1.txt')

end

end

end
end