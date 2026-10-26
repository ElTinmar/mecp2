clear all; close all; 

resultsFolder='/media/emiliano/CODES/Segmentation/Resultats/Comparation/Overlap/';
overlap=[0.45];
FileInFolder=sprintf('%s%s%s',resultsFolder,string(overlap),'/Normalized_vGlu1.1-Gad0.7-DAPI1.1/');
files=dir(sprintf('%s%s',FileInFolder,'*VerTodo*.mat'));
area=dir(sprintf('%s%s',FileInFolder,'*Mask_Area*.mat'));

%%  Difine variables  %%
ALL=struct();
ALL.Animal=[];
ALL.Treatment =[];
ALL.Side=[];
ALL.Color=[];
ALL.Plane=[];
ALL.TectalLine=[];
ALL.IfGlu=[];
ALL.IfGad=[];
ALL.Glu=[];
ALL.Gad=[];
ALL.Dapi=[];
ALL.X=[];
ALL.Y=[];
ALL.YClosest=[];
ALL.XClosest=[];
ALL.Dist=[];

AREA=struct();
AREA.Area=[];
AREA.Plane=[];
AREA.Animal=[];
AREA.Treatment=[];

%%  load the data of segmented cells   %%
for file=1:length(files)
     load(sprintf('%s%s',FileInFolder, files(file).name));
     load(sprintf('%s%s',FileInFolder, area(file).name));
     tempString=strsplit(files(file).name,{' ','-','_'});
     side=string(tempString{1});
     animal=tempString{4};
     treatment=string(tempString{5});
     cells_number=size(VerTodo,1);
     
     ALL.Side=vertcat(ALL.Side, repmat(side,cells_number,1));  
     ALL.Animal= vertcat(ALL.Animal, repmat(animal,cells_number,1));
     ALL.Treatment = vertcat(ALL.Treatment, repmat(treatment,cells_number,1));  
     ALL.Color=vertcat(ALL.Color, VerTodo(:,14));
     ALL.Dist=vertcat(ALL.Dist, VerTodo(:,13));
     ALL.Plane=vertcat(ALL.Plane, VerTodo(:,8));
     ALL.TectalLine=vertcat(ALL.TectalLine, VerTodo(:,7));
     ALL.IfGlu=vertcat(ALL.IfGlu, VerTodo(:,2));
     ALL.IfGad=vertcat(ALL.IfGad, VerTodo(:,3));
     ALL.Glu=vertcat(ALL.Glu, VerTodo(:,4));
     ALL.Gad=vertcat(ALL.Gad, VerTodo(:,5));
     ALL.Dapi=vertcat(ALL.Dapi, VerTodo(:,6));
     ALL.X=vertcat(ALL.X, VerTodo(:,9));
     ALL.Y=vertcat(ALL.Y, VerTodo(:,10));
     ALL.YClosest=vertcat(ALL.YClosest, VerTodo(:,11));
     ALL.XClosest=vertcat(ALL.XClosest, VerTodo(:,12));
     
     Plane_Number= size(Mask_Area,1);
     AREA.Animal=vertcat(AREA.Animal, repmat(animal,Plane_Number,1));
     AREA.Treatment=vertcat(AREA.Treatment, repmat(treatment,Plane_Number,1));
     AREA.Area=vertcat(AREA.Area, Mask_Area(:,2));
     AREA.Plane=vertcat(AREA.Plane, Mask_Area(:,1));
  
end

%% pictures of animals wt 1, 2 and 3 are captured bottom-top and the rest upside down then should be changed  %%
NewPlane=ALL.Plane;
new_plane=AREA.Plane;
for an=1:3
    for pl=0:86
        NewPlane(ALL.Animal==num2str(an) & ALL.Treatment=="wt" & ALL.Plane==pl)=86-pl;
        new_plane(AREA.Animal==num2str(an) & AREA.Treatment=="wt" & AREA.Plane==pl)=86-pl;
    end
end
ALL.Plane=NewPlane;
AREA.Plane=new_plane;

%% joint together both tectal lines (right and left) continuously in one axes  ##
ALL.TectalLine(ALL.Side=='Rigth')=-ALL.TectalLine(ALL.Side=='Rigth')+1;
%% Create a table with all content  ##
ALL.Table=table(ALL.Animal, ALL.Color, ALL.Dapi, ALL.Dist, ALL.Gad, ALL.Glu, ALL.IfGad, ALL.IfGlu, ALL.Plane, ALL.Side, ALL.TectalLine, ALL.Treatment, ALL.X, ALL.XClosest, ALL.Y, ALL.YClosest);
ALL.Table.Properties.VariableNames={'animal','color','dapi','dist','gad', 'glu', 'ifGad','ifGlu','plane','side','tectalLine', 'treatment', 'X','closetX','Y','closetY'};

AREA.Table=table(AREA.Animal, AREA.Treatment, AREA.Plane, AREA.Area);
AREA.Table.Properties.VariableNames={'animal','treatment','plane', 'area'};

%%  --- Cell type number ---  %%

treatments=unique(ALL.Table.treatment);
animals=unique(ALL.Table.animal);
temp=1;
clear Traeatmenss Cells_Numel Ratio_stat Glu_stat Gad_stat
for trea=1:length(treatments)
    for ani=1:length(animals) 
        ArEa= 1;%sum(AREA.Table.area(AREA.Table.animal==animals(ani) & AREA.Table.plane>=planes(pla) & AREA.Table.plane<planes(pla+1) & AREA.Table.treatment==string(treatments(trea))));
        Cells_Numel(temp,:) = [str2num(animals(ani)), numel(ALL.Table.glu(ALL.Table.animal==animals(ani) & ALL.Table.treatment==string(treatments(trea)) & ALL.Table.color==1))/ArEa, numel(ALL.Table.gad(ALL.Table.animal==animals(ani) & ALL.Table.treatment==string(treatments(trea)) &ALL.Table.color==2))/ArEa, numel(ALL.Table.dapi(ALL.Table.animal==animals(ani) & ALL.Table.treatment==string(treatments(trea)) &ALL.Table.color==3))/ArEa];
        Traeatmenss(temp,1)=treatments(trea);
        temp=temp+1;
    end
end

Cells_Numel_Table=table(Traeatmenss, Cells_Numel(:,1),Cells_Numel(:,2),Cells_Numel(:,3),Cells_Numel(:,4), Cells_Numel(:,2)./Cells_Numel(:,3)); %(Cells_Numel(:,4)-Cells_Numel(:,3))./(Cells_Numel(:,4)+Cells_Numel(:,3))
Cells_Numel_Table.Properties.VariableNames={'treatment','animal','glu', 'gad', 'dapi', 'ratio'};
Cells_Numel_Table(isnan(Cells_Numel_Table.ratio),:)=[];


Ratio_stat=ranksum(Cells_Numel_Table.ratio(Cells_Numel_Table.treatment=="MeCP2"),Cells_Numel_Table.ratio(Cells_Numel_Table.treatment=="wt"));
Glu_stat=ranksum(Cells_Numel_Table.glu(Cells_Numel_Table.treatment=="MeCP2"),Cells_Numel_Table.glu(Cells_Numel_Table.treatment=="wt"));
Gad_stat=ranksum(Cells_Numel_Table.gad(Cells_Numel_Table.treatment=="MeCP2"),Cells_Numel_Table.gad(Cells_Numel_Table.treatment=="wt"));


Mean=varfun(@mean, Cells_Numel_Table,'InputVariables',{'glu', 'gad','ratio'},'GroupingVariables',{'treatment'});
STD=varfun(@std, Cells_Numel_Table,'InputVariables',{'glu', 'gad','ratio'},'GroupingVariables',{'treatment'});
ES=table(STD.std_glu./sqrt(STD.GroupCount), STD.std_gad./sqrt(STD.GroupCount), STD.std_ratio./sqrt(STD.GroupCount));
ES.Properties.VariableNames={'es_glu', 'es_gad', 'es_ratio'};
Datas=join(Mean,STD);
Datas=[Datas,ES];

figure
col=get(gca,'colororder');col(3:4,:,:)=[];
clear h
set(gca,'linewidth',1.5,'XColor',[0.3 0.3 0.3],'YColor',[0.3 0.3 0.3],'Box','on','TickLength',[0.005 0.005],'FontSize',16)
hold on

m= Datas.mean_ratio;
se= Datas.std_ratio;

errorbar(1:2,m,se,'.','color',[0.3 0.3 0.3],'LineWidth',2);
for i=1:2
hold on
bar(i,m(i),'Facecolor',col(i,:,:),'EdgeColor',[0.3 0.3 0.3],'LineWidth',2,'Barwidth',0.8);
end
axis square
xlim([0.5 2.5]);
ylim([0 1]);
% ht=text(1.5,0.6-0.6*0.05,['p = ' sprintf('%.2E',stat)],'VerticalAlignment','middle','HorizontalAlignment','center');
box on
set(gca,'xtick',[1,2],'xticklabel',[{'wild-type'},{'mecp2'}])
ylabel('Ratio (Vglut2a/Gad1b)','FontSize',16);
axis square



%%% stadistics
TfE = varfun(@sum, Cells_Numel_Table,'InputVariables',{'glu', 'gad','ratio'},'GroupingVariables',{'treatment','animal'});
TfE.sum_ratio = TfE.sum_glu./TfE.sum_gad;
tfe = [TfE.treatment(TfE.treatment=="MeCP2"), TfE.animal(TfE.treatment=="MeCP2"), repelem("gad",length(TfE.animal(TfE.treatment=="MeCP2")))', TfE.sum_gad(TfE.treatment=="MeCP2"); TfE.treatment(TfE.treatment=="MeCP2"), TfE.animal(TfE.treatment=="MeCP2"), repelem("glu",length(TfE.animal(TfE.treatment=="MeCP2")))', TfE.sum_glu(TfE.treatment=="MeCP2"); TfE.treatment(TfE.treatment=="wt"), TfE.animal(TfE.treatment=="wt")+8, repelem("gad",length(TfE.animal(TfE.treatment=="wt")))', TfE.sum_gad(TfE.treatment=="wt"); TfE.treatment(TfE.treatment=="wt"), TfE.animal(TfE.treatment=="wt")+8, repelem("glu",length(TfE.animal(TfE.treatment=="wt")))', TfE.sum_glu(TfE.treatment=="wt")];
% tfe(:,4) = str2double(tfe(:,4));
% glmfit(tfe(:,[1 3]),str2double(tfe(:,4)))

Line = TfE.treatment;
Ratio = TfE.sum_ratio;
BD = table(Line, Ratio);
glme = fitglme(BD, 'Ratio ~ 1 + Line','Distribution', 'Gamma', 'FitMethod','Laplace', 'DummyVarCoding','effects')



Fish = tfe(:,2);
Line = tfe(:,1);
NT = tfe(:,3);
Num = round(str2double(tfe(:,4)));
BD = table(Fish, Line, NT, Num);
glme = fitglme(BD,...
		'Num ~ NT * Line + (1|Fish)',...
		'Distribution','Poisson','Link','log','FitMethod','Laplace',...
		'DummyVarCoding','effects')

figure
plotResiduals(glme,'histogram','ResidualType','Pearson')
figure
plotResiduals(glme,'fitted','ResidualType','Pearson')
figure
plotResiduals(glme,'lagged','ResidualType','Pearson')

% Suponiendo que glme es tu modelo ajustado
H1 = [0 1 0 0]; % Comparacion del numero de GABAergicas en wt y MeCP2 (ns)
H2 = [0 1 0 1]; % Comparacion del numero de glutamatergicas en wt y MeCP2 (ns)
H3 = [0 0 1 0]; % Comparacion del numero de GABAergicas y glutamatergicas en MeCP2
H4 = [0 0 1 -1]; % Comparacion del numero de GABAergicas y glutamatergicas en wt
[pValue1 FTest1,] = coefTest(glme, H1) %NS
[pValue2, FTest2] = coefTest(glme, H2) %NS
[pValue3, FTest3] = coefTest(glme, H3) %*
[pValue4, FTest4] = coefTest(glme, H4) %*

% Correccion de Bonferrroni-Holms (super innecesario para tus valores)
pValue1c=pValue1-0.05/(4+1-3)
pValue2c=pValue2-0.05/(4+1-4)
pValue3c=pValue3-0.05/(4+1-1)
pValue4c=pValue4-0.05/(4+1-2)





figure
col=get(gca,'colororder');col(3:4,:,:)=[];
clear h
set(gca,'linewidth',1.5,'XColor',[0.3 0.3 0.3],'YColor',[0.3 0.3 0.3],'Box','on','TickLength',[0.005 0.005],'FontSize',16)
hold on

m= [Datas.mean_glu; Datas.mean_gad]/1000;
se= [Datas.std_glu; Datas.std_gad]/1000;

errorbar(1:4,m,se,'.','color',[0.3 0.3 0.3],'LineWidth',2);
co = [1 2 1 2];
for i=1:4
hold on
bar(i,m(i),'Facecolor',col(co(i),:,:),'EdgeColor',[0.3 0.3 0.3],'LineWidth',2,'Barwidth',0.8);
end

% significativ=planes(Gad_stat>0.01 & Gad_stat<0.05);
% for ast=1:length(significativ)
%     ht=text(significativ(ast),max(Datas.mean_gad(Datas.treatment=="wt" & Datas.plane==significativ(ast))+Datas.es_gad(Datas.treatment=="wt" & Datas.plane==significativ(ast)),Datas.mean_gad(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast))+Datas.es_gad(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast)))+1,'*','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
% end
% significativ=planes(Gad_stat<=0.01);
% for ast=1:length(significativ)
%     ht=text(significativ(ast),max(Datas.mean_gad(Datas.treatment=="wt" & Datas.plane==significativ(ast))+Datas.es_gad(Datas.treatment=="wt" & Datas.plane==significativ(ast)),Datas.mean_gad(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast))+Datas.es_gad(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast)))+1,'**','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
% end


axis square
xlim([0.5 4.5]);
ylim([0 7]);
% ht=text(1.5,0.6-0.6*0.05,['p = ' sprintf('%.2E',stat)],'VerticalAlignment','middle','HorizontalAlignment','center');
box on
set(gca,'xtick',[1.5,3.3],'xticklabel',[{'Vglut2a'},{'gad1b'}]) %,[{'Vglut2a wild-type'},{'Vglut2a mecp2'},{'gad1b wild-type'},{'gad1b mecp2'}]
ylabel('Cell Number (x10^3)','FontSize',16);
axis square


















%%  --- Plane Distribution ---  %%

treatments=unique(ALL.Table.treatment);
animals=unique(ALL.Table.animal);
planes=unique(ALL.Table.plane);
planes=planes(1:3:end);
planes=planes(3:end-2);
temp=1;
clear Traeatmenss Cells_Numel Ratio_stat Glu_stat Gad_stat
for trea=1:length(treatments)
    for ani=1:length(animals) 
        for pla=1:(length(planes)-1)
            ArEa= 1;%sum(AREA.Table.area(AREA.Table.animal==animals(ani) & AREA.Table.plane>=planes(pla) & AREA.Table.plane<planes(pla+1) & AREA.Table.treatment==string(treatments(trea))));
            Cells_Numel(temp,:) = [str2num(animals(ani)), (planes(pla)+planes(pla+1))/2, numel(ALL.Table.glu(ALL.Table.animal==animals(ani) & ALL.Table.plane>=planes(pla) & ALL.Table.plane<planes(pla+1) & ALL.Table.treatment==string(treatments(trea)) & ALL.Table.color==1))/ArEa, numel(ALL.Table.gad(ALL.Table.animal==animals(ani) & ALL.Table.plane>=planes(pla) & ALL.Table.plane<planes(pla+1) & ALL.Table.treatment==string(treatments(trea)) &ALL.Table.color==2))/ArEa, numel(ALL.Table.dapi(ALL.Table.animal==animals(ani) & ALL.Table.plane>=planes(pla) & ALL.Table.plane<planes(pla+1) & ALL.Table.treatment==string(treatments(trea)) &ALL.Table.color==3))/ArEa];
            Traeatmenss(temp,1)=treatments(trea);
            temp=temp+1;
        end
    end
end

Cells_Numel_Table=table(Traeatmenss, Cells_Numel(:,1),Cells_Numel(:,2),Cells_Numel(:,3),Cells_Numel(:,4),Cells_Numel(:,5), Cells_Numel(:,3)./Cells_Numel(:,4)); %(Cells_Numel(:,4)-Cells_Numel(:,3))./(Cells_Numel(:,4)+Cells_Numel(:,3))
Cells_Numel_Table.Properties.VariableNames={'treatment','animal','plane','glu', 'gad', 'dapi', 'ratio'};
Cells_Numel_Table(isnan(Cells_Numel_Table.ratio),:)=[];

planes=unique(Cells_Numel_Table.plane);
for pla=1:length(planes)
    Ratio_stat(pla)=ranksum(Cells_Numel_Table.ratio(Cells_Numel_Table.treatment=="MeCP2" & Cells_Numel_Table.plane==planes(pla)),Cells_Numel_Table.ratio(Cells_Numel_Table.treatment=="wt" & Cells_Numel_Table.plane==planes(pla)));
    Glu_stat(pla)=ranksum(Cells_Numel_Table.glu(Cells_Numel_Table.treatment=="MeCP2" & Cells_Numel_Table.plane==planes(pla)),Cells_Numel_Table.glu(Cells_Numel_Table.treatment=="wt" & Cells_Numel_Table.plane==planes(pla)));
    Gad_stat(pla)=ranksum(Cells_Numel_Table.gad(Cells_Numel_Table.treatment=="MeCP2" & Cells_Numel_Table.plane==planes(pla)),Cells_Numel_Table.gad(Cells_Numel_Table.treatment=="wt" & Cells_Numel_Table.plane==planes(pla)));
end

Mean=varfun(@mean, Cells_Numel_Table,'InputVariables',{'glu', 'gad','ratio'},'GroupingVariables',{'treatment','plane'});
STD=varfun(@std, Cells_Numel_Table,'InputVariables',{'glu', 'gad','ratio'},'GroupingVariables',{'treatment','plane'});
ES=table(STD.std_glu./sqrt(STD.GroupCount), STD.std_gad./sqrt(STD.GroupCount), STD.std_ratio./sqrt(STD.GroupCount));
ES.Properties.VariableNames={'es_glu', 'es_gad', 'es_ratio'};
Datas=join(Mean,STD);
Datas=[Datas,ES];
Datas.plane=Datas.plane*2;
%stats
Cells_Numel_Table.animal(Cells_Numel_Table.treatment=="wt")=Cells_Numel_Table.animal(Cells_Numel_Table.treatment=="wt")+8;
glme = fitglme(Cells_Numel_Table,...
		'glu ~ plane + treatment + (1|animal)',...
		'Distribution','Poisson','Link','log','FitMethod','Laplace')

glme = fitglme(Cells_Numel_Table,...
		'gad ~ plane + treatment + (1|animal)',...
		'Distribution','Poisson','Link','log','FitMethod','Laplace')


col=get(gca,'colororder');col(3:4,:,:)=[];
figure;
set(gca,'linewidth',1.5,'XColor',[0.3 0.3 0.3],'YColor',[0.3 0.3 0.3],'Box','on','TickLength',[0.005 0.005],'FontSize',16)
hold on
h(2)=plot(Datas.plane(Datas.treatment=="MeCP2"),Datas.mean_ratio(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',2);
X=[Datas.plane(Datas.treatment=="MeCP2");flipud(Datas.plane(Datas.treatment=="MeCP2"))];
Y=[Datas.mean_ratio(Datas.treatment=="MeCP2")-Datas.es_ratio(Datas.treatment=="MeCP2");flipud(Datas.mean_ratio(Datas.treatment=="MeCP2")+Datas.es_ratio(Datas.treatment=="MeCP2"))];
h1=fill(X, Y, col(2,:,:),'EdgeColor','none');
set(h1,'FaceAlpha',0.3);
% errorbar(Datas.plane(Datas.treatment=="MeCP2"),Datas.mean_ratio(Datas.treatment=="MeCP2"), Datas.es_ratio(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',1);
h(1)=plot(Datas.plane(Datas.treatment=="wt"),Datas.mean_ratio(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',2);
X=[Datas.plane(Datas.treatment=="wt");flipud(Datas.plane(Datas.treatment=="wt"))];
Y=[Datas.mean_ratio(Datas.treatment=="wt")-Datas.es_ratio(Datas.treatment=="wt");flipud(Datas.mean_ratio(Datas.treatment=="wt")+Datas.es_ratio(Datas.treatment=="wt"))];
h2=fill(X, Y, col(1,:,:),'EdgeColor','none');
set(h2,'FaceAlpha',0.3);
% errorbar(Datas.plane(Datas.treatment=="wt"),Datas.mean_ratio(Datas.treatment=="wt"), Datas.es_ratio(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',1);
planes=unique(Datas.plane);
significativ=planes(Ratio_stat>0.01 & Ratio_stat<0.05);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_ratio(Datas.treatment=="wt" & Datas.plane==significativ(ast))+Datas.es_ratio(Datas.treatment=="wt" & Datas.plane==significativ(ast)),Datas.mean_ratio(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast))+Datas.es_ratio(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast)))+0.05,'*','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
significativ=planes(Ratio_stat<=0.01);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_ratio(Datas.treatment=="wt" & Datas.plane==significativ(ast))+Datas.es_ratio(Datas.treatment=="wt" & Datas.plane==significativ(ast)),Datas.mean_ratio(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast))+Datas.es_ratio(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast)))+0.05,'**','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
ylim([-0.5,2.5]);%([min(Datas.mean_ratio-Datas.es_ratio)-0.1 max(Datas.mean_ratio+Datas.es_ratio)+0.1]);
xlim([12,156]);%([min(Datas.plane) max(Datas.plane)]);
xlabel('Plane [um]','FontSize',16);
ylabel('Rario (Vglut2a/gad1b)','FontSize',16);
title(sprintf('%s', 'Ratio'));
l=legend(h,'wild-type','mecp2','Location','southeast');
l.Box='off';


figure
set(gca,'linewidth',1.5,'XColor',[0.3 0.3 0.3],'YColor',[0.3 0.3 0.3],'Box','on','TickLength',[0.005 0.005],'FontSize',16)
hold on
h(2)=plot(Datas.plane(Datas.treatment=="MeCP2"),Datas.mean_gad(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',2);
X=[Datas.plane(Datas.treatment=="MeCP2");flipud(Datas.plane(Datas.treatment=="MeCP2"))];
Y=[Datas.mean_gad(Datas.treatment=="MeCP2")-Datas.es_gad(Datas.treatment=="MeCP2");flipud(Datas.mean_gad(Datas.treatment=="MeCP2")+Datas.es_gad(Datas.treatment=="MeCP2"))];
h1=fill(X, Y, col(2,:,:),'EdgeColor','none');
set(h1,'FaceAlpha',0.3);
% errorbar(Datas.plane(Datas.treatment=="MeCP2"),Datas.mean_gad(Datas.treatment=="MeCP2"), Datas.es_gad(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',1);
h(1)=plot(Datas.plane(Datas.treatment=="wt"),Datas.mean_gad(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',2);
X=[Datas.plane(Datas.treatment=="wt");flipud(Datas.plane(Datas.treatment=="wt"))];
Y=[Datas.mean_gad(Datas.treatment=="wt")-Datas.es_gad(Datas.treatment=="wt");flipud(Datas.mean_gad(Datas.treatment=="wt")+Datas.es_gad(Datas.treatment=="wt"))];
h2=fill(X, Y, col(1,:,:),'EdgeColor','none');
set(h2,'FaceAlpha',0.3);
% errorbar(Datas.plane(Datas.treatment=="wt"),Datas.mean_gad(Datas.treatment=="wt"), Datas.es_gad(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',1);
significativ=planes(Gad_stat>0.01 & Gad_stat<0.05);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_gad(Datas.treatment=="wt" & Datas.plane==significativ(ast))+Datas.es_gad(Datas.treatment=="wt" & Datas.plane==significativ(ast)),Datas.mean_gad(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast))+Datas.es_gad(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast)))+1,'*','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
significativ=planes(Gad_stat<=0.01);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_gad(Datas.treatment=="wt" & Datas.plane==significativ(ast))+Datas.es_gad(Datas.treatment=="wt" & Datas.plane==significativ(ast)),Datas.mean_gad(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast))+Datas.es_gad(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast)))+1,'**','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
xlim([12,156]);%([min(Datas.plane) max(Datas.plane)]);
 ylim([0,450]);%([min(Datas.mean_gad-Datas.es_gad)-2 max(Datas.mean_gad+Datas.es_gad)+2]);
xlabel('Plane [um]','FontSize',16);
ylabel('Gad1b','FontSize',16);
title(sprintf('%s', 'Gad1b'));
l=legend(h,'wild-type','mecp2','Location','northeast');
l.Box='off';

figure
set(gca,'linewidth',1.5,'XColor',[0.3 0.3 0.3],'YColor',[0.3 0.3 0.3],'Box','on','TickLength',[0.005 0.005],'FontSize',16)
hold on
h(2)=plot(Datas.plane(Datas.treatment=="MeCP2"),Datas.mean_glu(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',2);
X=[Datas.plane(Datas.treatment=="MeCP2");flipud(Datas.plane(Datas.treatment=="MeCP2"))];
Y=[Datas.mean_glu(Datas.treatment=="MeCP2")-Datas.es_glu(Datas.treatment=="MeCP2");flipud(Datas.mean_glu(Datas.treatment=="MeCP2")+Datas.es_glu(Datas.treatment=="MeCP2"))];
h1=fill(X, Y, col(2,:,:),'EdgeColor','none');
set(h1,'FaceAlpha',0.3);
% errorbar(Datas.plane(Datas.treatment=="MeCP2"),Datas.mean_glu(Datas.treatment=="MeCP2"), Datas.es_glu(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',1);
h(1)=plot(Datas.plane(Datas.treatment=="wt"),Datas.mean_glu(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',2);
X=[Datas.plane(Datas.treatment=="wt");flipud(Datas.plane(Datas.treatment=="wt"))];
Y=[Datas.mean_glu(Datas.treatment=="wt")-Datas.es_glu(Datas.treatment=="wt");flipud(Datas.mean_glu(Datas.treatment=="wt")+Datas.es_glu(Datas.treatment=="wt"))];
h2=fill(X, Y, col(1,:,:),'EdgeColor','none');
set(h2,'FaceAlpha',0.3);
% errorbar(Datas.plane(Datas.treatment=="wt"),Datas.mean_glu(Datas.treatment=="wt"), Datas.es_glu(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',1);
significativ=planes(Glu_stat>0.01 & Glu_stat<0.05);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_glu(Datas.treatment=="wt" & Datas.plane==significativ(ast))+Datas.es_glu(Datas.treatment=="wt" & Datas.plane==significativ(ast)),Datas.mean_glu(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast))+Datas.es_glu(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast)))+1,'*','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
significativ=planes(Glu_stat<=0.01);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_glu(Datas.treatment=="wt" & Datas.plane==significativ(ast))+Datas.es_glu(Datas.treatment=="wt" & Datas.plane==significativ(ast)),Datas.mean_glu(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast))+Datas.es_glu(Datas.treatment=="MeCP2" & Datas.plane==significativ(ast)))+1,'**','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
xlim([12,156]);%([min(Datas.plane) max(Datas.plane)]);
 ylim([0,450]);%([min(Datas.mean_glu-Datas.es_glu)-2 max(Datas.mean_glu+Datas.es_glu)+2]);
xlabel('Plane [um]','FontSize',16);
ylabel('vGlut2a','FontSize',16);
title(sprintf('%s', 'vGlut2a'));
l=legend(h,'wild-type','mecp2','Location','northeast');
l.Box='off';





%%  --- Mean Line Distribution ---  %%

Line_segments=[0:51]/50; %[0:101 105:206]/100;
temp=1;
clear Traeatmenss Cells_Numel Ratio_stat Glu_stat Gad_stat
for trea=1:length(treatments)
    for ani=1:length(animals) 
        for lin=1:(length(Line_segments)-1)
            Cells_Numel(temp,:) = [str2num(animals(ani)), (Line_segments(lin)+Line_segments(lin+1))/2, numel(ALL.Table.glu(ALL.Table.animal==animals(ani) & ALL.Table.tectalLine>=Line_segments(lin) & ALL.Table.tectalLine<Line_segments(lin+1) & ALL.Table.treatment==string(treatments(trea)) & ALL.Table.color==1)), numel(ALL.Table.gad(ALL.Table.animal==animals(ani) & ALL.Table.tectalLine>=Line_segments(lin) & ALL.Table.tectalLine<Line_segments(lin+1) & ALL.Table.treatment==string(treatments(trea)) &ALL.Table.color==2)), numel(ALL.Table.dapi(ALL.Table.animal==animals(ani) & ALL.Table.tectalLine>=Line_segments(lin) & ALL.Table.tectalLine<Line_segments(lin+1) & ALL.Table.treatment==string(treatments(trea)) &ALL.Table.color==3))];
            Traeatmenss(temp,1)=treatments(trea);
            temp=temp+1;
        end
    end
end

Cells_Numel_Table=table(Traeatmenss, Cells_Numel(:,1),Cells_Numel(:,2),Cells_Numel(:,3),Cells_Numel(:,4),Cells_Numel(:,5), Cells_Numel(:,3)./Cells_Numel(:,4)); %(Cells_Numel(:,4)-Cells_Numel(:,3))./(Cells_Numel(:,4)+Cells_Numel(:,3))
Cells_Numel_Table.Properties.VariableNames={'treatment','animal','mean_line','glu', 'gad', 'dapi', 'ratio'};
Cells_Numel_Table(isnan(Cells_Numel_Table.ratio),:)=[];

mean_lines=unique(Cells_Numel_Table.mean_line);
for lin=1:length(mean_lines)
    Ratio_stat(lin)=ranksum(Cells_Numel_Table.ratio(Cells_Numel_Table.treatment=="MeCP2" & Cells_Numel_Table.mean_line==mean_lines(lin)),Cells_Numel_Table.ratio(Cells_Numel_Table.treatment=="wt" & Cells_Numel_Table.mean_line==mean_lines(lin)));
    Glu_stat(lin)=ranksum(Cells_Numel_Table.glu(Cells_Numel_Table.treatment=="MeCP2" & Cells_Numel_Table.mean_line==mean_lines(lin)),Cells_Numel_Table.glu(Cells_Numel_Table.treatment=="wt" & Cells_Numel_Table.mean_line==mean_lines(lin)));
    Gad_stat(lin)=ranksum(Cells_Numel_Table.gad(Cells_Numel_Table.treatment=="MeCP2" & Cells_Numel_Table.mean_line==mean_lines(lin)),Cells_Numel_Table.gad(Cells_Numel_Table.treatment=="wt" & Cells_Numel_Table.mean_line==mean_lines(lin)));
end

Mean=varfun(@mean, Cells_Numel_Table,'InputVariables',{'glu', 'gad','ratio'},'GroupingVariables',{'treatment','mean_line'});
STD=varfun(@std, Cells_Numel_Table,'InputVariables',{'glu', 'gad','ratio'},'GroupingVariables',{'treatment','mean_line'});
ES=table(STD.std_glu./sqrt(STD.GroupCount), STD.std_gad./sqrt(STD.GroupCount), STD.std_ratio./sqrt(STD.GroupCount));
ES.Properties.VariableNames={'es_glu', 'es_gad', 'es_ratio'};
Datas=join(Mean,STD);
Datas=[Datas,ES];

%stats
Cells_Numel_Table.animal(Cells_Numel_Table.treatment=="wt")=Cells_Numel_Table.animal(Cells_Numel_Table.treatment=="wt")+8;
glme = fitglme(Cells_Numel_Table,...
		'glu ~ mean_line + treatment + (1|animal)',...
		'Distribution','Poisson','Link','log','FitMethod','Laplace')

glme = fitglme(Cells_Numel_Table,...
		'gad ~ mean_line + treatment + (1|animal)',...
		'Distribution','Poisson','Link','log','FitMethod','Laplace')


figure
set(gca,'linewidth',1.5,'XColor',[0.3 0.3 0.3],'YColor',[0.3 0.3 0.3],'Box','on','TickLength',[0.005 0.005],'FontSize',16)
hold on
h(2)=plot(Datas.mean_line(Datas.treatment=="MeCP2"),Datas.mean_ratio(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',2);
X=[Datas.mean_line(Datas.treatment=="MeCP2");flipud(Datas.mean_line(Datas.treatment=="MeCP2"))];
Y=[Datas.mean_ratio(Datas.treatment=="MeCP2")-Datas.es_ratio(Datas.treatment=="MeCP2");flipud(Datas.mean_ratio(Datas.treatment=="MeCP2")+Datas.es_ratio(Datas.treatment=="MeCP2"))];
h1=fill(X, Y, col(2,:,:),'EdgeColor','none');
set(h1,'FaceAlpha',0.3);
% errorbar(Datas.mean_line(Datas.treatment=="MeCP2"),Datas.mean_ratio(Datas.treatment=="MeCP2"), Datas.es_ratio(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',1);
h(1)=plot(Datas.mean_line(Datas.treatment=="wt"),Datas.mean_ratio(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',2);
X=[Datas.mean_line(Datas.treatment=="wt");flipud(Datas.mean_line(Datas.treatment=="wt"))];
Y=[Datas.mean_ratio(Datas.treatment=="wt")-Datas.es_ratio(Datas.treatment=="wt");flipud(Datas.mean_ratio(Datas.treatment=="wt")+Datas.es_ratio(Datas.treatment=="wt"))];
h2=fill(X, Y, col(1,:,:),'EdgeColor','none');
set(h2,'FaceAlpha',0.3);
% errorbar(Datas.mean_line(Datas.treatment=="wt"),Datas.mean_ratio(Datas.treatment=="wt"), Datas.es_ratio(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',1);
significativ=mean_lines(Ratio_stat>0.01 & Ratio_stat<0.05);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_ratio(Datas.treatment=="wt" & Datas.mean_line==significativ(ast))+Datas.es_ratio(Datas.treatment=="wt" & Datas.mean_line==significativ(ast)),Datas.mean_ratio(Datas.treatment=="MeCP2" & Datas.mean_line==significativ(ast))+Datas.es_ratio(Datas.treatment=="MeCP2" & Datas.mean_line==significativ(ast)))+0.05,'*','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
significativ=mean_lines(Ratio_stat<=0.01);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_ratio(Datas.treatment=="wt" & Datas.mean_line==significativ(ast))+Datas.es_ratio(Datas.treatment=="wt" & Datas.mean_line==significativ(ast)),Datas.mean_ratio(Datas.treatment=="MeCP2" & Datas.mean_line==significativ(ast))+Datas.es_ratio(Datas.treatment=="MeCP2" & Datas.mean_line==significativ(ast)))+0.05,'**','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
ylim([-0.5,2.5]);%([min(Datas.mean_ratio-Datas.es_ratio)-0.1 max(Datas.mean_ratio+Datas.es_ratio)+0.1]);
xlim([-0.01,1.03]);%([min(Datas.mean_line) max(Datas.mean_line)]);
xlabel('Tectal Mean Line','FontSize',16);
ylabel('Rario (Vglut2a/gad1b)','FontSize',16);
title(sprintf('%s', 'Ratio'));
l=legend(h,'wild-type','mecp2','Location','southeast');
l.Box='off';



figure
set(gca,'linewidth',1.5,'XColor',[0.3 0.3 0.3],'YColor',[0.3 0.3 0.3],'Box','on','TickLength',[0.005 0.005],'FontSize',16)
hold on
h(2)=plot(Datas.mean_line(Datas.treatment=="MeCP2"),Datas.mean_gad(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',2);
X=[Datas.mean_line(Datas.treatment=="MeCP2");flipud(Datas.mean_line(Datas.treatment=="MeCP2"))];
Y=[Datas.mean_gad(Datas.treatment=="MeCP2")-Datas.es_gad(Datas.treatment=="MeCP2");flipud(Datas.mean_gad(Datas.treatment=="MeCP2")+Datas.es_gad(Datas.treatment=="MeCP2"))];
h1=fill(X, Y, col(2,:,:),'EdgeColor','none');
set(h1,'FaceAlpha',0.3);
% errorbar(Datas.mean_line(Datas.treatment=="MeCP2"),Datas.mean_gad(Datas.treatment=="MeCP2"), Datas.es_gad(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',1);
h(1)=plot(Datas.mean_line(Datas.treatment=="wt"),Datas.mean_gad(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',2);
X=[Datas.mean_line(Datas.treatment=="wt");flipud(Datas.mean_line(Datas.treatment=="wt"))];
Y=[Datas.mean_gad(Datas.treatment=="wt")-Datas.es_gad(Datas.treatment=="wt");flipud(Datas.mean_gad(Datas.treatment=="wt")+Datas.es_gad(Datas.treatment=="wt"))];
h2=fill(X, Y, col(1,:,:),'EdgeColor','none');
set(h2,'FaceAlpha',0.3);
% errorbar(Datas.mean_line(Datas.treatment=="wt"),Datas.mean_gad(Datas.treatment=="wt"), Datas.es_gad(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',1);
significativ=mean_lines(Gad_stat>0.01 & Gad_stat<0.05);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_gad(Datas.treatment=="wt" & Datas.mean_line==significativ(ast))+Datas.es_gad(Datas.treatment=="wt" & Datas.mean_line==significativ(ast)),Datas.mean_gad(Datas.treatment=="MeCP2" & Datas.mean_line==significativ(ast))+Datas.es_gad(Datas.treatment=="MeCP2" & Datas.mean_line==significativ(ast)))+1,'*','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
significativ=mean_lines(Gad_stat<=0.01);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_gad(Datas.treatment=="wt" & Datas.mean_line==significativ(ast))+Datas.es_gad(Datas.treatment=="wt" & Datas.mean_line==significativ(ast)),Datas.mean_gad(Datas.treatment=="MeCP2" & Datas.mean_line==significativ(ast))+Datas.es_gad(Datas.treatment=="MeCP2" & Datas.mean_line==significativ(ast)))+1,'**','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
xlim([-0.01,1.03]);%([min(Datas.mean_line) max(Datas.mean_line)]);
ylim([0,180]);%([min(Datas.mean_gad-Datas.es_gad)-2 max(Datas.mean_gad+Datas.es_gad)+2]);
xlabel('Tectal Mean Line','FontSize',16);
ylabel('Gad1b','FontSize',16);
title(sprintf('%s', 'Gad1b'));
l=legend(h,'wild-type','mecp2','Location','southeast');
l.Box='off';

figure
set(gca,'linewidth',1.5,'XColor',[0.3 0.3 0.3],'YColor',[0.3 0.3 0.3],'Box','on','TickLength',[0.005 0.005],'FontSize',16)
hold on
h(2)=plot(Datas.mean_line(Datas.treatment=="MeCP2"),Datas.mean_glu(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',2);
X=[Datas.mean_line(Datas.treatment=="MeCP2");flipud(Datas.mean_line(Datas.treatment=="MeCP2"))];
Y=[Datas.mean_glu(Datas.treatment=="MeCP2")-Datas.es_glu(Datas.treatment=="MeCP2");flipud(Datas.mean_glu(Datas.treatment=="MeCP2")+Datas.es_glu(Datas.treatment=="MeCP2"))];
h1=fill(X, Y, col(2,:,:),'EdgeColor','none');
set(h1,'FaceAlpha',0.3);
% errorbar(Datas.mean_line(Datas.treatment=="MeCP2"),Datas.mean_glu(Datas.treatment=="MeCP2"), Datas.es_glu(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',1);
h(1)=plot(Datas.mean_line(Datas.treatment=="wt"),Datas.mean_glu(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',2);
X=[Datas.mean_line(Datas.treatment=="wt");flipud(Datas.mean_line(Datas.treatment=="wt"))];
Y=[Datas.mean_glu(Datas.treatment=="wt")-Datas.es_glu(Datas.treatment=="wt");flipud(Datas.mean_glu(Datas.treatment=="wt")+Datas.es_glu(Datas.treatment=="wt"))];
h2=fill(X, Y, col(1,:,:),'EdgeColor','none');
set(h2,'FaceAlpha',0.3);
% errorbar(Datas.mean_line(Datas.treatment=="wt"),Datas.mean_glu(Datas.treatment=="wt"), Datas.es_glu(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',1);
significativ=mean_lines(Glu_stat>0.01 & Glu_stat<0.05);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_glu(Datas.treatment=="wt" & Datas.mean_line==significativ(ast))+Datas.es_glu(Datas.treatment=="wt" & Datas.mean_line==significativ(ast)),Datas.mean_glu(Datas.treatment=="MeCP2" & Datas.mean_line==significativ(ast))+Datas.es_glu(Datas.treatment=="MeCP2" & Datas.mean_line==significativ(ast)))+1,'*','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
significativ=mean_lines(Glu_stat<=0.01);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_glu(Datas.treatment=="wt" & Datas.mean_line==significativ(ast))+Datas.es_glu(Datas.treatment=="wt" & Datas.mean_line==significativ(ast)),Datas.mean_glu(Datas.treatment=="MeCP2" & Datas.mean_line==significativ(ast))+Datas.es_glu(Datas.treatment=="MeCP2" & Datas.mean_line==significativ(ast)))+1,'**','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
xlim([-0.01,1.03]);%([min(Datas.mean_line) max(Datas.mean_line)]);
ylim([0,180]);%([min(Datas.mean_glu-Datas.es_glu)-2 max(Datas.mean_glu+Datas.es_glu)+2]);
xlabel('Tectal Mean Line','FontSize',16);
ylabel('vGlut2a','FontSize',16);
title(sprintf('%s', 'vGlut2a'));
l=legend(h,'wild-type','mecp2','Location','southeast');
l.Box='off';




%%  --- Ortogonal to Mean Line Distribution ---  %%

% Distances=round(min(ALL.Table.dist))-1:length(round(min(ALL.Table.dist))-1:round(max(ALL.Table.dist)))/50:round(max(ALL.Table.dist))+1;
Distances=[-110:6:110];
temp=1;
clear Traeatmenss Cells_Numel Ratio_stat Glu_stat Gad_stat
for trea=1:length(treatments)
    for ani=1:length(animals) 
        for dis=1:(length(Distances)-1)
            Cells_Numel(temp,:) = [str2num(animals(ani)), (Distances(dis)+Distances(dis+1))/2, numel(ALL.Table.glu(ALL.Table.animal==animals(ani) & ALL.Table.dist>=Distances(dis) & ALL.Table.dist<Distances(dis+1) & ALL.Table.treatment==string(treatments(trea)) & ALL.Table.color==1)), numel(ALL.Table.gad(ALL.Table.animal==animals(ani) & ALL.Table.dist>=Distances(dis) & ALL.Table.dist<Distances(dis+1) & ALL.Table.treatment==string(treatments(trea)) &ALL.Table.color==2)), numel(ALL.Table.dapi(ALL.Table.animal==animals(ani) & ALL.Table.dist>=Distances(dis) & ALL.Table.dist<Distances(dis+1) & ALL.Table.treatment==string(treatments(trea)) &ALL.Table.color==3))];
            Traeatmenss(temp,1)=treatments(trea);
            temp=temp+1;
        end
    end
end

Cells_Numel_Table=table(Traeatmenss, Cells_Numel(:,1),Cells_Numel(:,2),Cells_Numel(:,3),Cells_Numel(:,4),Cells_Numel(:,5), Cells_Numel(:,3)./Cells_Numel(:,4)); %(Cells_Numel(:,4)-Cells_Numel(:,3))./(Cells_Numel(:,4)+Cells_Numel(:,3))
Cells_Numel_Table.Properties.VariableNames={'treatment','animal','distance','glu', 'gad', 'dapi', 'ratio'};
Cells_Numel_Table(isnan(Cells_Numel_Table.ratio),:)=[];

Distt=unique(Cells_Numel_Table.distance);
for dis=1:length(Distt)
    if isempty(Cells_Numel_Table.ratio(Cells_Numel_Table.treatment=="wt" & Cells_Numel_Table.distance==Distt(dis)))==1
        Ratio_stat(dis)=1;
        Glu_stat(dis)=1;
        Gad_stat(dis)=1;
    else
        Ratio_stat(dis)=ranksum(Cells_Numel_Table.ratio(Cells_Numel_Table.treatment=="MeCP2" & Cells_Numel_Table.distance==Distt(dis)),Cells_Numel_Table.ratio(Cells_Numel_Table.treatment=="wt" & Cells_Numel_Table.distance==Distt(dis)));
        Glu_stat(dis)=ranksum(Cells_Numel_Table.glu(Cells_Numel_Table.treatment=="MeCP2" & Cells_Numel_Table.distance==Distt(dis)),Cells_Numel_Table.glu(Cells_Numel_Table.treatment=="wt" & Cells_Numel_Table.distance==Distt(dis)));
        Gad_stat(dis)=ranksum(Cells_Numel_Table.gad(Cells_Numel_Table.treatment=="MeCP2" & Cells_Numel_Table.distance==Distt(dis)),Cells_Numel_Table.gad(Cells_Numel_Table.treatment=="wt" & Cells_Numel_Table.distance==Distt(dis)));
    end
end

Mean=varfun(@mean, Cells_Numel_Table,'InputVariables',{'glu', 'gad','ratio'},'GroupingVariables',{'treatment','distance'});
STD=varfun(@std, Cells_Numel_Table,'InputVariables',{'glu', 'gad','ratio'},'GroupingVariables',{'treatment','distance'});
ES=table(STD.std_glu./sqrt(STD.GroupCount), STD.std_gad./sqrt(STD.GroupCount), STD.std_ratio./sqrt(STD.GroupCount));
ES.Properties.VariableNames={'es_glu', 'es_gad', 'es_ratio'};
Datas=join(Mean,STD);
Datas=[Datas,ES];


%stats
Cells_Numel_Table.animal(Cells_Numel_Table.treatment=="wt")=Cells_Numel_Table.animal(Cells_Numel_Table.treatment=="wt")+8;
glme = fitglme(Cells_Numel_Table,...
		'glu ~ distance + treatment + (1|animal)',...
		'Distribution','Poisson','Link','log','FitMethod','Laplace')

glme = fitglme(Cells_Numel_Table,...
		'gad ~ distance + treatment + (1|animal)',...
		'Distribution','Poisson','Link','log','FitMethod','Laplace')




figure
set(gca,'linewidth',1.5,'XColor',[0.3 0.3 0.3],'YColor',[0.3 0.3 0.3],'Box','on','TickLength',[0.005 0.005],'FontSize',16)
hold on
h(2)=plot(Datas.distance(Datas.treatment=="MeCP2"),Datas.mean_ratio(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',2);
X=[Datas.distance(Datas.treatment=="MeCP2");flipud(Datas.distance(Datas.treatment=="MeCP2"))];
Y=[Datas.mean_ratio(Datas.treatment=="MeCP2")-Datas.es_ratio(Datas.treatment=="MeCP2");flipud(Datas.mean_ratio(Datas.treatment=="MeCP2")+Datas.es_ratio(Datas.treatment=="MeCP2"))];
h1=fill(X, Y, col(2,:,:),'EdgeColor','none');
set(h1,'FaceAlpha',0.3);
% errorbar(Datas.distance(Datas.treatment=="MeCP2"),Datas.mean_ratio(Datas.treatment=="MeCP2"), Datas.es_ratio(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',1);
h(1)=plot(Datas.distance(Datas.treatment=="wt"),Datas.mean_ratio(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',2);
X=[Datas.distance(Datas.treatment=="wt");flipud(Datas.distance(Datas.treatment=="wt"))];
Y=[Datas.mean_ratio(Datas.treatment=="wt")-Datas.es_ratio(Datas.treatment=="wt");flipud(Datas.mean_ratio(Datas.treatment=="wt")+Datas.es_ratio(Datas.treatment=="wt"))];
h2=fill(X, Y, col(1,:,:),'EdgeColor','none');
set(h2,'FaceAlpha',0.3);
% errorbar(Datas.distance(Datas.treatment=="wt"),Datas.mean_ratio(Datas.treatment=="wt"), Datas.es_ratio(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',1);
significativ=Distt(Ratio_stat>0.01 & Ratio_stat<0.05);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_ratio(Datas.treatment=="wt" & Datas.distance==significativ(ast))+Datas.es_ratio(Datas.treatment=="wt" & Datas.distance==significativ(ast)),Datas.mean_ratio(Datas.treatment=="MeCP2" & Datas.distance==significativ(ast))+Datas.es_ratio(Datas.treatment=="MeCP2" & Datas.distance==significativ(ast)))+0.05,'*','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
significativ=Distt(Ratio_stat<=0.01);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_ratio(Datas.treatment=="wt" & Datas.distance==significativ(ast))+Datas.es_ratio(Datas.treatment=="wt" & Datas.distance==significativ(ast)),Datas.mean_ratio(Datas.treatment=="MeCP2" & Datas.distance==significativ(ast))+Datas.es_ratio(Datas.treatment=="MeCP2" & Datas.distance==significativ(ast)))+0.05,'**','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
ylim([-0.5,2.5]);%([min(Datas.mean_ratio-Datas.es_ratio)-0.1 max(Datas.mean_ratio+Datas.es_ratio)+0.1]);
xlim([-111,107]);%([min(Datas.distance) max(Datas.distance)]);
xlabel('Distance to Tectal Mean Line','FontSize',16);
ylabel('Rario (Vglut2a/gad1b)','FontSize',16);
title(sprintf('%s', 'Ratio'));
l=legend(h,'wild-type','mecp2','Location','southeast');
l.Box='off';



figure
set(gca,'linewidth',1.5,'XColor',[0.3 0.3 0.3],'YColor',[0.3 0.3 0.3],'Box','on','TickLength',[0.005 0.005],'FontSize',16)
hold on
h(2)=plot(Datas.distance(Datas.treatment=="MeCP2"),Datas.mean_gad(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',2);
X=[Datas.distance(Datas.treatment=="MeCP2");flipud(Datas.distance(Datas.treatment=="MeCP2"))];
Y=[Datas.mean_gad(Datas.treatment=="MeCP2")-Datas.es_gad(Datas.treatment=="MeCP2");flipud(Datas.mean_gad(Datas.treatment=="MeCP2")+Datas.es_gad(Datas.treatment=="MeCP2"))];
h1=fill(X, Y, col(2,:,:),'EdgeColor','none');
set(h1,'FaceAlpha',0.3);
% errorbar(Datas.distance(Datas.treatment=="MeCP2"),Datas.mean_gad(Datas.treatment=="MeCP2"), Datas.es_gad(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',1);
h(1)=plot(Datas.distance(Datas.treatment=="wt"),Datas.mean_gad(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',2);
X=[Datas.distance(Datas.treatment=="wt");flipud(Datas.distance(Datas.treatment=="wt"))];
Y=[Datas.mean_gad(Datas.treatment=="wt")- Datas.es_gad(Datas.treatment=="wt");flipud(Datas.mean_gad(Datas.treatment=="wt")+ Datas.es_gad(Datas.treatment=="wt"))];
h2=fill(X, Y, col(1,:,:),'EdgeColor','none');
set(h2,'FaceAlpha',0.3);
% errorbar(Datas.distance(Datas.treatment=="wt"),Datas.mean_gad(Datas.treatment=="wt"), Datas.es_gad(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',1);
significativ=Distt(Gad_stat>0.01 & Gad_stat<0.05);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_gad(Datas.treatment=="wt" & Datas.distance==significativ(ast))+Datas.es_gad(Datas.treatment=="wt" & Datas.distance==significativ(ast)),Datas.mean_gad(Datas.treatment=="MeCP2" & Datas.distance==significativ(ast))+Datas.es_gad(Datas.treatment=="MeCP2" & Datas.distance==significativ(ast)))+1,'*','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
significativ=Distt(Gad_stat<=0.01);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_gad(Datas.treatment=="wt" & Datas.distance==significativ(ast))+Datas.es_gad(Datas.treatment=="wt" & Datas.distance==significativ(ast)),Datas.mean_gad(Datas.treatment=="MeCP2" & Datas.distance==significativ(ast))+Datas.es_gad(Datas.treatment=="MeCP2" & Datas.distance==significativ(ast)))+1,'**','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
xlim([-111,107]);%([min(Datas.distance) max(Datas.distance)]);
ylim([0,450]);%([min(Datas.mean_gad-Datas.es_gad)-2 max(Datas.mean_gad+Datas.es_gad)+2]);
xlabel('Distance to Tectal Mean Line','FontSize',16);
ylabel('Gad1b','FontSize',16);
title(sprintf('%s', 'Gad1b'));
l=legend(h,'wild-type','mecp2','Location','northeast');
l.Box='off';

figure
set(gca,'linewidth',1.5,'XColor',[0.3 0.3 0.3],'YColor',[0.3 0.3 0.3],'Box','on','TickLength',[0.005 0.005],'FontSize',16)
hold on
h(2)=plot(Datas.distance(Datas.treatment=="MeCP2"),Datas.mean_glu(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',2);
X=[Datas.distance(Datas.treatment=="MeCP2");flipud(Datas.distance(Datas.treatment=="MeCP2"))];
Y=[Datas.mean_glu(Datas.treatment=="MeCP2")-Datas.es_glu(Datas.treatment=="MeCP2");flipud(Datas.mean_glu(Datas.treatment=="MeCP2")+Datas.es_glu(Datas.treatment=="MeCP2"))];
h1=fill(X, Y, col(2,:,:),'EdgeColor','none');
set(h1,'FaceAlpha',0.3);
% errorbar(Datas.distance(Datas.treatment=="MeCP2"),Datas.mean_glu(Datas.treatment=="MeCP2"), Datas.es_glu(Datas.treatment=="MeCP2"), 'color',col(2,:,:),'LineWidth',1);
h(1)=plot(Datas.distance(Datas.treatment=="wt"),Datas.mean_glu(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',2);
X=[Datas.distance(Datas.treatment=="wt");flipud(Datas.distance(Datas.treatment=="wt"))];
Y=[Datas.mean_glu(Datas.treatment=="wt")-Datas.es_glu(Datas.treatment=="wt");flipud(Datas.mean_glu(Datas.treatment=="wt")+Datas.es_glu(Datas.treatment=="wt"))];
h2=fill(X, Y, col(1,:,:),'EdgeColor','none');
set(h2,'FaceAlpha',0.3);
% errorbar(Datas.distance(Datas.treatment=="wt"),Datas.mean_glu(Datas.treatment=="wt"), Datas.es_glu(Datas.treatment=="wt"), 'color',col(1,:,:),'LineWidth',1);
significativ=Distt(Glu_stat>0.01 & Glu_stat<0.05);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_glu(Datas.treatment=="wt" & Datas.distance==significativ(ast))+Datas.es_glu(Datas.treatment=="wt" & Datas.distance==significativ(ast)),Datas.mean_glu(Datas.treatment=="MeCP2" & Datas.distance==significativ(ast))+Datas.es_glu(Datas.treatment=="MeCP2" & Datas.distance==significativ(ast)))+1,'*','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
significativ=Distt(Glu_stat<=0.01);
for ast=1:length(significativ)
    ht=text(significativ(ast),max(Datas.mean_glu(Datas.treatment=="wt" & Datas.distance==significativ(ast))+Datas.es_glu(Datas.treatment=="wt" & Datas.distance==significativ(ast)),Datas.mean_glu(Datas.treatment=="MeCP2" & Datas.distance==significativ(ast))+Datas.es_glu(Datas.treatment=="MeCP2" & Datas.distance==significativ(ast)))+1,'**','VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
end
xlim([-111,107]);%([min(Datas.distance) max(Datas.distance)]);
ylim([0,450]);%([min(Datas.mean_glu-Datas.es_glu)-2 max(Datas.mean_glu+Datas.es_glu)+2]);
xlabel('Distance to Tectal Mean Line','FontSize',16);
ylabel('vGlut2a','FontSize',16);
title(sprintf('%s', 'vGlut2a'));
l=legend(h,'wild-type','mecp2','Location','northeast');
l.Box='off';






