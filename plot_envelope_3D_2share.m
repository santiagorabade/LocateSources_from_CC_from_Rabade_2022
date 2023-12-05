clear all
close all
load MyColormap
paIn=[pwd,'/matfiles/']; %% path to header and valueZ files
pathFigs=[pwd,'/figs_location/']; %% path to save figs
paEq=['/uufs/chpc.utah.edu/common/home/flin-group3/srabade/MINES/Ver_2/Coordinates/']; %% path to earthquake data
load([paEq,'datesEq_R2.mat'])
load([paEq,'locEq_R2.mat'])
depEq=load([paEq,'depthR2.txt']);
lonEq=lon;
latEq=lat;
clear lat lon
t1=((0:10000)*0.01)-50;
t2=((0:5000)*0.02)-50;
period='1-5';
plotMax=0;
slow=1;
cross=0;
load geometry.mat
sourceLa=39.663:0.0005:39.686;
sourceLo=-111.235:-0.0005:-111.258;
daylist=load('timelist0726');
for iD=1:length(daylist)   
    timeCounter=1;
    day=[daylist(iD,:)]
    allcoords=combvec(sourceLo,sourceLa); %combine lats ands lons
    sourceLoP=sourceLa-39.6752; %% to make it relative to one location
    sourceLaP=sourceLo+111.247; %% to make it relative to one location
    [A B]=meshgrid(sourceLoP,sourceLaP);
    zdummie=ones(length(lat),1)*100000; %% to plot stations on top
    depth=0%[0,0.2,0.4,0.6,0.8,1.0,1.2,1.4]; %depth vector
    for idep=1:length(depth)
        depAux=depth(idep);
        load ([paIn,'valueZ_',period,'_slow_',num2str(slow),'_day_',num2str(day),'_dep_',num2str(depAux),'_newCoord_all_16.mat']) %%read Amp value
        load ([paIn,'header_',period,'_slow_',num2str(slow),'_day_',num2str(day),'_dep_',num2str(depAux),'_newCoord_all_16.mat']) %%reaz headers
        accept=1;
        for iEnv=1:length(header)
            if header(iEnv).rmsIn>header(iEnv).rmsOut*2 %% check "SNR" prob this needs to be updated
                newZ(:,accept,idep)=valueZ(:,iEnv);
                for iZ=1:size(newZ,1)
                    if newZ(iZ,accept,idep)<header(iEnv).rmsIn
                        newZ(iZ,accept,idep)=0;
                    end
                end
            else
                for ifail=1:size(valueZ,1)
                    newZ(ifail,accept,idep)=NaN;
                end
            end
            accept=accept+1;
        end
        
        maxbyColum(:,idep)=max(newZ(:,:,idep));
        acceptSave(idep)=accept;
    end
    maxbycolum2=max(maxbyColum,[],2); %% to normalize 
    
    %%%find earthquake in time window from catalogue
    yy='2018';
    dayAux=num2str(day);
    hourAux=str2num(dayAux(9:end));
    hour=['0', num2str(floor((hourAux/3600)))];
    min=['0', num2str(rem(hourAux,3600)/60)];
    min2=['0', num2str(rem(hourAux+300,3600)/60)];
    tlower=datetime([yy,'-',dayAux(5:6),'-',dayAux(7:8),' ',hour(end-1:end),':',min(end-1:end),':00.000'],...
        'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    tupper=datetime([yy,'-',dayAux(5:6),'-',dayAux(7:8),' ',hour(end-1:end),':',min2(end-1:end),':00.000'],...
        'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
    tf = isbetween(datetime(datstr,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'),tlower,tupper);
    tfIn=find(tf);
    zeq=ones(1,size(tfIn,1))*100;
    %%%
    
    
    
    for i=1: length(depth)
        suma(:,:,i)=nansum(newZ(:,:,i)./maxbycolum2',2); %% stack
        suma3D(:,:,i)=reshape(suma(:,:,i)/max(max(suma)),length(sourceLo),length(sourceLa)); %%3D stack
        zbis=reshape(suma(:,:,i)/max(max(suma)),length(sourceLo),length(sourceLa));
        
        h=surfc(B,A,zbis,'LineStyle','none');
        shading interp
        colormap parula
        set(h(2),'LineColor','k');
        view(2);
        ylim([-0.006 0.006])
        xlim([-0.006 0.006])
        xticks([-0.0055 -0.00275  0 0.00275  0.0055 ])
        xticklabels({'-500','-250','0','250','500'})
        yticks([-0.0055 -0.00275  0 0.00275  0.0055 ])
        yticklabels({'-500','-250','0','250','500'})
        xlabel(['Distance (m)'])
        ylabel(['Distance (m)'])
        set(gca,'Children',flip(h),'sortmethod','childorder')
        c=colorbar;
        c.Label.String = 'Normalized amplitude';
        title([num2str(day),' CC ',num2str(accept),' Eq ', num2str(size(tfIn,1)),' Dp ', num2str(depth(i)) ])
        caxis([0.4 1])
        hold on
        plot3(lon(1:16)+111.247,lat(1:16)-39.6752,zdummie(1:16),'k^','MarkerFaceColor', 'w','MarkerSize',7,'LineWidth',1.5)
        hold on
        plot3(lonEq(tfIn)+111.247,latEq(tfIn)-39.6752,zeq,'wo','MarkerFaceColor', 'none','MarkerSize',3.5,'LineWidth',1.0)
        print([pathFigs,num2str(day),'_',num2str(depth(i)),'_2D.png'],'-dpng','-r0')
        set(gcf,'Renderer','Painter')
        hgexport(gcf,[pathFigs,num2str(day),'_',num2str(depth(i)),'_2D.eps']);
        close   
        
    end
   if cross==1 
    for iCross=[21,22,23,24,25,26,27,28,29,30,31] %index of desired cross section
        %%%%%%% Plot across Latitude cross-seccion
        figure
        tosurf3D= reshape(suma3D(:,iCross,:),size(suma3D,1),size(suma3D,3));
        [depth3d lon3d]=meshgrid(depth,sourceLoP);
        h=surfc(lon3d,depth3d,tosurf3D,'LineStyle','none');
        shading interp
        colormap parula
        set(h(2),'LineColor','k');
        set(gca,'Children',flip(h),'sortmethod','childorder')
        view(2)
        set(gca, 'YDir','reverse')
        caxis([0.4 1])
        c=colorbar;
        c.Label.String = 'Normalized amplitude';
        hold on
        plot3(latEq(tfIn)-39.6752,(2710-abs(depEq(tfIn)))/1000,zeq,'wo','MarkerFaceColor', 'none','MarkerSize',3.5,'LineWidth',1.5)
        xlim([-0.0055,0.0055])
        %xlim([39.669 39.682])
        ylim([0 1])
        xlabel(['Distance(m)'])
        xticks([-0.00832 -0.0055 -0.00275  0 0.00275  0.0055 0.00832])
        xticklabels({'-750','-500','-250','0','250','500','750'})
        ylabel(['Depth (m)'])
        yticks([0 0.250 0.5 0.750 1. 1.25])
        yticklabels({'0','250','500','750','1000','1250'})
        title([num2str(day),' y ', num2str((sourceLo(iCross)+111.247)*111),' Eq ', num2str(size(tfIn,1)) ])
        set(gca, 'Layer', 'top');
        grid off
        set(gcf,'Renderer','Painter')
        %hgexport(gcf,[pathFigs,num2str(day),'_y_',num2str(sourceLaP(iCross)*111),'_2D.eps']);
        close 
    end

    for iCross=[19,20,21,22,23,24,25,26] %index of desired cross section
        %%%%%% plot across Longitude cross-seccion
        figure
        tosurf3D= reshape(suma3D(iCross,:,:),size(suma3D,1),size(suma3D,3));
        [depth3d lat3d]=meshgrid(depth,sourceLaP);
        h=surfc(lat3d,depth3d,tosurf3D,'LineStyle','none');
        shading interp
        colormap parula
        set(h(2),'LineColor','k');
        set(gca,'Children',flip(h),'sortmethod','childorder')
        view(2)
        set(gca, 'YDir','reverse')
        caxis([0.4 1])
        c=colorbar;
        hold on
        plot3(lonEq(tfIn)+111.247,(2710-abs(depEq(tfIn)))/1000,zeq,'wo','MarkerFaceColor', 'none','MarkerSize',3.5,'LineWidth',1.5)
        c.Label.String = 'Normalized amplitude';
        xlim([-0.0055,0.0055])
        ylim([0 1])
        xlabel(['Distance(m)'])
        xticks([-0.00832 -0.0055 -0.00275  0 0.00275  0.0055 0.00832])
        xticklabels({'-750','-500','-250','0','250','500','750'})
        ylabel(['Depth (m)'])
        yticks([0 0.250 0.5 0.750 1. 1.25])
        yticklabels({'0','250','500','750','1000','1250'})
        title([num2str(day),' x ', num2str((sourceLa(iCross)-39.6752)*111),' Eq ', num2str(size(tfIn,1)) ])
         set(gca, 'Layer', 'top');
         grid off
        set(gcf,'Renderer','Painter')
        %hgexport(gcf,[pathFigs,num2str(day),'_x_',num2str(sourceLoP(iCross)*111),'_2D.eps']);
        %close  

    end
   end
    clear newZ tfIn zeq maxbyColum maxbycolum2 suma
end