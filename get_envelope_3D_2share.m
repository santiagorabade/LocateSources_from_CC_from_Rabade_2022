clear all
close all
period='1-5';
plotPair=0;
format long
alt_file=load('stAlt.txt'); %% Load station altitude file
stAlt=2.700-alt_file(:,2); %% 
zdummie=ones(17,1)*10000000;
slow=1;
for depth=0%[0,0.2,0.4,0.6,0.8,1.0,1.2,1.4]; %%% Define array of depths 
    sourceLa=39.663:0.0005:39.686; %% Create array of latitutes (initial:step:final)
    sourceLo=-111.235:-0.0005:-111.258; %% Create array of longitudes (initial:step:final)
    sources=combvec(sourceLo,sourceLa); %% Create all the combinations of lat lon
    sources=sources';
    load geometry.mat
    load MyColormap
    daylist=load('timelist0726'); %% load timelist, it can be a 5-min, hourly, day, it's used to read the files
    stlistAux=[1:5]; %% station list number
    stlist=combvec(stlistAux,stlistAux); %% Create all the combinations for source-receiver, this can be replaced with a file
    for iDay=1:length(daylist)
        day=[daylist(iDay,:)]; 
        num2str(day)
        paOut=[pwd,'/matfiles/']; %% save result
        paFigs=[pwd,'/figures/'] %% save figs
        iAux=1;
        for iFile=1:length(stlist) %% for cycle to extract info from sac file headers
            stlist(:,iFile);
            if stlist(1,iFile)~=stlist(2,iFile)
                st1=['00',num2str(stlist(1,iFile))];
                st2=['00',num2str(stlist(2,iFile))];
                pathData=[pwd,'/data/',st1(end-2:end),'-',st2(end-2:end),'/'];
                file=dir([pathData,'cor*',st1(end-2:end),'-',st2(end-2:end),'_',num2str(day),'*.sac_bp1to5Hz']);
                hd = rdSacHead([pathData,file(1).name]);
                if hd.dist>0.01
                    signal=rdSac([pathData,file(1).name]);
                    env=abs(hilbert(signal));
                    env=decimate(env,2); %% decimate envelope for smoothness
                    [maxAm maxTAux]=max(env);
                    dist(1,iAux)=hd.dist;
                    maxT(1,iAux)=(maxTAux-(hd.npts-1)/2)*hd.delta;
                    hd.maxT=(maxTAux-(hd.npts-1)/2)*hd.delta;
                    hd.fileName=file(1).name;
                    hd.path=pathData;
                    header(iAux)=hd;
                    iAux=iAux+1;
                end
            end            
        end
        t1=hd.b:hd.delta:hd.b*-1; %% time array signal
        t2=hd.b:hd.delta*2:(hd.b)*-1; %% time array decimate envelope
        [A B]=meshgrid(sourceLa,sourceLo); %% meshgrid for plotting
        for iEnv=1:length(header)
            signal=rdSac([header(iEnv).path,header(iEnv).fileName]);
            env=abs(hilbert(signal));
            env=decimate(env,2);
            scName=str2num(header(iEnv).fileName(5:7));
            rcName=str2num(header(iEnv).fileName(9:11));
            for iSource=1:size(sources,1)
                [distSc] = get_dist(sources(iSource,1),sources(iSource,2),header(iEnv).evlo,header(iEnv).evla); %% dist from Sc
                [distRc] = get_dist(sources(iSource,1),sources(iSource,2),header(iEnv).stlo,header(iEnv).stla); %% dist from Rc
                deltaDist(iSource,1)=sqrt(distRc^2+(depth-stAlt(rcName))^2)-sqrt(distSc^2+(depth-stAlt(scName))^2); %% dist to depth
                deltaTime(iSource,1)=deltaDist(iSource,1)*slow; %% calculate time 
                valueZIndex(iSource,1)=round(deltaTime(iSource,1)*(0.5/(header(iEnv).delta)))+(header(iEnv).npts-1)/4; %%find index for time
                valueZ(iSource,iEnv)=env(valueZIndex(iSource,1)); %% Extact amplitude
            end
            dif=max(valueZIndex)-min(valueZIndex); %% index of time window of signal
            header(iEnv).dif=max(valueZIndex)-min(valueZIndex);
            header(iEnv).rmsIn=rms(env(min(valueZIndex):max(valueZIndex))); %% rms for the window
            header(iEnv).rmsFull=rms(env); %% rms outside
            if min(valueZIndex)<dif
                header(iEnv).rmsOut=rms([env(1:min(valueZIndex))', env(max(valueZIndex):end)']);
            else
                header(iEnv).rmsOut=rms([env(min(valueZIndex)-dif:min(valueZIndex))', env(max(valueZIndex):max(valueZIndex)+dif)']);
            end
            toSurf=reshape(valueZ(:,iEnv)/max(valueZ(:,iEnv)),length(sourceLo),length(sourceLa)); %% reshape for plotting
            if plotPair==1 %% plot if plotPai==1
                figure
                subplot(1,2,1)
                s=surf(B,A,toSurf);
                view(2)
                c=colorbar;
                c.Label.String = 'Normalized amplitude';
                caxis([0.4 1])
                shading interp
                ylim([sourceLa(1) sourceLa(end)])
                xlim([sourceLo(end) sourceLo(1)])
                s.EdgeColor = 'none';
                hold on
                plot3(header(iEnv).stlo,header(iEnv).stla,zdummie(1,1)+100,'k^','MarkerFaceColor', 'w','MarkerSize',8)
                hold on
                box on
                plot3(header(iEnv).evlo, header(iEnv).evla,zdummie(1,1)+100,'kp','MarkerFaceColor', 'w','MarkerSize',8)
                title([num2str(day),' ',header(iEnv).fileName(5:11)])             
                [signal,~]=rdSac([header(iEnv).path,header(iEnv).fileName]);
                envelope=abs(hilbert(signal));
                subplot(1,2,2)
                plot(t1,signal/max(signal),'k')
                hold on
                plot(t1,envelope/max(envelope),'b')
                hold on
                plot([min(deltaTime) min(deltaTime)],[-1.5, 1.5],'--k')
                hold on
                plot([max(deltaTime) max(deltaTime)],[-1.5, 1.5],'--k')
                hold on
                title(['maxT= ',num2str(max(deltaTime))])
                xlim([-3 3])
                ylim([-0.6 1.1])
                box on
                xlabel('Time (s)')
                fig = gcf;
                fig.PaperUnits = 'inches';
                fig.PaperPosition = [0 0 9 3];
                %print([paFigs,header(iEnv).fileName(5:11),'_',num2str(day),'.png'],'-dpng','-r0')
                %print('-painters','-dpdf',[paFigs,header(iEnv).fileName(5:11),'_',num2str(day)])
                close
            end
        end
        save([paOut,'valueZ_',period,'_slow_',num2str(slow),'_day_',num2str(day),'_dep_',num2str(depth),'_newCoord_all_16.mat'], 'valueZ') %%save header
        save([paOut,'header_',period,'_slow_',num2str(slow),'_day_',num2str(day),'_dep_',num2str(depth),'_newCoord_all_16.mat'], 'header') %%save matrix with location indivudal
        clear header
        clear valueZ
    end
    
end

