%% DRXNUC read in track csv data batch and extract in which case a reduction in migration speed was detected
clear;
clc;
digits(32);
format long;

%get folder content
dbl.prefix = 'H:/Paper/Paper11_BambachDislDensModel/bb_simulation/DRXNUC/build/mpi20/run/';
dbl.fns = dir([dbl.prefix '*.csv']); %'*.RA.12.RB.16.*.bin']);
dbl.n = length(dbl.fns);

for f=1:dbl.n
    %% read DRXNUC ascii dump file
    filename = [dbl.prefix dbl.fns(f).name];
    delimiter = ';';
    startRow = 3;
    formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    t = dataArray{:, 1};
    x = dataArray{:, 2};
    v = dataArray{:, 3};
    p = dataArray{:, 4};
    rhoA = dataArray{:, 5};
    rhoB = dataArray{:,6};
    rhoBA = dataArray{:, 7};
    
    %at least there are data logged additional to header
    nv = length(v(:,1));
    if nv > 1   
        %get evolution GB time-speed profile by piecewise difference values
        dv(1,1) = nan; %error flag
        dv(2:nv,1) = diff(v(:,1)); %v_ti+1 - v_ti
        
        %which values are negative indicating the sought after reduction of velocity?
        where = find(dv(:,1) < eps); %eps numerically zero
        
        %where is somewhere along the profile, even small dips...
        
        %what is the absolute minimal difference?
        minn = min(dv(:,1));
        
        %if any reasonable where is given to me, otherwise add dummy
        if ~isempty(where)
            RES(f,1) = v(where(1),1); %the velocity at which acceleration firstly (in time) becomes negative
            RES(f,2) = where(1)/nv; %when is this
        else
            RES(f,2) = nan;
        end
        %plot(t,x,'.')
        %plot(t,v,'.')
        %plot(t,dv,'.')
%         figure
%         hold on 
%         plot(log10(t),log10(rhoA))
%         hold on
%         plot(log10(t),log10(rhoB))
%         hold on
%         plot(log10(t),log10(rhoBA))
%         legend('rA','rB','{\Delta BA}')
%         ylim([17 19])
%         rr = rhoB-rhoA;
%         rrr = rr-rhoBA;
%         rrr = rrr ./ rhoB;
%         min(rrr)
        %plot(t,rhoBA,'.')
        %RES(f,1) = min(diff(v));
    else
        %no data at all, add dummies
        RES(f,1) = nan; %error value
        RES(f,2) = nan;  
    end
    %how many logs in total? (at most 1401 expected...)
    RES(f,3) = nv;
    clearvars -except dbl f RES;
    disp(['Processing ' num2str(f) '---> ' dbl.fns(f).name])
end
disp('Batch processed')

% f = 460;
% filename = [prefix fns(f).name];
% delimiter = ';';
% startRow = 3;
% formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% t = dataArray{:, 1};
% %x = dataArray{:, 2};
% v = dataArray{:, 3};
% plot(t,v);
% legend(fns(f).name);

%% filter interesting cases
targ = find(RES(:,2) > 0.01 & RES(:,2) < 0.99);  %those with intermediate reduction or constancy in velocity
%targ = find(RES(:,2) > 0.00 & RES(:,2) < 1.00); %almost all
scaler = 1.0;
figure('Position',[100 100 scaler*1920 scaler*1080]);
colormap('jet');
fontszcap = 26;
fontszax = 26;
fontszlg = 20;
fontnm = 'Calibri';
set(gca,'FontSize',fontszcap,'FontName',fontnm,'LineWidth',2);
set(gcf,'PaperUnits','Inches');
set(gcf,'PaperSize',[19.2 10.8]);
set(gcf,'color','w','visible','on');
xlabel({'log_{10} t (s)'},'FontSize',fontszax,'FontName',fontnm);
ylabel({'log_{10} v (m/s)'},'FontSize',fontszax,'FontName',fontnm);
hold on
%'.','MarkerSize',45,'color',[255/255 201/255 14/255]); %log(max(rho(1,1:rangelimit)))./log(10),'.','MarkerSize',20,'color',[1 0 0])
%xlim([0 15]); xt = [0:1:15]; xticks(xt);
%ylim([13 18.0]);
%yt = [13:1:18];
%yticks(yt);
pbaspect([1.92 1.08 1.08]); box on; grid on;


%get colormap
cmap = jet;
i = 1;
for T = 800:50:1200
    temp = ((T-800)/(1200-800)*64)+1;
    disp(temp)
    if temp == 65
        temp = 64;
    end
    T2COLOR(i,1) = T+273;
    T2COLOR(i,2:4) = cmap(temp,1:3);
    disp(T)
    i = i+1;
end

for trg=1:length(targ(:,1))
    %trg=19;
    %if int32(targ(trg,1)) ~= 460
    %    continue;
    %end
    hold on
    filename = [prefix fns(int32(targ(trg))).name];
    %get properties from filename
    fn.EPSA = str2double(extractBetween(filename,'EPSA.','.EPSB'));
    fn.EPSB = str2double(extractBetween(filename,'EPSB.','.T'));
    fn.T = str2double(extractBetween(filename,'T.','.RA'));
    fn.rhoA = str2double(extractBetween(filename,'RA.','.RB'));
    fn.rhoB = str2double(extractBetween(filename,'RA.','.RB'));
    
    delimiter = ';';
    startRow = 3;
    formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    t = dataArray{:, 1};
    x = dataArray{:, 2};
    v = dataArray{:, 3};
    hold on
    plot(log10(t),log10(v),'LineWidth',2,'Color',T2COLOR(find(T2COLOR(:,1)==fn.T),2:4));
    %plot(x,log10(v),'LineWidth',1,'Color',T2COLOR(find(T2COLOR(:,1)==fn.T),2:4)); 
    disp([num2str(trg) ' ---> ' num2str(targ(trg,1)) ' ---> ' fns(int32(targ(trg,1))).name])
    clearvars filename delimiter startRow formatSpec fileID dataArray t v fn;
end
colorbar 
legend(num2str(targ))

for trg=1:length(targ(:,1))
    disp([num2str(trg) ' ---> ' num2str(targ(trg,1)) ' ---> ' fns(int32(targ(trg,1))).name])
end
