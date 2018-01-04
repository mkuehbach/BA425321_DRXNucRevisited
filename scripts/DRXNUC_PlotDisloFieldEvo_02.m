%% DRXNUC read in fields binary files batch and visualize in animated gif
clear;
clc;
digits(32);
format long;

%get folder content
prefix = 'H:/Paper/Paper11_BambachDislDensModel/bb_simulation/DRXNUC/build/mpi01/Fields/';
fns = dir([prefix '*.bin']); %'*.RA.12.RB.16.*.bin']);

for f=1:length(fns(:,1))
    clearvars -except prefix fns f;
    nrows = 7551;
    ncols = str2double(extractBetween(fns(f).name,'NC.','.RhoFields'));
    %fnprefix = 'D:/DataTransfer/RELAX.tar/RELAX.SimID.1.EPSA.-300.EPSB.-300.T.1273.RA.13.RB.15';

    %% read DRXNUC binary dump file
    filename = [prefix fns(f).name]; %[prefix '.NR.' num2str(nrows) '.NC.' num2str(ncols) '.RhoFields.bin'];

    %adaptive NC
    %str2double(extractBetween(str,'NC.','.RhoFields'));

    fileID = fopen(filename);
    BIN = fread(fileID,[nrows,ncols],'float32');
    fclose(fileID);
    clearvars  filename fileID;

    %% extract data
    %first column is x spacing of simulation grid
    X = BIN(:,1)*1.0e6; %meter to micron
    %first row is time
    T = BIN(1,:);
    %delete unnecessary data
    BIN(:,1) = [];
    BIN(1,:) = [];
    X(1,:) = [];
    T(:,1) = [];
    if ~isempty(BIN)
        nrows = length(BIN(:,1));
        ncols = length(BIN(1,:));
        [fns(f).name ' loaded']

        fontszcap = 26;
        fontszax = 26;
        fontszlg = 20;
        fontnm = 'Calibri';

        fr = 1;
        fname = [prefix fns(f).name '.gif'];
        incr = int32(ncols/100);
        for c = 1:incr:ncols
            if c > ncols
                break;
            end
            %if fr > 5
            %    break;
            %end

            figure('Position',[100 100 1920 1080]);

            set(gca,'FontSize',fontszcap,'FontName',fontnm,'LineWidth',2);
            set(gcf,'PaperUnits','Inches');
            set(gcf,'PaperSize',[19.2 10.8]);
            if ( max(BIN(:,c)) < 1.0e18 ) 
                set(gcf,'color',[1 1 1],'visible','off');
            else
                set(gcf,'color',[1 0 0],'visible','off');
            end
            xlabel({'x ({\mu}m)'},'FontSize',fontszax,'FontName',fontnm);
            ylabel({'log_{10} \rho (m^{-2})'},'FontSize',fontszax,'FontName',fontnm);
            hold on
            plot(X,log10(BIN(:,c)),'LineWidth',2,'Color',[0 0 1]);
            delta = diff(BIN(:,c));
            xb = X(find(delta == max(delta)),1);
            hold on
            plot(xb, 13,'.','MarkerSize',45,'color',[255/255 201/255 14/255]); %log(max(rho(1,1:rangelimit)))./log(10),'.','MarkerSize',20,'color',[1 0 0])

            %left y axis limits
            xlim([0 15])
            xt = [0:1:15];
            xticks(xt);
            ylim([13 18.0]);
            yt = [13:1:18];
            yticks(yt);
            pbaspect([1.920 1.080 1.080])

            box on
            grid on
            legend([{'\rho(x,t) '} fns(f).name],'FontSize',fontszlg,'Location','northeast','Color',[1 1 1])
            legend('boxoff') %,'TextColor','blue')
            drawnow
            frame = getframe(gcf);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,65536,'nodither');
            if fr ~= 1 %all other frames
                imwrite(imind,cm,fname,'gif','WriteMode','append', 'DelayTime', 0.05);
            else %first time
                imwrite(imind,cm,fname,'gif', 'DelayTime',0.05,'Loopcount',inf);  
            end
            close(gcf);
            fr = c;
            %fr
        end

        [fns(f).name ' visualized']
        disp(['Processing file --> ' num2str(f)])
    else
        disp(['NO FIELD DATA FOR ' num2str(f)])
    end
end













%% read DRXNUC ascii dump file
filename = [fnprefix '.csv'];
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
rhoB = dataArray{:, 6};
rhoBA = dataArray{:, 7};

clearvars filename delimiter startRow formatSpec fileID dataArray ans;







plotfields = false;
if plotfields == true
    figure('Position',[100 100 500 500]);
    grid('on')
    box('on')
    view(0,90)
    fontszcap = 26;
    fontszax = 26;
    fontnm = 'Calibri';
    xlabel({'x (nm)'},'FontSize',fontszcap,'FontName',fontnm);
    ylabel({'\rho (m^{-2})'},'FontSize',fontszcap,'FontName',fontnm);
    %zlabel({'nill (m)'},'FontSize',fontszcap,'FontName',fontnm);
    pbaspect([2 1 1])
    %set(gca, 'YScale', 'log')
    set(gca,'FontSize',fontszcap,'FontName',fontnm,'LineWidth',1.5)
    set(gcf,'PaperUnits','Inches')
    set(gcf,'PaperSize',[30 30])
    set(gcf,'color','w')
    xlim([0 15]);
    xt = [0:1:15];
    xticks(xt);
    %ylim([1.0e14 1.0e16]);
    %yt = [1.0e14 1.0e15 1.0e16];
    %yticks(yt);
    %zlim([-1e-9 66e-9]);
    %zt = [-1e-9:10e-9:66e-9];
    %zticks(zt);
    for t=1:100 %ncols-1
        if mod(t,10) ~= 0
            continue;
        end
        hold on
        plot(X,BIN(:,t))
    end
end

%% plot v(x)
figure('Position',[100 100 500 500]);
grid('on')
box('on')
view(0,90)
fontszcap = 26;
fontszax = 26;
fontnm = 'Calibri';
xlabel({'x (m)'},'FontSize',fontszcap,'FontName',fontnm);
ylabel({'v (m/s})'},'FontSize',fontszcap,'FontName',fontnm);
%zlabel({'nill (m)'},'FontSize',fontszcap,'FontName',fontnm);
pbaspect([2 1 1])
%set(gca, 'YScale', 'log')
set(gca,'FontSize',fontszcap,'FontName',fontnm,'LineWidth',1.5)
set(gcf,'PaperUnits','Inches')
set(gcf,'PaperSize',[30 30])
set(gcf,'color','w')
%xlim([0 15]);
%xt = [0:1:15];
%xticks(xt);
%ylim([1.0e14 1.0e16]);
%yt = [1.0e14 1.0e15 1.0e16];
%yticks(yt);
%zlim([-1e-9 66e-9]);
%zt = [-1e-9:10e-9:66e-9];
%zticks(zt);
plot(x,v);



%% DEBUG
   %set(gca, 'YScale', 'log')
% hold on
% plot(x(1,1:rangelimit).*m2micron,log(rho(1,1:rangelimit))./log(10),'LineWidth',2)
% hold on
% plot(x(1,1:rangelimit).*m2micron,log(rho0(1,1:rangelimit))./log(10),'--','LineWidth',2)
% hold on
% %mark position of the boundary
% plot(gb.x*m2micron, ylow,'.','MarkerSize',30,'color',[255/255 201/255 14/255]); %log(max(rho(1,1:rangelimit)))./log(10),'.','MarkerSize',20,'color',[1 0 0])

% if ( mversion == 2016 )
%     hold on
%     yyaxis right
%     %len = length(res.x(1,1:fr));
%     plot(res.x(1,1:fr).*m2micron,log(res.v(1,1:fr)) ./ log(10),'LineWidth',1.5,'Color',[255/255 201/255 14/25])
%     ylabel({'lg Migration velocity (m/s)'},'FontSize',14,'FontName','CMU Bright')
%     ylim([-7 -4])
% else
%     %plot qualitative evolution of v into the same diagram
%     %fr=length(res.x(:));
%     hold on
%     yremapped = log(res.v(1:fr))./log(10) + ylow + abs(-10.0);
%     plot(res.x(1:fr).*m2micron,yremapped,'.','LineWidth',4.0,'Color',[255/255 201/255 14/25])
% end
    %store as png and read in again
    %print('tmp', '-dpng', '-r600');
    %[im,~] = imread('tmp.png');
    %}
    %work with (low)res gif directly

