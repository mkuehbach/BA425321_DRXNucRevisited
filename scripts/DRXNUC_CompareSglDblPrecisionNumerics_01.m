%% DRXNUC comparing predictions from simulation with double against single precision
clear;
clc;
digits(32);
format long;

%% get folder contents
sgl.prefix = 'H:/Paper/Paper11_BambachDislDensModel/bb_simulation/DRXNUC/build/mpi01/Track/SinglePrecision/';
dbl.prefix = 'H:/Paper/Paper11_BambachDislDensModel/bb_simulation/DRXNUC/build/mpi01/Track/DoublePrecision/';
sgl.fns = dir([sgl.prefix '*.csv']);
dbl.fns = dir([dbl.prefix '*.csv']);

%% compare results case by case
sgl.n = length(sgl.fns(:,1));
dbl.n = length(dbl.fns(:,1));
if sgl.n ~= dbl.n
    disp('Different number of files');
end

for f=1:sgl.n
    %% read DRXNUC ascii results file from single precision simulation
    filename = [sgl.prefix sgl.fns(f).name];
    delimiter = ';';
    startRow = 3;
    formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    s.t = dataArray{:, 1};
    s.x = dataArray{:, 2};
    s.v = dataArray{:, 3};
    %at least one data log additional to header?
    s.n = length(s.t(:,1));
    if s.n > 1
        s.p = dataArray{:, 4};
        s.rhoA = dataArray{:, 5};  %%6 for rhoA and col 5 for rhoB -->UNFORTUNATELY RHOA AND RHOB IN OUTPUT JUMBLED UP... has been correct in code 2018/01/03
        s.rhoB = dataArray{:,6};
        s.rhoBA = dataArray{:, 7};
    end
    clearvars -except sgl dbl f s RES;
    
    %% read DRXNUC ascii results file from double precision simulation
    filename = [dbl.prefix dbl.fns(f).name];
    delimiter = ';';
    startRow = 3;
    formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    d.t = dataArray{:, 1};
    d.x = dataArray{:, 2};
    d.v = dataArray{:, 3};
    %at least one data log additional to header?
    d.n = length(d.t(:,1));
    if d.n > 1
        d.p = dataArray{:, 4};
        d.rhoA = dataArray{:, 5};
        d.rhoB = dataArray{:,6};
        d.rhoBA = dataArray{:, 7};
    end
    clearvars -except sgl dbl f s d RES;
    
    %% two files the same?
    
    
%     yi = interp1(x,Y,xi) returns vector yi containing elements corresponding to the elements of 
%     xi and determined by interpolation within vectors x and Y. The vector x specifies the 
%     points at which the data Y is given. If Y is a matrix, then the interpolation 
%     is performed for each column of Y and yi is length(xi)-by-size(Y,2).
     
    if s.n > 1 & d.n > 1
        %% interpolation
        s.xi = interp1(d.t,d.x,s.t);
        s.xdff = s.x - s.xi;
        RES(f).sxerr = max(abs(s.xdff));

        s.vi = interp1(d.t,d.v,s.t);
        s.vdff = s.v - s.vi;
        RES(f).sverr = max(abs(s.vdff));

        s.pi = interp1(d.t,d.p,s.t);
        s.pdff = s.p - s.pi;
        RES(f).perr = max(abs(s.pdff));

        s.ra = interp1(d.t,d.rhoA,s.t);
        s.radff = s.rhoA - s.ra;
        RES(f).raerr = max(abs(s.radff));

        s.rb = interp1(d.t,d.rhoB,s.t);
        s.rbdff = s.rhoB - s.rb;
        RES(f).rberr = max(abs(s.rbdff));

        s.rba = interp1(d.t,d.rhoBA,s.t);
        s.rbadff = s.rhoBA - s.rba;
        RES(f).rbaerr = max(abs(s.rbadff));
%         figure
%         hold on
%         plot(s.t,s.rhoBA,'.','MarkerSize',4)
%         hold on
%         plot(d.t,d.rhoBA)
    end
    clearvars -except sgl dbl f RES;
    disp(f);
end

tmp = nan(1,sgl.n);
for f=1:sgl.n
    if ~isempty(RES(f).raerr)
        tmp(1,f) = RES(f).raerr;
    end
end
plot(1:1:sgl.n,log10(tmp))