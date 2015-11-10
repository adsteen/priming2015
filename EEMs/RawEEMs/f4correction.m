function [ifile, FI, HIX, FrI, abs254, maxI, maxEm, maxEx, microfulvic] = f4correction(Afile,uvfile,ifile,rfile,bfile,dilution_factor,correctedpath,uvlength)

%Make sure you hit save any time you make changes or else it won't yet be
%integrated into the code. 
%USER INPUT
%This is where you input the scan parameters that you used
eminc = 2.5; %the increment of the emission spectra you collected
exinc = 5; %the increment of the excitation spectra you collected
em = 250:eminc:550; %Emission start wavelength:eminc:emission ending wavelength
ex = 240:exinc:450; %Excitation start wavelength:exinc:excitation ending wavelength
RamanEnd = 450; %This is the end of your Raman scan, or where you want the 
                %raman scan to end the integration area
RamanBegin = 365; %Where you want the Raman scan to start integration
                   %Usually 370 but make sure its after your scan starts
                   % This should match the range for the mcRaman correction file input as "RC" in
                   % the RunCorrectionsII code
RamanInc = 1; %The increment on your raman scan  
SynchInc = 1; %Increment of Synch Scan
inc_em = eminc/0.5; % Increment of UV-Vis Absorbance scan
inc_ex = exinc/0.5; % Increment of UV-Vis Absorbance scan    

%ExTop=[0,ex]; %Recreates the row of excitation wavelengths cut off when inputting the .dat files
ExTop=ex; % My files don't seem to have an extra column of wavelengths as the first colunn

%CODE
emlen = length(em);
exlen = length(ex);

%Read in Raman file, instrument correct, Calculate Area under curve
R=importdata(rfile,'\t');%R=importdata(rfile,'\t');
R=R.data;

%Trims the old raman files if it's needed from 370 to scan end
Rfind = find(R == RamanBegin);
Rfindend = find(R == RamanEnd);
Rlen = length(R);
Raman = R(Rfind:Rfindend,2);


%The section below calculates the area under the raman curve.
y = Raman; %y is the trimmed Raman data vector
x = R(Rfind:Rfindend,1); %x is the wavelengths for the Raman vector
xlen = length(x) - 1; %One less than the max x value, because we integrate from x1 to xlen+1
summation = 0; 
iteration=1; %Unnecessary, I think

for i=1:xlen %This integrates from RamanBegin to RamanEnd.
    y0 = y(i); % The current Raman value
    y1 = y(i + 1); % The subsequent y value
    dx = x(i+1) - x(i); % dx
    summation = summation + dx * (y0 + y1)/2; %This is the right way to integrate (but a bit silly given that the x vector is defined as having integer values
    iteration = iteration+1; %Unneccessary, I think
end
%%%%%%%% CODE ERROR? %%%%%%%%%%%
BaseRect = (y(1)+y(xlen))/2*(x(xlen)-x(1)); %IS THIS SUPPOSED TO BE y(xlen) + y(xlen)?
RamanArea = summation - BaseRect;

%Read in blank file, instrument correct, Raman normalize
B=importdata(bfile,'\t');
B=B.data;

B=[ExTop; B]; %Adds the excitation wavelengths back into the Blank file

Bsize = size(B);

%emfind = B(:,1);
%emstart = find(emfind == em(1));
emstart = 1;
%emend = Bsize(1);
emend = Bsize(1);

%exfind = B(emstart-1,:);
%exstart = find(exfind == ex(1));
exstart = 1;
%exend = Bsize(2);
exend = Bsize(2);

B = B(emstart:emend,exstart:exend); % Cut B (blank) down to the desired size

Br=B/RamanArea; %This raman normalizes the corrected blank file.


%Read in sample file, instrument correct, IFE, Raman normalize, Blank subtract
% Reads in raw EEM file in .dat format.
%A=importdata(Afile,'\t');
A=importdata(Afile,'\t', 2);
A=A.data; %MY ADDITION

%Calculates the emission wavelegnth where the Max Intensity occurs
% maxEm = em(find(A(ex370, :) == max(A(ex370, :))));
%A=[ExTop ; A]; %Adds the excitation wavelengths back into the EEM file
A=[ExTop ; A(:, 2:end)]; %Cuts of emission wavelengths
maxA = max(A(:));
[row column]=find(A==maxA);
A(row,column);
maxEm= A(row,1);
maxEx= A(1, column);
Asize = size(A);
emfind = A(:,1);
%emstart = find(emfind == em(1));
emstart=2; % MY ADDITION
emend = Asize(1);
exfind = A(emstart-1,:);
%exstart = find(exfind == ex(1));
exstart = 1;
exend = Asize(2);
A = A(emstart:emend,exstart:exend);

Abs = xlsread(uvfile,1);  % Reads in the UV absorance file that has been transferred to csv file format
waves = Abs(:,1); %Waves=wavelengths
wave254 = find(waves == 254);
abs254 = Abs(wave254, 2);
exabsstart = find(waves == ex(1));
%%%%%%%%%%
% I MAY NEED TO MODIFY THE NEXT 3 LINES TO DEAL WITH DIFFERENT INSTRUMENTS
%%%%%%%%%%
exabsend = find(waves == ex(exlen)); %cuts 240 to 450 by 1 our scan is 240 to 450 by 10
emabsstart = find(waves == em(1)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THIS RETURNS 721
emabsend = find(waves == em(emlen)); %cuts 300 to  550 by 1 our scan is 300 to 550 by 2
absint = Abs(:,2)/uvlength;

% Performs Inner Filter Calculation using the UV absorbance spectrum.
ex_abs=absint(exabsstart:-exinc:exabsend,:);
em_abs=absint(emabsstart:-eminc:emabsend,:); % THrows a warning about integer indices - not sure if this is a problem
%ex_abs = absint(exabsstart:inc_ex:exabsend); %%%%%%%%THESE ARE NOW EMPTY 
%em_abs = absint(emabsstart:inc_em:emabsend);
IFC=zeros(length(em_abs),length(ex_abs));
for i=1:length(em_abs)
    for j=1:(length(ex_abs))
        IFC(i,j)=ex_abs(j)+em_abs(i);
    end
end

% WHy is IFC matrix almost twice as big as the A (sample) matrix?
Aci = A.*10.^(0.5*IFC); %This is the IFC.  See p56 of Principles of Fluorescence Spectroscopy by Lakowitz for the calculations

Acir = Aci/RamanArea; %This raman normalizes the instrument and IFE corrected EEM file.

Asub = Acir - Br; %This blank subtracts the corrected EEM file.

Adil = Asub.*dilution_factor; %This applies the dilution factor normalization.

%%%%Cutting out Rayleigh reflection
temp = ones(size(Adil));
for n = 1:size(em,2)
    temp(n,em(n)<=ex+10) = 0;
    temp(n,em(n)>=ex*2-20) = 0;
end

Adil = Adil.*temp;
%Anorm = Adil.*temp;    

Afinal= Adil./(max(Adil(:)));
%Afinal = Adil./(max(Anorm(:)));
maxI = max(Asub(:));



% Save the raman normalized and correceted EEM matrix (inner filter too).
pathname = correctedpath;

for i=1:length(ifile)

    pathname(length(pathname) + 1) = ifile(i);

end

pathnamelength = length(pathname);

pathname(pathnamelength + 1: pathnamelength + 5) = '.xlsx';

xlswrite(pathname, Adil);

% This next part calculates the fluorescence index

ex370 = find(ex == 370); %Index where excitation is 370
em470 = find(em == 470); %Index where emission is 470
em520 = find(em == 520); %Index where emission is 520

A=Afinal'; %Transposes corrected matrix for plotting and FI.

FI = A(ex370, em470)./A(ex370, em520); %Calculates the fluorescence index.

%This part calculates the humification index

%Specifies size of matrix of your EEM
Asize = size(A);
ylen = Asize(1);
xlen = Asize(2);
y = ex;
x = em;
xend = x(xlen);
yend = y(ylen);

%Creates a grid of the EEM pairs you need to interpolate for
[xi yi] = meshgrid(x(1):1:xend,y(1):1:yend);

zi = interp2(x, y, A, xi, yi, 'spline');

%micro:fulvic ratio
em370 = find(xi(1,:) == 370);
line370 = zi( : ,em370);
ex240= find(yi(:,1) == 240);
ex250= find(yi(:,1) ==250);
microline= line370(ex240:ex250, 1);
setzeros= zeros(21, 1); 
microtrapz= cumtrapz(microline);
microarea=sum(microtrapz);

em425 = find(xi(1,:) == 425);
line425 = zi( : ,em425);
fulvicline= line425(ex240:ex250, 1); 
fulvictrapz= cumtrapz(fulvicline); 
fulvicarea=sum(fulvictrapz);

microfulvic= microarea/(fulvicarea+microarea); 

%Finds the row with excitation of 254 nm and extracts that data
ex254 = find(yi(:,1) == 254);
line254=zi(ex254,:);

% %Find boundaries for red humic peaks
em435 = find(xi(1,:) == 435);
em480 = find(xi(1,:) == 480);
RedHum = line254(1,em435:em480);

% %find boundaries for blue humic peaks
em300 = find(xi(1,:) == 300);
em345 = find(xi(1,:) == 345);
BlueHum = line254(1,em300:em345);

% %Finds the area under the peaks
RedA = trapz(RedHum);
BlueA = trapz(BlueHum);

%This first code calculates HIX with the updated formula from Ohno(2002)
HIX = RedA/(RedA+BlueA);

%This calculates HIX using the original formula from Zsolnay(1999).  Ensure
%you have low concentrations of DOC to not have inner filter or
%concentration affects
% HIX = RedA/BlueA;

%This part calculates the freshness index

ex310 = find(ex == 310); %index where excitation is 310 nm
em380 = find(em == 380); %index where emission is 380 nm
em420 = find(em == 420); %index where emission is 420 nm
em435 = find(em == 435); %index where emission is 436 nm
%Calculates the Freshness Index
FrI = A(ex310,em380)/max(A(ex310,em420:em435));


% Graph ...

figure(1), clf; %, subplot(2,1,1); This subplot puts the 370ex emission curve on the bottom half of the EEM tiff file.

% set colormap

% set min and max intensity values (to be mapped to min and max colors in

% colormap)

% prevent auto-setting of caxis by changing caxis to manual control

% now draw the graph

contourf(em,ex,A,20); % with 20 contour lines

handle = gca;

set(handle,'fontsize', 14);

colormap(jet);

% caxis sets range for plotting. Change as necessary.

caxis([0, max(max(A))]);

caxis('manual');

H = colorbar('vert');

set(H,'fontsize',14);

xlabel('Emission Wavelength, nm','fontsize',12)

ylabel('Excitation Wavelength, nm','fontsize',12)

title (ifile, 'fontsize', 12); 

% Saves current object, this won't work if you close the figure first.
% This command saves the current object only.
% Enter the path name for the corrected EEM figure to be placed.

pathname = correctedpath;

for i=1:length(ifile)

    pathname(length(pathname) + 1) = ifile(i);

end

pathnamelength = length(pathname);

pathname(pathnamelength + 1: pathnamelength + 4) = '.tif';

saveas(gcf, pathname, 'tif');

close;



