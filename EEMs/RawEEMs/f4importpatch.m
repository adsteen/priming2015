
%%%%%%%%
% Assumes Current Folder is /Users/asteen1/Dropbox/Drew/Priming
% effect/data/EEMS/Priming/RawEEMs
%%%%%%%%

% Check that working directory is correct
if ~strcmpi(pwd, '/Users/asteen1/Dropbox/Drew/Priming effect/data/EEMS/Priming/RawEEMs')
    error('Change your working directory to /Users/asteen1/Dropbox/Drew/Priming effect/data/EEMS/Priming/RawEEMs')
end
    


%RunCorrections.m
%Originally written by Laurel Larsen, 5/25/07
%Update for the 6/09 fluorescence workshop by Kaelin Cawley, 6/18/09
%Update for the 4/2010 fluorescence workshop by Rachel Gabor, Bailey Simone,
%and Mike SanClements
%Updated significantly by Kathleen Brannen 8/1/2013

%Make sure you hit save any time you make changes or else it won't be
%integrated into the code. 
%Run this code 1st and have CorrectFunII also open, as it will call the code.
%Set the current directory to the folder with your raw eems. Also have this
%code and CorrectfunII in that same folder. Also must create an input text 
%file that the code reads in (see CorrectionInputsIINFO.xls for a guide of what
%info needs to be included). 
%It saves the corrected EEMs.  The code also calculates and saves FI.
%Users should pay attention to the USER INPUT section in this code and in 
%CorrectFunII code!


%Notes on directories: Put this code, CorrectFun.m, and the instrument 
%correction files in the directory with the raw EEM files prior to running the code. 

clear
clc
close all

%USER INPUT
%%%%%%%
%
% There is a problem with theser paths
%
%%%%%%%
uvpath = '../UVFiles/'; %Enter the path for location of your UV files in .xlsx format
correctedpath = '../CorrectedEEMs/'; %~/Dropbox/Drew/Priming effect/data/EEMS/Priming/UVFiles/%Enter the path where the corrected EEMs will be saved
blankpath = '../FLBlank/'; %~/Dropbox/Drew/Priming effect/data/EEMS/Priming/FLBlank/%Enter the path where the blank files from the fluorometer are
ramanpath = '../FLRaman/'; %~/Dropbox/Drew/Priming effect/data/EEMS/Priming/FLRaman/%Enter the path where the raman files from the fluorometer are
EEMpath = ''; %~/Dropbox/Drew/Priming effect/data/EEMS/Priming/RawEEMs/
% All path locations should end with either "/" or "\" depending on your
% computer operating system

   

    %Read input file
    % for the input file, make an excel file (use the CorrectionInputsIINFO.xls
    % to figure out what goes in what column, then save it as a text file and 
    % input the path below.)
fid = fopen('CorrectionInputs.txt'); %path for correction inputs spreadsheet
text = textscan(fid,'%s%s%s%s%s%f32%f32%s%s','Delimiter','\t'); %or try %b
textlen = length(text{1,1}); %This should be the number of samples you have
n = 1;
%END USER INPUT

%END USER INPUT


%CODE
%RUN CORRECTION CODE ON EACH FILE
for n = 1:textlen

    Afile = sprintf('%s%s', EEMpath, char(text{1,1}(n)), '.dat'); %Name of exported raw EEM (first column)

    uvfile = sprintf('%s%s%s', uvpath, char(text{1,2}(n)), '.xlsx'); %Name of UV file WHICH IS READ WITH XLS

    ifile = sprintf('%s%s', char(text{1,3}(n))); %Name of corrected file to save
    
    rfile = sprintf('%s%s%s', ramanpath, char(text{1,4}(n)), '.dat'); %Name of raman file
    
    bfile = sprintf('%s%s%s', blankpath, char(text{1,5}(n)), '.dat'); %Name of blank file
    
    dilution_factor = text{1,6}(n); %Input dilution factor
    
    uvlength = text{1,7}(n); %Input uv pathlength (in cm)
   
    [FI(n), HIX(n), FrI(n), abs254(n), maxI(n), maxEm(n), maxEx(n), microfulvic(n)] = f4correction(Afile, uvfile,ifile, rfile, bfile, dilution_factor, correctedpath, uvlength); %Save out corrected EEM and fluorescence index as well as abs254
    fprintf('Progress: File number ')
    n;
end

%SAVE FLUORESCENCE INDEX TABLE

path = correctedpath;
filename = 'Indicespatch';
data =[FI' HIX' FrI' abs254' maxI' maxEm' maxEx' microfulvic']; %(FI; maxEm; abs254; HIX; FrI; SI75)
Ifile = sprintf('%s%s', correctedpath, filename, '.xlsx')
xlswrite(Ifile,data);

