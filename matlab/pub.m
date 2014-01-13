% pub: A script for creating html documentation
% Executing this script via MATLAB will create html documentation in the
% ./doc/matlab directory for the files listed in the fname cell array. The
% publish command utilizes a custom xsl file for rendering any equations
% with MathJax instead of LaTeX. For this reason, equations entered into
% the MATLAB documentation should be surrounded with \[ and \] for equation
% blocks and \( \) for inline equations. 
%
% This command is automatically run with the make matlab-doc command, see
% the root CMakeLists.txt file.

% Display a startup message
disp('Publishing m-files to html...');

% List the files to publish
fname{1} = [cd,filesep,'fem/volume_average/elem_length.m'];
fname{2} = [cd,filesep,'fem/volume_average/quad4.m'];
fname{3} = [cd,filesep,'fem/volume_average/tau_1.m'];
fname{4} = [cd,filesep,'fem/volume_average/energy.m'];
fname{5} = [cd,filesep,'fem/volume_average/parameters.m'];
fname{6} = [cd,filesep,'fem/volume_average/thermodynamics.m'];


% Define the publishing settings
opt.stylesheet = '../doc/matlab.xsl';   % include MathJax script to render equations
opt.outputDir  = '../doc/matlab';       % html output directory
opt.evalCode   = false;                 % disable code evaluatoin

% Loop through the files and publish them
for i = 1:length(fname);
    if exist(fname{i},'file')
        publish(fname{i}, opt);   
    else
        errordlg(['ERROR: The file does not exist: ',fname{i}]);
    end
end

% Display a startup message
disp('MATLAB documentation creation complete.');
