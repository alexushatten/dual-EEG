function Out = Rscript( Command, Options )
% Rscript: run R command through Rscript
%
% Rscript( Command, Options )
%
%  Inputs
% --------
% Command : R command to evaluate
% Options : Rscript options, with respect to the following syntax:
% 
%   /path/to/Rscript [--options] [-e expr [-e expr2 ...] | file] [args] 
%
%   --options accepted are 
%       --help              Print usage and exit 
%       --version           Print version and exit 
%       --verbose           Print information on progress 
%       --default-packages=list 
%                       Where 'list' is a comma-separated set 
%                         of package names, or 'NULL' 
%   or options to R, in addition to --slave --no-restore, such as 
%       --save              Do save workspace at the end of the session 
%       --no-environ        Don't read the site and user environment files 
%       --no-site-file      Don't read the site-wide Rprofile 
%       --no-init-file      Don't read the user R profile 
%       --restore           Do restore previously saved objects at startup 
%       --vanilla           Combine --no-save, --no-restore, --no-site-file 
%                         --no-init-file and --no-environ 
%  
%-------------------------------------------------------------------------
% Nota bene:
%       change PathToRscript such that it fits your local settings
%
%       'file' may contain spaces but not shell metacharacters 
%       Expressions (one or more '-e <expr>') may be used *instead* of 'file' 
%       See also  ?Rscript  from within R 
%
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, June 2018
%-------------------------------------------------------------------------

PathToRscript = '"C:\Program Files\R\R-3.4.1\bin\Rscript.exe"';

if nargin < 2
    Options = '';
end
CommandFull = ['!',PathToRscript,' ',Options,' -e ',Command];
Out = evalc(CommandFull);
% Out = regexprep(Rscript(Out),'\[\d\]','');

end

