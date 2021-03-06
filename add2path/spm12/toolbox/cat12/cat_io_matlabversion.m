function year = cat_io_matlabversion(varargin)
% return Matlab version year
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: cat_io_matlabversion.m 1791 2021-04-06 09:15:54Z gaser $


  if strcmpi(spm_check_version,'octave')
    year = 20201;
    return
  end
  
  vers = version;
  [t1,t2,t3,year] = regexp(vers,'\(R.....\)');
  switch year{1}(end-1)
    case 'a',  year = [year{1}(3:end-2) '1'];
    case 'b',  year = [year{1}(3:end-2) '2'];
    otherwise, year = [year{1}(3:end-2) '0'];
  end
  year = str2double(year);
  
  if nargin==1
    year = varargin{1}==year;
  end
  if nargin==2
    year = varargin{1}<=year && year<=varargin{1};
  end
end