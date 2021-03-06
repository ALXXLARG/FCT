function MODIFIEDSTR = cat_io_strrep(ORIGSTR,OLDSUBSTR,NEWSUBSTR)
% _______________________________________________________________________
% cat_io_strrep replace strings by other strings. It based on the strrep 
% and allows to use cellstrings that were replace by another string or a
% similar number of cellstrings depending on their input order.
%
% claim{1} = 'This is a good example';
% claim{2} = 'This is a bad example';
% new_claimA = cat_io_strrep(claim,{' good',' bad'},'n')
% new_claimB = cat_io_strrep(claim,{'good','bad'},{'great','acceptable'})
%
% See also strrep, strfind, regexprep.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: cat_io_strrep.m 1791 2021-04-06 09:15:54Z gaser $

  if nargin==0, help cat_io_strrep; return; end

  if iscell(ORIGSTR)
    MODIFIEDSTR = ORIGSTR; 
    for i=1:numel(ORIGSTR)
      MODIFIEDSTR{i} = cat_io_strrep(ORIGSTR{i},OLDSUBSTR,NEWSUBSTR);
    end
  else
    if iscell(OLDSUBSTR)
      if iscell(NEWSUBSTR)
        if numel(OLDSUBSTR)==numel(NEWSUBSTR)
          MODIFIEDSTR = ORIGSTR; 
          for i=1:numel(OLDSUBSTR) 
            MODIFIEDSTR = strrep(MODIFIEDSTR,OLDSUBSTR{i},NEWSUBSTR{i});
          end
        else
          error('cat_io_strrep:input',...
            'If multiple new strings were used, their number must be equal to the number of old strings.\n'); 
        end
      else
        MODIFIEDSTR = ORIGSTR; 
        for i=1:numel(OLDSUBSTR) 
          MODIFIEDSTR = strrep(MODIFIEDSTR,OLDSUBSTR{i},NEWSUBSTR);
        end
      end
    else
      MODIFIEDSTR = strrep(ORIGSTR,OLDSUBSTR,NEWSUBSTR);
    end
  end
end
