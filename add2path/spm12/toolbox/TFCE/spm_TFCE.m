function spm_TFCE
% TFCE Toolbox wrapper to call TFCE functions
%_______________________________________________________________________
% $Id: spm_TFCE.m 103 2016-09-22 14:10:38Z gaser $

rev = '$Rev: 103 $';

SPMid = spm('FnBanner',mfilename,rev);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','TFCE Toolbox');
spm_help('!Disp','TFCE.man','',Fgraph,'Threshold-Free Cluster Enhancement Toolbox');

fig = spm_figure('GetWin','Interactive');
h0  = uimenu(fig,...
	'Label',	'TFCE',...
	'Separator',	'on',...
	'Tag',		'TFCE',...
	'HandleVisibility','on');
h1  = uimenu(h0,...
	'Label',	'Estimate',...
	'Separator',	'off',...
	'Tag',		'Estimate TFCE statistic',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.tfce_estimate'');',...
	'HandleVisibility','on');
h2  = uimenu(h0,...
	'Label',	'Results',...
	'Separator',	'off',...
	'Tag',		'Non-parametric inference',...
	'CallBack','[hReg xSPM SPM] = cg_tfce_results(''Setup'');',...
	'HandleVisibility','on');
h3  = uimenu(h0,...
	'Label',	'Check for updates',...
	'Separator',	'off',...
	'Tag',		'Check for updates',...
	'CallBack','cg_tfce_update(1);',...
	'HandleVisibility','on');
