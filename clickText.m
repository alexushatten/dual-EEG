function clickText(h,txt,interpreter)

%function clickText(h,txt,interpreter)

%

% Adds code to a figure object so that when you click on it, a text box

% appears with the desired text in it. When you click on the text box, it

% disappears.

%

% Required Inputs:

%  h   - figure object handle (vector or singleton)

%  txt - string of cell array of strings containing object labels

%

% Optional Inputs:

%  interpreter - The value of the text object's "interpreter" property.

%                {default: 'tex'}

%

% Author: David Groppe

% Mehtalab, 2012

% Renaud Marquis, FBMlab: added thickening of linewidth when clicked



if nargin<3,

   interpreter='tex'; 

end



hTparams=['set(ht,''backgroundcolor'',''w'',''horizontalalignment'',''center'',''verticalalignment'',''middle'',''interpreter'',''' interpreter ''',''buttondownfcn'',''delete(gcbo);'');'];

if iscell(txt)

    if length(h)~=length(txt)

        error('To the number of elements of h and txt are different.');

    end

    for a=1:length(h),

        set(h(a),'userdata',txt{a});
        set(h(a),'linewidth',0.5); % #RM@FBMlab,December 2017

%         bdfcn=['Cp = get(gca,''CurrentPoint''); '  'Xp=Cp(1,1);',  'Yp=Cp(1,2);',   'Zp=Cp(1,3);',  'dat=get(gcbo,''userdata'');',   'ht=text(Xp,Yp,Zp,sprintf(''%s'',dat));',  hTparams]; %#RM@FBMlab,December 2017
        bdfcn=['Cp = get(gca,''CurrentPoint''); '  'Xp=Cp(1,1);',  'Yp=Cp(1,2);',   'Zp=Cp(1,3);',  'dat=get(gcbo,''userdata'');',   'ht=text(Xp,Yp,Zp,sprintf(''%s'',dat)); if get(gco,''linewidth'')==0.5,set(gco,''linewidth'',3),else set(gco,''linewidth'',0.5),end; ',  hTparams]; %#RM@FBMlab,December 2017

        set(h(a),'buttondownfcn',bdfcn);

    end

else

    set(h,'userdata',txt);
    set(h,'linewidth',0.5); %#RM@FBMlab,December 2017
    
%     bdfcn=['Cp = get(gca,''CurrentPoint''); '  'Xp=Cp(1,1);',  'Yp=Cp(1,2);',   'Zp=Cp(1,3);', 'dat=get(gcbo,''userdata'');', 'ht=text(Xp,Yp,Zp,sprintf(''%s'',dat));',  hTparams]; %#RM@FBMlab,December 2017
    bdfcn=['Cp = get(gca,''CurrentPoint''); '  'Xp=Cp(1,1);',  'Yp=Cp(1,2);',   'Zp=Cp(1,3);', 'dat=get(gcbo,''userdata'');', 'ht=text(Xp,Yp,Zp,sprintf(''%s'',dat)); if get(gco,''linewidth'')==0.5,set(gco,''linewidth'',3),else set(gco,''linewidth'',0.5),end; ',  hTparams]; %#RM@FBMlab,December 2017

    set(h,'buttondownfcn',bdfcn);

end