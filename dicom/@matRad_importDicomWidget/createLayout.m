function  this = createLayout(this,handleParent) 

h1 = this.widgetHandle;

%Create all handles
% LOGO 
h2 = axes(...
    'Parent',h1,...
    'CameraPosition',[0.5 0.5 9.16025403784439],...
    'CameraTarget',[0.5 0.5 0.5],...
    'CameraViewAngle',6.60861036031192,...
    'PlotBoxAspectRatio',[1 0.204545454545455 0.204545454545455],...
    'Colormap',[0.2422 0.1504 0.6603;0.250390476190476 0.164995238095238 0.707614285714286;0.257771428571429 0.181780952380952 0.751138095238095;0.264728571428571 0.197757142857143 0.795214285714286;0.270647619047619 0.21467619047619 0.836371428571429;0.275114285714286 0.234238095238095 0.870985714285714;0.2783 0.255871428571429 0.899071428571429;0.280333333333333 0.278233333333333 0.9221;0.281338095238095 0.300595238095238 0.941376190476191;0.281014285714286 0.322757142857143 0.957885714285714;0.279466666666667 0.344671428571429 0.971676190476191;0.275971428571429 0.366680952380952 0.982904761904762;0.269914285714286 0.3892 0.9906;0.260242857142857 0.412328571428571 0.995157142857143;0.244033333333333 0.435833333333333 0.998833333333333;0.220642857142857 0.460257142857143 0.997285714285714;0.196333333333333 0.484719047619048 0.989152380952381;0.183404761904762 0.507371428571429 0.979795238095238;0.178642857142857 0.528857142857143 0.968157142857143;0.176438095238095 0.549904761904762 0.952019047619048;0.168742857142857 0.570261904761905 0.935871428571428;0.154 0.5902 0.9218;0.146028571428571 0.609119047619048 0.907857142857143;0.13802380952381 0.627628571428572 0.897290476190476;0.124814285714286 0.645928571428571 0.888342857142857;0.111252380952381 0.6635 0.876314285714286;0.0952095238095238 0.679828571428571 0.859780952380952;0.0688714285714285 0.694771428571429 0.839357142857143;0.0296666666666667 0.708166666666667 0.816333333333333;0.00357142857142857 0.720266666666667 0.7917;0.00665714285714287 0.731214285714286 0.766014285714286;0.0433285714285715 0.741095238095238 0.739409523809524;0.096395238095238 0.75 0.712038095238095;0.140771428571429 0.7584 0.684157142857143;0.1717 0.766961904761905 0.655442857142857;0.193766666666667 0.775766666666667 0.6251;0.216085714285714 0.7843 0.5923;0.246957142857143 0.791795238095238 0.556742857142857;0.290614285714286 0.797290476190476 0.518828571428572;0.340642857142857 0.8008 0.478857142857143;0.3909 0.802871428571428 0.435447619047619;0.445628571428572 0.802419047619048 0.390919047619048;0.5044 0.7993 0.348;0.561561904761905 0.794233333333333 0.304480952380953;0.617395238095238 0.787619047619048 0.261238095238095;0.671985714285714 0.779271428571429 0.2227;0.7242 0.769842857142857 0.191028571428571;0.773833333333333 0.759804761904762 0.164609523809524;0.820314285714286 0.749814285714286 0.153528571428571;0.863433333333333 0.7406 0.159633333333333;0.903542857142857 0.733028571428571 0.177414285714286;0.939257142857143 0.728785714285714 0.209957142857143;0.972757142857143 0.729771428571429 0.239442857142857;0.995647619047619 0.743371428571429 0.237147619047619;0.996985714285714 0.765857142857143 0.219942857142857;0.995204761904762 0.789252380952381 0.202761904761905;0.9892 0.813566666666667 0.188533333333333;0.978628571428571 0.838628571428572 0.176557142857143;0.967647619047619 0.8639 0.164290476190476;0.961009523809524 0.889019047619048 0.153676190476191;0.959671428571429 0.913457142857143 0.142257142857143;0.962795238095238 0.937338095238095 0.126509523809524;0.969114285714286 0.960628571428571 0.106361904761905;0.9769 0.9839 0.0805],...
    'XTick',[0 0.2 0.4 0.6 0.8 1],...
    'XTickLabel',{  '0'; '0.2'; '0.4'; '0.6'; '0.8'; '1' },...
    'YTick',[0 0.5 1],...
    'YTickLabel',{  '0'; '0.5'; '1' },...
    'Position',[0.659070191431176 0.0833333333333333 0.320875113947129 0.13953488372093],...
    'ActivePositionProperty','position',...
    'LooseInset',[0.21882384176291 0.326703619171829 0.159909730519049 0.222752467617156],...
    'FontSize',8,...
    'SortMethod','childorder',...
    'Tag','axesDKFZLogo');

h3 = get(h2,'title');
set(h3,...
    'Parent',h2,...
    'Units','data',...
    'FontUnits','points',...
    'Color',[0 0 0],...
    'Position',[0.500002438371832 1.03819444444444 0.5],...
    'PositionMode','auto',...
    'String',blanks(0),...
    'Interpreter','tex',...
    'Rotation',0,...
    'RotationMode','auto',...
    'FontName','Helvetica',...
    'FontSize',11,...
    'FontAngle','normal',...
    'FontWeight','bold',...
    'HorizontalAlignment','center',...
    'HorizontalAlignmentMode','auto',...
    'VerticalAlignment','bottom',...
    'VerticalAlignmentMode','auto',...
    'EdgeColor','none',...
    'LineStyle','-',...
    'LineWidth',0.5,...
    'BackgroundColor','none',...
    'Margin',3,...
    'Clipping','off',...
    'XLimInclude','on',...
    'YLimInclude','on',...
    'ZLimInclude','on',...
    'Visible','on',...
    'HandleVisibility','off',...
    'BusyAction','queue',...
    'Interruptible','on',...
    'HitTest','on',...
    'PickableParts','visible');

h4 = get(h2,'xlabel');
set(h4,...
    'Parent',h2,...
    'Units','data',...
    'FontUnits','points',...
    'Color',[0.15 0.15 0.15],...
    'Position',[0.500000476837158 -0.224074068279178 0],...
    'PositionMode','auto',...
    'String',blanks(0),...
    'Interpreter','tex',...
    'Rotation',0,...
    'RotationMode','auto',...
    'FontName','Helvetica',...
    'FontSize',11,...
    'FontAngle','normal',...
    'FontWeight','normal',...
    'HorizontalAlignment','center',...
    'HorizontalAlignmentMode','auto',...
    'VerticalAlignment','top',...
    'VerticalAlignmentMode','auto',...
    'EdgeColor','none',...
    'LineStyle','-',...
    'LineWidth',0.5,...
    'BackgroundColor','none',...
    'Margin',3,...
    'Clipping','off',...
    'XLimInclude','on',...
    'YLimInclude','on',...
    'ZLimInclude','on',...
    'Visible','on',...
    'HandleVisibility','off',...
    'BusyAction','queue',...
    'Interruptible','on',...
    'HitTest','on',...
    'PickableParts','visible');

h5 = get(h2,'ylabel');

set(h5,...
    'Parent',h2,...
    'Units','data',...
    'FontUnits','points',...
    'Color',[0.15 0.15 0.15],...
    'Position',[-0.052462119942136 0.500000476837158 0],...
    'PositionMode','auto',...
    'String',blanks(0),...
    'Interpreter','tex',...
    'Rotation',90,...
    'RotationMode','auto',...
    'FontName','Helvetica',...
    'FontSize',11,...
    'FontAngle','normal',...
    'FontWeight','normal',...
    'HorizontalAlignment','center',...
    'HorizontalAlignmentMode','auto',...
    'VerticalAlignment','bottom',...
    'VerticalAlignmentMode','auto',...
    'EdgeColor','none',...
    'LineStyle','-',...
    'LineWidth',0.5,...
    'BackgroundColor','none',...
    'Margin',3,...
    'Clipping','off',...
    'XLimInclude','on',...
    'YLimInclude','on',...
    'ZLimInclude','on',...
    'Visible','on',...
    'HandleVisibility','off',...
    'ButtonDownFcn',blanks(0),...
    'BusyAction','queue',...
    'Interruptible','on',...
    'HitTest','on',...
    'PickableParts','visible');

h6 = get(h2,'zlabel');

set(h6,...
    'Parent',h2,...
    'Units','data',...
    'FontUnits','points',...
    'Color',[0.15 0.15 0.15],...
    'Position',[0 0 0],...
    'PositionMode','auto',...
    'String',blanks(0),...
    'Interpreter','tex',...
    'Rotation',0,...
    'RotationMode','auto',...
    'FontName','Helvetica',...
    'FontSize',10,...
    'FontAngle','normal',...
    'FontWeight','normal',...
    'HorizontalAlignment','left',...
    'HorizontalAlignmentMode','auto',...
    'VerticalAlignment','middle',...
    'VerticalAlignmentMode','auto',...
    'EdgeColor','none',...
    'LineStyle','-',...
    'LineWidth',0.5,...
    'BackgroundColor','none',...
    'Margin',3,...
    'Clipping','off',...
    'XLimInclude','on',...
    'YLimInclude','on',...
    'ZLimInclude','on',...
    'Visible','off',...
    'HandleVisibility','off',...
    'BusyAction','queue',...
    'Interruptible','on',...
    'HitTest','on',...
    'PickableParts','visible');

h7 = axes(...
    'Parent',h1,...
    'CameraPosition',[0.5 0.5 9.16025403784439],...
    'CameraTarget',[0.5 0.5 0.5],...
    'CameraViewAngle',6.60861036031192,...
    'PlotBoxAspectRatio',[1 0.301587301587302 0.301587301587302],...
    'Colormap',[0.2422 0.1504 0.6603;0.250390476190476 0.164995238095238 0.707614285714286;0.257771428571429 0.181780952380952 0.751138095238095;0.264728571428571 0.197757142857143 0.795214285714286;0.270647619047619 0.21467619047619 0.836371428571429;0.275114285714286 0.234238095238095 0.870985714285714;0.2783 0.255871428571429 0.899071428571429;0.280333333333333 0.278233333333333 0.9221;0.281338095238095 0.300595238095238 0.941376190476191;0.281014285714286 0.322757142857143 0.957885714285714;0.279466666666667 0.344671428571429 0.971676190476191;0.275971428571429 0.366680952380952 0.982904761904762;0.269914285714286 0.3892 0.9906;0.260242857142857 0.412328571428571 0.995157142857143;0.244033333333333 0.435833333333333 0.998833333333333;0.220642857142857 0.460257142857143 0.997285714285714;0.196333333333333 0.484719047619048 0.989152380952381;0.183404761904762 0.507371428571429 0.979795238095238;0.178642857142857 0.528857142857143 0.968157142857143;0.176438095238095 0.549904761904762 0.952019047619048;0.168742857142857 0.570261904761905 0.935871428571428;0.154 0.5902 0.9218;0.146028571428571 0.609119047619048 0.907857142857143;0.13802380952381 0.627628571428572 0.897290476190476;0.124814285714286 0.645928571428571 0.888342857142857;0.111252380952381 0.6635 0.876314285714286;0.0952095238095238 0.679828571428571 0.859780952380952;0.0688714285714285 0.694771428571429 0.839357142857143;0.0296666666666667 0.708166666666667 0.816333333333333;0.00357142857142857 0.720266666666667 0.7917;0.00665714285714287 0.731214285714286 0.766014285714286;0.0433285714285715 0.741095238095238 0.739409523809524;0.096395238095238 0.75 0.712038095238095;0.140771428571429 0.7584 0.684157142857143;0.1717 0.766961904761905 0.655442857142857;0.193766666666667 0.775766666666667 0.6251;0.216085714285714 0.7843 0.5923;0.246957142857143 0.791795238095238 0.556742857142857;0.290614285714286 0.797290476190476 0.518828571428572;0.340642857142857 0.8008 0.478857142857143;0.3909 0.802871428571428 0.435447619047619;0.445628571428572 0.802419047619048 0.390919047619048;0.5044 0.7993 0.348;0.561561904761905 0.794233333333333 0.304480952380953;0.617395238095238 0.787619047619048 0.261238095238095;0.671985714285714 0.779271428571429 0.2227;0.7242 0.769842857142857 0.191028571428571;0.773833333333333 0.759804761904762 0.164609523809524;0.820314285714286 0.749814285714286 0.153528571428571;0.863433333333333 0.7406 0.159633333333333;0.903542857142857 0.733028571428571 0.177414285714286;0.939257142857143 0.728785714285714 0.209957142857143;0.972757142857143 0.729771428571429 0.239442857142857;0.995647619047619 0.743371428571429 0.237147619047619;0.996985714285714 0.765857142857143 0.219942857142857;0.995204761904762 0.789252380952381 0.202761904761905;0.9892 0.813566666666667 0.188533333333333;0.978628571428571 0.838628571428572 0.176557142857143;0.967647619047619 0.8639 0.164290476190476;0.961009523809524 0.889019047619048 0.153676190476191;0.959671428571429 0.913457142857143 0.142257142857143;0.962795238095238 0.937338095238095 0.126509523809524;0.969114285714286 0.960628571428571 0.106361904761905;0.9769 0.9839 0.0805],...
    'XTick',[0 0.2 0.4 0.6 0.8 1],...
    'XTickLabel',{  '0'; '0.2'; '0.4'; '0.6'; '0.8'; '1' },...
    'YTick',[0 0.5 1],...
    'YTickLabel',{  '0'; '0.5'; '1' },...
    'Position',[0.296160877513711 0.836202830188679 0.229433272394881 0.147058823529412],...
    'ActivePositionProperty','position',...
    'LooseInset',[0.265885262648529 0.319135999073457 0.194300768858541 0.217592726640994],...
    'FontSize',8,...
    'SortMethod','childorder',...
    'Tag','axesMatRadLogo' );

h8 = get(h7,'title');
set(h8,...
    'Parent',h7,...
    'Units','data',...
    'FontUnits','points',...
    'Color',[0 0 0],...
    'Position',[0.500001710558694 1.03618421052632 0.5],...
    'PositionMode','auto',...
    'Interpreter','tex',...
    'Rotation',0,...
    'RotationMode','auto',...
    'FontName','Helvetica',...
    'FontSize',11,...
    'FontAngle','normal',...
    'FontWeight','bold',...
    'HorizontalAlignment','center',...
    'HorizontalAlignmentMode','auto',...
    'VerticalAlignment','bottom',...
    'VerticalAlignmentMode','auto',...
    'EdgeColor','none',...
    'LineStyle','-',...
    'LineWidth',0.5,...
    'BackgroundColor','none',...
    'Margin',3,...
    'Clipping','off',...
    'XLimInclude','on',...
    'YLimInclude','on',...
    'ZLimInclude','on',...
    'Visible','on',...
    'HandleVisibility','off',...
    'BusyAction','queue',...
    'Interruptible','on',...
    'HitTest','on',...
    'PickableParts','visible');

h9 = get(h7,'xlabel');
set(h9,...
    'Parent',h7,...
    'Units','data',...
    'FontUnits','points',...
    'Color',[0.15 0.15 0.15],...
    'Position',[0.500000476837158 -0.212280696264485 0],...
    'PositionMode','auto',...
    'Interpreter','tex',...
    'Rotation',0,...
    'RotationMode','auto',...
    'FontName','Helvetica',...
    'FontSize',11,...
    'FontAngle','normal',...
    'FontWeight','normal',...
    'HorizontalAlignment','center',...
    'HorizontalAlignmentMode','auto',...
    'VerticalAlignment','top',...
    'VerticalAlignmentMode','auto',...
    'EdgeColor','none',...
    'LineStyle','-',...
    'LineWidth',0.5,...
    'BackgroundColor','none',...
    'Margin',3,...
    'Clipping','off',...
    'XLimInclude','on',...
    'YLimInclude','on',...
    'ZLimInclude','on',...
    'Visible','on',...
    'HandleVisibility','off',...
    'BusyAction','queue',...
    'Interruptible','on',...
    'HitTest','on',...
    'PickableParts','visible');

h10 = get(h7,'ylabel');

set(h10,...
    'Parent',h7,...
    'Units','data',...
    'FontUnits','points',...
    'Color',[0.15 0.15 0.15],...
    'Position',[-0.0732804215064755 0.500000476837158 0],...
    'PositionMode','auto',...
    'Interpreter','tex',...
    'Rotation',90,...
    'RotationMode','auto',...
    'FontName','Helvetica',...
    'FontSize',11,...
    'FontAngle','normal',...
    'FontWeight','normal',...
    'HorizontalAlignment','center',...
    'HorizontalAlignmentMode','auto',...
    'VerticalAlignment','bottom',...
    'VerticalAlignmentMode','auto',...
    'EdgeColor','none',...
    'LineStyle','-',...
    'LineWidth',0.5,...
    'BackgroundColor','none',...
    'Margin',3,...
    'Clipping','off',...
    'XLimInclude','on',...
    'YLimInclude','on',...
    'ZLimInclude','on',...
    'Visible','on',...
    'HandleVisibility','off',...
    'BusyAction','queue',...
    'Interruptible','on',...
    'HitTest','on',...
    'PickableParts','visible');

h11 = get(h7,'zlabel');

set(h11,...
    'Parent',h7,...
    'Units','data',...
    'FontUnits','points',...
    'Color',[0.15 0.15 0.15],...
    'Position',[0 0 0],...
    'PositionMode','auto',...
    'String',blanks(0),...
    'Interpreter','tex',...
    'Rotation',0,...
    'RotationMode','auto',...
    'FontName','Helvetica',...
    'FontSize',10,...
    'FontAngle','normal',...
    'FontWeight','normal',...
    'HorizontalAlignment','left',...
    'HorizontalAlignmentMode','auto',...
    'VerticalAlignment','middle',...
    'VerticalAlignmentMode','auto',...
    'EdgeColor','none',...
    'LineStyle','-',...
    'LineWidth',0.5,...
    'BackgroundColor','none',...
    'Margin',3,...
    'Clipping','off',...
    'XLimInclude','on',...
    'YLimInclude','on',...
    'ZLimInclude','on',...
    'Visible','off',...
    'HandleVisibility','off',...
    'BusyAction','queue',...
    'Interruptible','on',...
    'HitTest','on',...
    'PickableParts','visible');

h12 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'String','RT Plan (SOPInstanceUID)',...
    'Style','text',...
    'Position',[0.319051959890611 0.697961527418892 0.152067622441369 0.0391903531438416],...
    'BackgroundColor',[0.502 0.502 0.502],...
    'Children',[],...
    'Tag','text11');

h13 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'String','RT Dose (SOPInstanceUID)',...
    'Style','text',...
    'Position',[0.660893345487694 0.697961527418892 0.126377724372255 0.0391903531438416],...
    'BackgroundColor',[0.502 0.502 0.502],...
    'Children',[],...
    'Tag','text10' );

h14 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'String','Directory',...
    'Style','text',...
    'Position',[0.0438756855575868 0.810661764705882 0.0722120658135283 0.0386029411764706],...
    'BackgroundColor',[0.502 0.502 0.502],...
    'Children',[],...
    'Tag','directory_label');

h15 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'String','Browse',...
    'Position',[0.913162705667276 0.768382352941176 0.0630712979890311 0.0404411764705882],...
    'Callback',@(hObject,event) browse_button_Callback(this,hObject,event),...
    'Children',[],...
    'Tag','browse_button' );

h16 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'String','Patient ID',...
    'Style','text',...
    'Position',[0.0447897623400366 0.674632352941177 0.0557586837294333 0.059917312661499],...
    'BackgroundColor',[0.502 0.502 0.502],...
    'Children',[],...
    'Tag','patient_label' );

h17 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'Style','listbox',...
    'Value',1,...
    'Position',[0.0446672743846855 0.511627906976744 0.22971741112124 0.186046511627907],...
    'BackgroundColor',[1 1 1],...
    'Callback',@(hObject,event) patient_listbox_Callback(this,hObject,event),...
    'Children',[],...
    'Tag','patient_listbox');

h18 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'String','CT ( ',...
    'Style','text',...
    'Position',[0.318223253501284 0.451622164800459 0.0464075578022706 0.0335917312661499],...
    'BackgroundColor',[0.502 0.502 0.502],...
    'Children',[],...
    'Tag','ct_label' );

h19 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'Style','listbox',...
    'Value',1,...
    'Position',[0.318140382862352 0.267441860465116 0.320875113947128 0.186046511627907],...
    'BackgroundColor',[1 1 1],... 'Callback',@(hObject,event) ctseries_listbox_Callback(this,hObject,event),...
    'Children',[],...
    'Tag','ctseries_listbox');

h20 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'String','RT Structure Set (SOPInstanceUID)',...
    'Style','text',...
    'Position',[0.66172205187702 0.446023542922768 0.179000580094473 0.0391903531438415],...
    'BackgroundColor',[0.502 0.502 0.502],...
    'Children',[],...
    'Tag','struct_label');

h21 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'Style','listbox',...
    'Value',1,...
    'Position',[0.660893345487694 0.265503875968992 0.320875113947129 0.186046511627907],...
    'BackgroundColor',[1 1 1],... 'Callback',@(hObject,event) rtseries_listbox_Callback(this,hObject,event),...
    'Children',[],...
    'Tag','rtseries_listbox');

h22 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'String','Import',...
    'Position',[0.181403828623519 0.253875968992248 0.0628988149498633 0.0406976744186047],...
    'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
    'Callback',@(hObject,event) import_button_Callback(this,hObject,event),...
    'Children',[],...
    'Enable','off',...
    'Tag','import_button');

h23 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'String','Cancel',...
    'Position',[0.181403828623519 0.195736434108527 0.0628988149498633 0.0406976744186047],...
    'Callback',@(hObject,event) cancel_button_Callback(this,hObject,event),...
    'Children',[],...
    'Tag','cancel_button' );

h24 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'String','Series Instance UID /',...
    'Style','radiobutton',...
    'Value',1,...
    'Position',[0.339283500455788 0.454421475739305 0.113118422143035 0.0345248349124318],...
    'BackgroundColor',[0.502 0.502 0.502],...
    'Callback',@(hObject,event) SeriesUID_radiobutton_Callback(this,hObject,event),...
    'Children',[],...
    'Tag','SeriesUID_radiobutton' );

h25 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'String','Series Number )',...
    'Style','radiobutton',...
    'Position',[0.44872826717494 0.452555268446741 0.0999864092152151 0.0345248349124318],...
    'BackgroundColor',[0.502 0.502 0.502],...
    'Callback',@(hObject,event) SeriesNumber_radiobutton_Callback(this,hObject,event),...
    'Children',[],...
    'Tag','SeriesNumber_radiobutton' );

h26 = uipanel(...
    'Parent',h1,...
    'Title','Resolution',...
    'BackgroundColor',[0.502 0.502 0.502],...
    'Tag','uipanel1',...
    'Clipping','off',...
    'Position',[0.0446672743846855 0.151162790697674 0.103919781221513 0.174418604651163] );

h27 = uicontrol(...
    'Parent',h26,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'String','y step:',...
    'Style','text',...
    'Position',[0.0727272727272727 0.337349397590362 0.718181818181818 0.253012048192771],...
    'BackgroundColor',[0.502 0.502 0.502],...
    'Tag','text8' );

h28 = uicontrol(...
    'Parent',h26,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'String','z step:',...
    'Style','text',...
    'Position',[0.0727272727272727 0.0602409638554217 0.718181818181818 0.253012048192771],...
    'BackgroundColor',[0.502 0.502 0.502],...
    'Tag','text9' );

h29 = uicontrol(...
    'Parent',h26,...
    'Units','normalized',...
    'Style','edit',...
    'Position',[0.463636363636364 0.653846153846154 0.4 0.217948717948718],...
    'BackgroundColor',[1 1 1],... 'Callback',@(hObject,event) resx_edit_Callback(this,hObject,event),...
    'Tag','resx_edit');

h30 = uicontrol(...
    'Parent',h26,...
    'Units','normalized',...
    'Style','edit',...
    'Position',[0.463636363636364 0.385542168674699 0.4 0.216867469879518],...
    'BackgroundColor',[1 1 1],... 'Callback',@(hObject,event) resy_edit_Callback(this,hObject,event),...
    'Tag','resy_edit');

h31 = uicontrol(...
    'Parent',h26,...
    'Units','normalized',...
    'Style','edit',...
    'Position',[0.463636363636364 0.120481927710843 0.4 0.216867469879518],...
    'BackgroundColor',[1 1 1],... 'Callback',@(hObject,event) resz_edit_Callback(this,hObject,event),...
    'Tag','resz_edit');

h32 = uicontrol(...
    'Parent',h26,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'String','x step:',...
    'Style','text',...
    'Position',[0.0727272727272727 0.614457831325301 0.381818181818182 0.253012048192771],...
    'BackgroundColor',[0.502 0.502 0.502],...
    'Tag','text7' );

h33 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'String','DICOM Import',...
    'Style','text',...
    'Position',[0.510968921389397 0.828849889012209 0.209323583180987 0.0808823529411765],...
    'BackgroundColor',[0.502 0.502 0.502],...
    'Children',[],...
    'ForegroundColor',[0.0980392156862745 0.305882352941176 0.615686274509804],...
    'Tag','text12',...
    'FontSize',18,...
    'FontName','Century',...
    'FontWeight','bold' );

h34 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'Max',32,...
    'Style','listbox',...
    'Value',1,...
    'Position',[0.660893345487694 0.513565891472868 0.320875113947129 0.186046511627907],...
    'BackgroundColor',[1 1 1],...
    'Callback',@(hObject,event) doseseries_listbox_Callback(this,hObject,event),...
    'Tag','doseseries_listbox');

h35 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
   'HorizontalAlignment','left',...
    'Max',2,...
    'Style','listbox',...
    'Value',1,...
    'Position',[0.318140382862352 0.513565891472868 0.320875113947128 0.186046511627907],...
    'BackgroundColor',[1 1 1],...
    'Callback',@(hObject,event) rtplan_listbox_Callback(this,hObject,event),...
    'Tag','rtplan_listbox');

h36 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'Style','edit',...
    'Position',[0.0447897623400366 0.762867647058823 0.857404021937843 0.0477941176470588],...
    'BackgroundColor',[1 1 1],...
    'Callback',@(hObject,event) dir_path_field_Callback(this,hObject,event),...
    'Tag','dir_path_field' );

h37 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
     'String',{  'Import patient name & DICOM meta information'; '                                         '; '                                         ' },...
    'Style','checkbox',...
    'Value',1,...
    'Position',[0.0446672743846855 0.415661785816824 0.228805834092981 0.0484496124031008],...
    'BackgroundColor',[0.502 0.502 0.502],...
    'Callback',@(hObject,event) checkPatientName_Callback(this,hObject,event),...
    'Tag','checkPatientName' );

h38 = uicontrol(...
    'Parent',h1,...
    'Units','normalized',...
    'String',{  'Use RT Dose grid'; '                                         '; '                                         ' },...
    'Style','checkbox',...
    'Position',[0.0446672743846855 0.369150157909848 0.228805834092981 0.0484496124031008],...
    'BackgroundColor',[0.502 0.502 0.502],...
    'Callback',@(hObject,event) checkbox3_Callback(this,hObject,event),...
    'Enable','off',...
    'Tag','checkbox3' );

  this.createHandles();
  
end

