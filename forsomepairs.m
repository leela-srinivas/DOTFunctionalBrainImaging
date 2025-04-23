function reconXYslice_gui()
    %% 0) Setup
    clear; close all;
    load('NeuroDOT_Data_Sample_CCW1.mat','data','info');
    load('A.mat','A');               % measurements × voxels
    lambda = 0.1;                    % regularization
    
    %% 1) Compute log‑ratio data
    y = -log(bsxfun(@times, data, 1./mean(data,2)));
    
    %% 2) Define voxel grid (must match imageReconstruction)
    xBnds = [-70  70];   yBnds = [-30  30];   zBnds = [1 10];
    mmX=2; mmY=2; mmZ=2;
    [Y, X, Z] = meshgrid( ...
      yBnds(1):mmY:yBnds(2), ...
      xBnds(1):mmX:xBnds(2), ...
      zBnds(1):mmZ:zBnds(2) );
    [nX,nY,nZ] = size(X);
    
    %% 3) Ask user for a time interval
    maxT = size(y,2);
    prompt = sprintf('Enter start and end time indices [1 … %d] as [t0 t1]: ', maxT);
    tRange = input(prompt);
    tStart = max(1, min(tRange(1), maxT));
    tEnd   = max(1, min(tRange(2), maxT));
    if tEnd < tStart, tmp=tStart; tStart=tEnd; tEnd=tmp; end
    Tsel = tStart:tEnd;
    nFrames = numel(Tsel);
    
    %% 4) Precompute reconstructions over that interval
    fprintf('Reconstructing %d frames (this may take a moment)...\n', nFrames);
    xVol = zeros(nX,nY,nZ,nFrames);
    for k=1:nFrames
        y_t = y(:, Tsel(k));
        x3  = imageReconstruction(y_t, A, lambda);  % nX×nY×nZ
        xVol(:,:,:,k) = x3;
    end
    fprintf('Done.\n');
    
    %% 5) Build GUI
    hFig = figure('Name','DOT/fNIRS Reconstruction Browser','Color','w', ...
                  'Position',[200 200 700 600]);
    hAx  = axes('Parent',hFig,'Position',[0.10 0.25 0.85 0.70]);
    
    % Initial display
    curFrame = 1;  curZ = 1;
    imagesc(X(:,1,1), Y(1,:,1), squeeze(xVol(:,:,curZ,curFrame))','Parent',hAx);
    axis image; set(hAx,'YDir','normal');
    xlabel(hAx,'X (mm)'); ylabel(hAx,'Y (mm)');
    title(hAx,sprintf('t = %d / %d,  z = %.1f mm', Tsel(curFrame), maxT, Z(1,1,curZ)));
    colorbar('peer',hAx);
    
    % Time slider
    uicontrol('Style','text',  'Position',[100 65 120 15], 'String','Frame (time)');
    hTime = uicontrol('Style','slider','Min',1,'Max',nFrames,'Value',1, ...
        'SliderStep',[1/(nFrames-1) 5/(nFrames-1)], ...
        'Position',[100 40 200 20], ...
        'Callback',@updateImage);
    
    % Z‑slice slider
    uicontrol('Style','text','Position',[380 65 120 15], 'String','Z slice');
    hZ = uicontrol('Style','slider','Min',1,'Max',nZ,'Value',1, ...
        'SliderStep',[1/(nZ-1) 1/(nZ-1)], ...
        'Position',[380 40 200 20], ...
        'Callback',@updateImage);

    % Nested callback to update the image on slider moves
    function updateImage(~,~)
        curFrame = round(get(hTime,'Value'));
        curZ     = round(get(hZ,  'Value'));
        slice    = squeeze(xVol(:,:,curZ,curFrame))';
        imagesc(X(:,1,1), Y(1,:,1), slice, 'Parent', hAx);
        axis image; set(hAx,'YDir','normal');
        title(hAx,sprintf('t = %d,  z = %.1f mm', Tsel(curFrame), Z(1,1,curZ)));
        drawnow;
    end

end
