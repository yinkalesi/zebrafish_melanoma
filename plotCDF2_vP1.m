function fit = plotCDF2_vP1(res1,res2,data,select,axis_lims,height,width,xlab,ylab,fontsizes,legpos,xtic,ytic)
%% plotCDF2_vP1
%  Version 2.0
%  Author: Adeyinka Lesi
%  Date: 10/31/18
%  Project: Tumor Growth, Logarithmic Continuum Form
%  function fit = measureFit_v1_0(res,data)
%  res: struct of model results
%  data: struct of experimental data
%  fit: 1xC vector with fit measure at each comparison time
%  (this version plots results)
% fitting at intermediate points rather than at data points

%% Version History
%  1.0: This is just measureFit_v3_6 with very minor changes
%  2.0: plotting two model distributions

rt = res1.t2;
% rxim = res.max_index;
% rx = res.xp(1:rxim);
rkey = res1.key;
rkey2 = res2.key;
% ra = rkey.CONV.v2a(rx);
rf = rkey.FIELD;
rf2 = rkey2.FIELD;
% rdx = rf.dxs(1:rxim);
% rdis = res.dist2(1:rxim,:);

% et = data.t;

if(isempty(rt))
    fit = Inf;
else
    fit = zeros(1,length(rt));
end

weights_E = [1 1/2 1/6 1/24]';

legend_text = cell(1,length(select)*3);
lowx = data.x(1);
h = figure;
h.Position = [962 42 width height];
colors = hsv(length(select));
% semilogx(1,1);
% hold on;

% analyze each time separately
for i = 1:length(rt)
    ea = data.a(data.selector{i});
    ex = data.x(data.selector{i});
    edis = data.dist(data.selector{i},i);
    ecdf = data.cdf(data.selector{i},i);
    
    % get cummulative plot
    imin = rf.indexOf(1,rkey);
    imin2 = rf2.indexOf(1,rkey2);
    imax = max(res1.max_index,rf.indexOf(ex(end),rkey));
    imax2 = max(res2.max_index,rf2.indexOf(ex(end),rkey2));
    gx = res1.x(imin:imax+1);
    gx2 = res2.x(imin2:imax2+1);
    gxp = res1.xp(imin:imax);
    gxp2 = res2.xp(imin2:imax2);
    gdis = res1.dist2(imin:imax,i);
    gdis2 = res2.dist2(imin2:imax2,i);
    gdx = rf.dxs(imin:imax);
    gdx2 = rf2.dxs(imin2:imax2);
    gcdf = zeros(size(gx));
    gcdf2 = zeros(size(gx2));
    for k = 1:length(gcdf)-1
        if(gdis(k)>0)
            gcdf(k+1) = gdis(k)*gdx(k)+gcdf(k);
        else
            gcdf(k+1) = gcdf(k);
        end
    end
    for k = 1:length(gcdf2)-1
        if(gdis2(k)>0)
            gcdf2(k+1) = gdis2(k)*gdx2(k)+gcdf2(k);
        else
            gcdf2(k+1) = gcdf2(k);
        end
    end
    
    if(isfield(data,'cdf_adj'))
        p_adj = data.cdf_adj(i);
    else
        p_adj = 0;
    end
    
    [bin] = discretize(gcdf(end)-[ecdf(end)+p_adj 0],gcdf);
    [bin2] = discretize(gcdf2(end)-[ecdf(end)+p_adj 0],gcdf2);
    
%     ielow = rf.indexOf(ex(1),rkey);
    ielow = 1;
    if(isnan(bin(1)))
        iglow = rf.indexOf(data.lowx,rkey)-imin+1;
    else
        iglow = bin(1);
    end
    irlow = min(ielow,iglow);
    irhigh = bin(2);

    if(isnan(irlow))
        xlow = data.lowx;
        irlow = rf.indexOf(xlow,rkey)-imin+1;
    else
        xlow = gx(irlow);
    end
    if(isnan(irhigh))
        xhigh = gx(end);
        irhigh = rf.indexOf(xhigh,rkey)-imin+1;
    else
        xhigh = gx(irhigh+1);
    end
    
    ielow2 = 1;
    if(isnan(bin2(1)))
        iglow2 = rf2.indexOf(data.lowx,rkey)-imin2+1;
    else
        iglow2 = bin2(1);
    end
    irlow2 = min(ielow2,iglow2);
    irhigh2 = bin2(2);

    if(isnan(irlow2))
        xlow2 = data.lowx;
        irlow2 = rf.indexOf(xlow2,rkey)-imin2+1;
    else
        xlow2 = gx2(irlow2);
    end
    if(isnan(irhigh2))
        xhigh2 = gx2(end);
        irhigh2 = rf2.indexOf(xhigh2,rkey2)-imin2+1;
    else
        xhigh2 = gx(irhigh2+1);
    end
    
    % first order correction of binning
%     cshift = 0;
%     eshift = 0;
%     xshift = gdx(irlow)/(gcdf(irlow+1)-gcdf(irlow))*cshift;
    xshift = 0; % may not be useful anymore (may not have ever been useful)
    
    lowx = min(lowx,xlow);
    rxp = gxp(irlow:irhigh);
    rx = gx(irlow:irhigh+1);
    rdx = gdx(irlow:irhigh);
    rdis = gdis(irlow:irhigh);
   
    lowx2 = min(lowx,xlow2);
    rxp2 = gxp2(irlow2:irhigh2);
    rx2 = gx2(irlow2:irhigh2+1);
    rdx2 = gdx2(irlow2:irhigh2);
    rdis2 = gdis2(irlow2:irhigh2);
    
    % for squares error calculation
    ras = rkey.CONV.v2a([rx(1)+xshift rx(2:end)]);
    rcdfs = gcdf(irhigh+1)-gcdf(irlow:irhigh+1);
    ecdfs = p_adj+ecdf(end)-ecdf;
    ecdfs1 = p_adj+ecdf(end)-[0; ecdf(1:end-1)];
    
    ras2 = rkey2.CONV.v2a([rx2(1)+xshift rx2(2:end)]);
    rcdfs2 = gcdf2(irhigh2+1)-gcdf2(irlow2:irhigh2+1);
    
    %         figure;
    if(any(select==i))
        eac = zeros(1,2*length(ea));
        eac(1:2:end-1) = ea;
        eac(2:2:end) = ea;
        ecdfc = zeros(1,2*length(ea));
        ecdfc(1:2:end-1) = ecdfs1;
        ecdfc(2:2:end) = ecdfs;
        %         semilogx(ea,ecdfs,'o',ras,rcdfs,ras2,rcdfs2,'--','MarkerEdgeColor',colors(i,:),...
        %             'Color',0.9*colors(i,:),'LineWidth',1.5);
%         semilogx(eac,ecdfc,'o:','MarkerEdgeColor',colors(i,:),...
%             'Color',0.9*colors(i,:),'LineWidth',0.75);
%         semilogx(ras,rcdfs,ras2,rcdfs2,'--','Color',0.9*colors(i,:),'LineWidth',1.5);
%         hold on;
        ii = find(select==i,1,'first');
        h(3*i-2) = semilogx(ea,ecdfs,'o','MarkerEdgeColor',colors(ii,:),...
            'Color',0.9*colors(ii,:),'LineWidth',1.5,'MarkerSize',6);
        hold on;
        h(3*i-1) = semilogx(ras,rcdfs,'Color',colors(ii,:),'LineWidth',1.5);
        h(3*i) = semilogx(ras2,rcdfs2,'--','Color',[1 0.5 1].*colors(ii,:)+[0 0.5 0].*(1-colors(ii,:)),'LineWidth',2.5);
        legend_text{3*i-2} = ['Data, ' num2str(data.t(i)) ' Days'];
        legend_text{3*i-1} = ['Model 1, ' num2str(rt(i)) ' Days'];
        legend_text{3*i} = ['Model 2, ' num2str(rt(i)) ' Days'];
    end
    if(length(ex)>2)

%         ormax = 1;
%         mom_ints = zeros(ormax+1,2);
%         for or = 0:ormax
%             mom_ints(or+1,1) = getDataMomentIntegral(edis,ex,or);
%             mom_ints(or+1,2) = getModelMomentIntegral(rdis,rxp,rdx,or);
%         end
        
        inds = rf.indexOf(0.5*(ex(2:end)+ex(1:end-1)),rkey) - imin - irlow + 2;
        if(any(inds<=0))
            warning('Data values out of bounds');
            err = Inf;
        else
            diff = ecdfs(1:end-1) - rcdfs(inds)';
%             semilogx(ras(inds),rcdfs(inds),'d',ras(inds),abs(diff),'-');
%             hold on;
%             lea = (ea);
%             dea1 = [lea(2)-lea(1); lea(2:end)-lea(1:end-1)];
%             dea = 0.5*(dea1(2:end)+dea1(1:end-1))/dea1(end);
%             edis_ave = edis(1:end-1)./dea;
            weights_diff = ones(size(diff));
            weights_diff = weights_diff/sum(weights_diff);
            err = diff'.^2*weights_diff/length(diff);
        end
        % penalty for not having enough metasatsis
%         penalty = (min(0,gcdf(end)-ecdf(end)-p_adj-gcdf(irlow))/ecdf(end)/0.05)^2;
%         err = err + penalty;
    else
        
        ormax = 1;
        mom_ints = zeros(ormax+1,2);
        for or = 0:ormax
            mom_ints(or+1,1) = getDataMomentIntegral(edis,ex,or);%edis*ex^or;
            mom_ints(or+1,2) = getModelMomentIntegral(rdis,rxp,rdx,or,rx(1)+xshift,rx(end));
        end
    end

    if(length(ex)>2)
        fit(i) = err;
    else
        % calculate moments
        tot1 = mom_ints(1,1);
        tot2 = mom_ints(1,2);
        
        mom_ave = [mom_ints(:,1)/tot1 mom_ints(:,2)/tot2];
        
        xave1 = mom_ave(2,1);
        xave2 = mom_ave(2,2);
        
        moms = zeros(size(mom_ave));
        moms(1,:) = [tot1,tot2];
        moms(2,:) = [xave1,xave2];
        for or = 2:ormax
            moms(or+1,:) = mom_ave(or+1,:) - [xave1^or xave2^or];
        end
        
        % calculate difference measure
        % when the model behaves badly-can get negative numbers and E
        % becomes imaginary - fixing that with abs
        E = abs(log(moms(:,2)./moms(:,1))).*weights_E(1:ormax+1);
        fit(i) = prod(1+E.^2)-1;
    end
    disp(['Finished ' num2str(i) ': ' num2str(fit(i))]);
end

hold off;
xlabel(xlab,'FontSize',fontsizes(2));
xticks(xtic);
ylabel(ylab,'FontSize',fontsizes(3));
yticks(ytic);
legend(h([3*select-1 3*select]),legend_text([3*select-1 3*select]),'Location',legpos,'FontSize',fontsizes(4));
legend boxoff;
axis(axis_lims);
set(gca,'Box','off');
set(gca,'FontSize',fontsizes(1));
end

function mom = getDataMomentIntegral(Nx,x,or)
    mom = Nx'*(x.^or);
end

function mom = getModelMomentIntegral2(p,xp,dx,or)
    f = p'.*(xp-0.5).^or;
    mom = f*dx';
end

function mom = getModelMomentIntegral(p,xp,dx,or,xlow,xhigh)
    x0 = xp(1)-0.5*dx(1);
    x1 = xp(end)-0.5*dx(end);
    
    ext1 = (xlow-x0)/dx(1);
    ext2 = (xhigh-x1)/dx(end);

    f = p'.*(xp-0.5).^or;
    if(length(f)>1)
        df1 = 2*(f(2)-f(1))/(dx(1)+dx(2));
        df2 = 2*(f(end)-f(end-1))/(dx(end)+dx(end-1));
        
        mom = f(1)*dx(1)*(1-ext1)+df1*dx(1)^2*ext1*(1-ext1)+...
            f(2:end-1)*dx(2:end-1)'+...
            f(end)*dx(end)*ext2+df2*dx(end)^2*ext2*(ext2-1);
    else
        mom = (ext2-ext1)*f*dx;
    end
end

%% Version History

%  1.1: plot results
%  1.2: getting xlow and xhigh from range that gives the proper number to
%  tumors - leads to better measure when tested agaist Iwata literature
%  values
%  1.3: no plots
%  1.4: wrote search algorithm to avoid using discretize method - it did
%  not end up being faster though
%  1.5: version doesn't need dist_cum2 pre-calculated - some errors fixed
%  here - USE THIS or 1.6 over previous versions
%  1.6: no plot, adjusted getModelMomentIntegral to 
%  2.3: removing lower limit at data.lowx - better for fish data
%  2.7: penalty for not having enough metastasis
%  2.8: penalty for not having enough metastasis + linked-data plotting
%  3.0: removes penalty for not having enough metastasis by moving data cdf down
%  3.1: linkeddata version of 3.0
%  3.2: changing fit function
%  3.4: unweighted fit and N-cdf plot