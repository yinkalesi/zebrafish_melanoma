%% getFullChain_v2_0
%  Version 2.0
%  Author: Adeyinka Lesi
%  Date: 7/19/19
%  Project: Tumor Growth, Logarithmic Continuum Form

% model for cdf

%% Version History
%  2.0: using getTransfomredDistribution_CTC_v1_0

function [chains] = getFullChain_v2_0(res,data,xcomps)
% chains{par}{t} = d(ccdf)/dx(1:xcomps)

% initialize
key = res.key;
pars = key.PARAMETERS;
parnames = fieldnames(key.PARAMETERS);
npar = length(parnames);
rates = key.RATES;
f0 = rates.size1_ratio;
nt = length(data.t);
chains = cell(npar,nt);

% calculate trajectories for primaries and largest possible meta
x0s = key.TIME_ZERO_SIZES;
c0s = key.TIME_ZERO_DIST;
xcompi = cell(1,length(xcomps));
for i = 1:length(xcomps)
    xcompi{i} = key.FIELD.indexOf(xcomps{i},key);
end
Np = length(x0s);
if(key.USING_METAGEN)
    t = res.t3(res.metagen_end_index:end);
else
    t = res.t3;
end
ntt = length(t)-1;
dt = t(end)/ntt;

% factors encapsulating time dependence
if(key.USING_METAGEN)
    cf = res.capacity_factor(res.metagen_end_index:end);
    cf = min(1,cf);
    cfm = 0.5*(cf(1:end-1)+cf(2:end));
    cr = res.complementary_ramp(res.metagen_end_index:end);
    crm = 0.5*(cr(1:end-1)+cr(2:end));
    nm = res.num_meta(res.metagen_end_index:end);
    nmm = 0.5*(nm(1:end-1)+nm(2:end));
    ntum = res.num_tumors(res.metagen_end_index:end);
else
    cf = res.capacity_factor;
    cf = min(1,cf);
    cfm = 0.5*(cf(1:end-1)+cf(2:end));
    cr = res.complementary_ramp;
    crm = 0.5*(cr(1:end-1)+cr(2:end));
    nm = res.num_meta;
    nmm = 0.5*(nm(1:end-1)+nm(2:end));
    ntum = res.num_tumors;
end
nmtum = max(0,ntum-Np);
nmtumm = 0.5*(nmtum(1:end-1)+nmtum(2:end));
nmtnot0 = find(nmtumm~=0);
lnnmtm = zeros(size(nmtumm));
lnnmtm(nmtnot0) = log(nmtumm(nmtnot0));

xps = getSizeTrajectory_v3_4(key,[x0s; 1],dt,t(end),cf);
xp1 = xps(end,:);
xps = xps(1:end-1,:);
xpm = 0.5*(xps(:,1:end-1)+xps(:,2:end));
delxp = xps(:,2:end)-xps(:,1:end-1);
xp1m = 0.5*(xp1(1:end-1)+xp1(2:end));
delxp1 = xp1(:,2:end)-xp1(:,1:end-1);
tm = 0.5*(t(1:end-1)+t(2:end));

% use interpolation to make piecewise continuous pdf and cdf
pxcs = cell(1,nt);
pxca = cell(1,nt);
hxcs = cell(1,nt);
for it = 1:nt
    % identify appropriate x and t values
    xcomp = xcomps{it};
    ixc = xcompi{it};
    % need to get pdf values around xcomps - but need to account for
    % distribution in both directions by averaging
    pxc2 = res.dist2(ixc,it);
    xc2 = res.xp(ixc)';
    ixcgt1 = find(ixc>1,1,'first');
    pxc1 = [zeros(length(1:ixcgt1-1),1); res.dist2(ixc(ixcgt1:end)-1,it)];
    pxc1a = [res.dist2(1,it)*ones(length(1:ixcgt1-1),1); res.dist2(ixc(ixcgt1:end)-1,it)];
    xc1 = [(res.x(1)-0.5)*ones(1,length(1:ixcgt1-1)) res.xp(ixc(ixcgt1:end)-1)]';
    ixcltm = find(ixc<length(res.xp),1,'last');
    pxc3 = [res.dist2(ixc(1:ixcltm)+1,it);zeros(length(ixcltm+1:length(ixc)),1)];
    xc3 = [res.xp(ixc(1:ixcltm)+1) (res.x(end)+0.5)*ones(1,length(ixcltm+1:length(ixc)))]';
    
    % calculate best pdf values through interpolation
    qxp = xcomp<xc2;
    pxc = qxp.*(pxc1+(pxc2-pxc1).*(xcomp-xc1)./(xc2-xc1))+...
        (1-qxp).*(pxc2+(pxc3-pxc2).*(xcomp-xc2)./(xc3-xc1));
    pxcs{it} = pxc;
    pxca{it} = qxp.*(pxc1a+(pxc2-pxc1a).*(xcomp-xc1)./(xc2-xc1))+...
        (1-qxp).*(pxc2+(pxc3-pxc2).*(xcomp-xc2)./(xc3-xc1));
    
    % calculate best cdf values
    hx = res.dist_cum2(end,it)-res.dist_cum2(:,it);
    hxc = hx(ixc)+(hx(ixc+1)-hx(ixc)).*(xcomp-res.x(ixc)')./(res.x(ixc+1)-res.x(ixc))';
    hxcs{it} = hxc;
end

% calculate weights
wp = cell(1,nt);
wm = cell(1,nt);
wj = zeros(1,ntt); % portion of metastasis generation due to secondary mets
wD = cell(1,nt);
pmetgen = zeros(1,ntt);

% calculate diffusivity and diffusion length for size 1
Din1 = 0.5*((1-cfm).*rates.growth(xp1m)...
    +(1-crm).*rates.death(xp1m)...
    +rates.shed(xp1m));
sqrDlen = zeros(1,length(xp1));
for i = 1:length(tm)
    sqrDlen(i+1) = sqrDlen(i)+4*Din1(i)*dt;
end
Dlen = sqrt(sqrDlen)/4; % 1/4 factor seems to be useful...

for i = 1:nt
    it = round(data.t(i)/dt);
    if(it == 0) % treat t=0 differently
        it = 1;
    end
    xcomp = xcomps{i};
    % calculate g(x-xp)
    gs = zeros(length(xcomp),Np);
    dgdp = zeros(length(xcomp),Np);
    sumpp = zeros(length(xcomp),1);
    for ip = 1:Np
        xpi = xpm(ip,it);
        Dpi = 0.5*((1-cfm(it)).*rates.growth(xpi)...
            +(1-crm(it)).*rates.death(xpi)...
            +rates.shed(xpi));
        gs(:,ip) = c0s(ip)*(exp(-0.25*(xcomp-xpi).^2/Dpi/(t(it)+dt))*0.5/sqrt(pi*Dpi*(t(it)+dt))...
            +exp(-0.25*(xcomp+xpi-2).^2/Dpi/(t(it)+dt))*0.5/sqrt(pi*Dpi*(t(it)+dt)));
        dgdp(:,ip) = c0s(ip)*(-(xcomp-xpi).*exp(-0.25*(xcomp-xpi).^2/Dpi/(t(it)+dt))*0.5/sqrt(pi*Dpi*(t(it)+dt))...
            -(xcomp+xpi-2).*exp(-0.25*(xcomp+xpi-2).^2/Dpi/(t(it)+dt))*0.5/sqrt(pi*Dpi*(t(it)+dt)));
        sumpp = sumpp+0.5*c0s(ip)*(erfc(0.5*(xcomp-xpi)./sqrt(Dpi*(t(it)+dt)))...
            +erfc(0.5*(xcomp+xpi-2)./sqrt(Dpi*(t(it)+dt))));
    end
    gsum = max(key.LIMIT_RESET_THRESHOLD,sum(gs,2));
    
    % calculate wp
    wp{i} = gs./gsum;
    wD{i} = dgdp./sumpp;
    
    % calculate wm
    p_xcomp = max(key.LIMIT_RESET_THRESHOLD,pxca{i});
    wm{i} = (1 - gsum./p_xcomp).*(xcomp<=xp1m(it)+Dlen(it));
    wm{i} = max(0,wm{i});
end

% calculate velocity and metastasis generation for each primary and size 1
% metastasis
kmins = zeros(Np,size(xpm,2));
vins = zeros(Np,size(xpm,2));
Dins = zeros(Np,size(xpm,2));
vine = zeros(Np,size(xps,2));
sf1 = sum(c0s.*rates.seed(xpm),1);
sf2 = sum(c0s.*rates.seed_parderiv(xpm,'seed1'),1);
sf3 = sum(c0s.*rates.seed_parderiv(xpm,'seed2'),1);
sf = struct('seed',sf1,'seed_deriv1',sf2,'seed_deriv2',sf3);
for ix = 1:Np
    % metastasis generation by primary
    if(key.USING_EXPLICIT_CTC)
        kmins(ix,:) = rates.meta_ss(xpm(ix,:),sf);
    else
        kmins(ix,:) = rates.meta(xpm(ix,:));
    end
    % this calculation is to get weights later (portion of metastasis
    % generation due to secondary metastasis (the actual calculation is
    % actually the tumors generated by primary))
    pmetgen = pmetgen+c0s(ix)*kmins(ix,:);
    vins(ix,:) = (1-cfm).*rates.growth(xpm(ix,:))...
        -(1-crm).*rates.death(xpm(ix,:))...
        -rates.shed(xpm(ix,:));
    Dins(ix,:) = 0.5*((1-cfm).*rates.growth(xpm(ix,:))...
        +(1-crm).*rates.death(xpm(ix,:))...
        +rates.shed(xpm(ix,:)));
    vine(ix,:) = (1-cf).*rates.growth(xps(ix,:))...
        -(1-cr).*rates.death(xps(ix,:))...
        -rates.shed(xps(ix,:));
end
if(key.USING_EXPLICIT_CTC)
    kmin1 = rates.meta_ss(xp1m,sf);
else
    kmin1 = rates.meta(xp1m);
end
vin1 = (1-cfm).*rates.growth(xp1m)...
    -(1-crm).*rates.death(xp1m)...
    -rates.shed(xp1m);
vine1 = (1-cf).*rates.growth(xp1)...
    -(1-cr).*rates.death(xp1)...
    -rates.shed(xp1);    
% calculate actual weights
inot0 = find(nmm>0);
wj(inot0) = max(0,1 - dt*pmetgen(inot0)./nmm(inot0));


% calculate integrands for parameter specific derivatives

vdins = cell(1,npar);
Ddins = cell(1,npar);
dxdps = cell(1,npar);
kmdins = cell(1,npar);
vdin1 = cell(1,npar);
Ddin1 = cell(1,npar);
dx1dp = cell(1,npar);
kmdin1 = cell(1,npar);
for ip = 1:npar
    parname = parnames{ip}; % will select appropriate derivatives using this name
    % calculate integrands for each parameter
    vdini = zeros(Np,size(xpm,2));
    Ddini = zeros(Np,size(xpm,2));
    vdrvvi = zeros(Np,size(xpm,2));
    dxdpi = zeros(Np,size(xps,2));
    kmdini = zeros(Np,size(xpm,2));
    for ix = 1:Np
        vxm = vins(ix,:);
        vxe = vine(ix,:);
        vd = (1-cfm).*rates.growth_parderiv(xpm(ix,:),parname)...
            -(1-crm).*rates.death_parderiv(xpm(ix,:),parname)...
            -rates.shed_parderiv(xpm(ix,:),parname);
        Dd = 0.5*((1-cfm).*rates.growth_parderiv(xpm(ix,:),parname)...
            +(1-crm).*rates.death_parderiv(xpm(ix,:),parname)...
            +rates.shed_parderiv(xpm(ix,:),parname));
        if(key.USING_EXPLICIT_CTC)
            kmd = rates.meta_ss_parderiv(xpm(ix,:),parname,sf);
        else
            kmd = rates.meta_parderiv(xpm(ix,:),parname);
        end
        if(strcmp(parname,'carrying_capacity') && ~isinf(pars.carrying_capacity))
            vdcc = rates.growth(xpm(ix,:)).*cfm/pars.carrying_capacity;
            vdcc2 = -rates.growth(xpm(ix,:)).*vxm/pars.carrying_capacity*Np.*(dxdpi(1,1:end-1)+dxdpi(1,2:end))*0.5; % bad approximation
        else
            vdcc = zeros(1,size(xpm,2));
            vdcc2 = zeros(1,size(xpm,2));
        end
        if(strcmp(parname,'death_ramp_center'))
            sfc = 2.29756;
            rw = max(pars.death_ramp_width,10*dt);
            rc = pars.death_ramp_center;
            tterm = 2*sfc/rw*(tm-rc);
            vdrc = rates.death(xpm(ix,:)).*vxm.*(1-tanh(tterm).^2)*sfc/rw;
        else
            vdrc = zeros(1,size(xpm,2));
        end
        if(strcmp(parname,'death_ramp_width'))
            sfc = 2.29756;
            rw = max(pars.death_ramp_width,10*dt);
            rc = pars.death_ramp_center;
            tterm = 2*sfc/rw*(tm-rc);
            vdrw = rates.death(xpm(ix,:)).*vxm.*(1-tanh(tterm).^2).*(tm-rc)*sfc/rw^2;
        else
            vdrw = zeros(1,size(xpm,2));
        end
        if(strcmp(parname,'death_start_time'))
            % treating death_start_time effect as a sharp ramp
            sfc = 2.29756;
            rw = dt*10;
            rc = pars.death_start_time;
            tterm = 2*sfc/rw*(tm-rc);
            vdst = rates.death(xpm(ix,:)).*vxm.*(1-tanh(tterm).^2)*sfc/rw;
        else
            vdst = zeros(1,size(xpm,2));
        end
        % finish calculations
        vdini(ix,:) = vd+vdcc+vdcc2+vdrc+vdrw+vdst;
        Ddini(ix,:) = Dd+0.5*(vdcc+vdcc2-vdrc-vdrw-vdst);
        vdrvvi(ix,:) = vdini(ix,:)./vxm.^2.*delxp(ix,:);
        sumi = zeros(1,size(xps,2));
        for is = 1:length(sumi)-1
            sumi(is+1) = sumi(is)+vdrvvi(ix,is);
        end
        dxdpi(ix,:) = sumi.*vxe;
        kmdini(ix,:) = kmd;
    end
    
    % calculate integrands for tumors starting at size 1
    
    vxm = vin1;
    vxe = vine1;
    vd = (1-cfm).*rates.growth_parderiv(xp1m,parname)...
        -(1-crm).*rates.death_parderiv(xp1m,parname)...
        -rates.shed_parderiv(xp1m,parname);
    Dd = 0.5*((1-cfm).*rates.growth_parderiv(xp1m,parname)...
        +(1-crm).*rates.death_parderiv(xp1m,parname)...
        +rates.shed_parderiv(xp1m,parname));
    if(key.USING_EXPLICIT_CTC)
        kmd = rates.meta_ss_parderiv(xpm(ix,:),parname,sf);
    else
        kmd = rates.meta_parderiv(xp1m,parname);
    end
    if(strcmp(parname,'carrying_capacity') && ~isinf(pars.carrying_capacity))
        vdcc = rates.growth(xp1m).*cfm/pars.carrying_capacity;
        vdcc2 = -rates.growth(xp1m).*vxm/pars.carrying_capacity*Np.*(dxdpi(1,1:end-1)+dxdpi(1,2:end))*0.5; % bad approximation
    else
        vdcc = zeros(1,size(xp1m,2));
        vdcc2 = zeros(1,size(xp1m,2));
    end
    if(strcmp(parname,'death_ramp_center'))
        sfc = 2.29756;
        rw = max(pars.death_ramp_width,10*dt);
        rc = pars.death_ramp_center;
        tterm = 2*sfc/rw*(tm-rc);
        vdrc = rates.death(xp1m).*vxm.*(1-tanh(tterm).^2)*sfc/rw;
    else
        vdrc = zeros(1,size(xp1m,2));
    end
    if(strcmp(parname,'death_ramp_width'))
        sfc = 2.29756;
        rw = max(pars.death_ramp_width,10*dt);
        rc = pars.death_ramp_center;
        tterm = 2*sfc/rw*(tm-rc);
        vdrw = rates.death(xp1m).*vxm.*(1-tanh(tterm).^2).*(tm-rc)*sfc/rw^2;
    else
        vdrw = zeros(1,size(xp1m,2));
    end
    if(strcmp(parname,'death_start_time'))
        % treating death_start_time effect as a sharp ramp
        sfc = 2.29756;
        rw = dt*10;
        rc = pars.death_start_time;
        tterm = 2*sfc/rw*(tm-rc);
        vdst = rates.death(xp1m).*vxm.*(1-tanh(tterm).^2)*sfc/rw;
    else
        vdst = zeros(1,size(xp1m,2));
    end
    % finish calculations
    vdin1i = vd+vdcc+vdcc2+vdrc+vdrw+vdst;
    Ddin1i = Dd+0.5*(vdcc+vdcc2-vdrc-vdrw-vdst);
    vdrvv1i = vdin1i./vxm.^2.*delxp1;
    sumi = zeros(1,size(xp1,2));
    for is = 1:length(sumi)-1
        sumi(is+1) = sumi(is)+vdrvv1i(is);
    end
    dx1dpi = sumi.*vxe;
    kmdin1i = kmd;
    
    % save results
    vdins{ip} = vdini;
    Ddins{ip} = Ddini;
    dxdps{ip} = dxdpi;
    kmdins{ip} = kmdini;
    vdin1{ip} = vdin1i;
    Ddin1{ip} = Ddin1i;
    dx1dp{ip} = dx1dpi;
    kmdin1{ip} = kmdin1i;
end

% linear interpolation to get N1 and N2

N1i = zeros(1,ntt+1);
N1i(1) = res.key.FIELD.init(1);
N2i = zeros(1,ntt+1);
N2i(1) = res.key.FIELD.init(2);
if(length(res.t2)>=length(res.t))
    % t2 are points time points from the data
    N1t = res.dist2(1,:)*(2-res.x(1))/(res.x(2)-res.x(1));
    N2t = res.dist2(2,:)*(3-res.x(2))/(res.x(3)-res.x(2));
    nrt = length(res.t2);
    it2i = [1 1+round(res.t2/dt)];
else
    % find appropriate initial index
    iz = find(res.t>=t(1),1,'first');
    % t are saved time points
    N1t = res.dist(1,iz:end)*(2-res.x(1))/(res.x(2)-res.x(1));
    N2t = res.dist(2,iz:end)*(3-res.x(2))/(res.x(3)-res.x(2));
    nrt = length(res.t(iz:end));
    it2i = [1 1+round(res.t(iz:end)/dt)];
end

for it3 = 1:nrt
    denom = max(dt,it2i(it3+1)-it2i(it3));
    N1i(it2i(it3)+1:it2i(it3+1)) = N1i(it2i(it3))+...
        (N1t(it3)-N1i(it2i(it3)))*(1:it2i(it3+1)-it2i(it3))/denom;
    N2i(it2i(it3)+1:it2i(it3+1)) = N2i(it2i(it3))+...
        (N2t(it3)-N2i(it2i(it3)))*(1:it2i(it3+1)-it2i(it3))/denom;
end
% if there are some edge values not interpolated:
N1i(it2i(end)+1:end) = N1i(it2i(end));
N2i(it2i(end)+1:end) = N2i(it2i(end));
N1m = 0.5*(N1i(1:end-1)+N1i(2:end));
N2m = 0.5*(N2i(1:end-1)+N2i(2:end));

% prepare meta variables
if(key.USING_EXPLICIT_CTC)
    shed_exp = pars.shed2;
    meta_par = pars.ctc1*pars.shed1./(pars.ctc1+pars.ctc2+sf.seed);
else
    shed_exp = pars.meta2;
    meta_par = pars.meta1;
end
% need to account for number of metasis that are greater than 1
N1vm = N1m-N1t(1)*exp(-(rates.growth(1)+rates.death(1))*tm);
metm = max(0,(nmtumm.^(1+shed_exp)-N1vm));
nmdsci = zeros(1,ntt);
not0m = find(metm>0);
nmrat = zeros(1,ntt);
nmrat(not0m) = nmtumm(not0m).^(1+shed_exp)./metm(not0m);
    
% caculate derivatives at each time
nmdpa = cell(1,npar); %change in primary metastasis due to change in meta parameters
nmdpax = cell(1,npar); % individual primaries separately
nmdpb = cell(1,npar); %change in primary metastasis due to change in non-meta parameters
nmdsa = cell(1,npar); %change in secondary metastasis due to change in meta parameters
nmdsb = cell(1,npar); %change in secondary metastasis due parameters that change total number of metastasis
nmdsc = cell(1,npar); %change in metastasis due to change in size1 metastasis
dnmpdp = cell(1,npar); %change in meta generation due to primaries, and due to meta params increasing secondary metastasis
dnmsdp = cell(1,npar); %change in meta generation due to secondaries

it2 = round(data.t(nt)/dt);
if(it2 == 0) % treat t=0 differently
    it2 = 1;
end
for ip = 1:npar
    parname = parnames{ip};
    % integrands for changes in metastasis distribution
    % caculate indices for t-delt(x)
    
    df0i = rates.size1_ratio_parderiv(parname);
    if(strcmp(parname,'meta1') && pars.meta1~=0)
        dlnkm1i = 1/pars.meta1;
    elseif(strcmp(parname,'shed1') && pars.shed1~=0)
        dlnkm1i = 1/pars.shed1;
    elseif(strcmp(parname,'ctc1') && pars.ctc1~=0)
        dlnkm1i = 1/pars.ctc1-1./(pars.ctc1+pars.ctc1+sf.seed);
    elseif(strcmp(parname,'ctc2') && pars.ctc2~=0)
        dlnkm1i = -1./(pars.ctc1+pars.ctc1+sf.seed);
    else
        dlnkm1i = 0;
    end
    if(strcmp(parname,'meta2') || strcmp(parname,'shed2'))
        dkm2i = 1;
    else
        dkm2i = 0;
    end
    
    % calculate changes due to each primary
    nmdpaxi = zeros(Np,it2);
    for ix = 1:Np
        nmdpaxi(ix,:) = kmdins{ip}(ix,1:it2)./vins(ix,1:it2).*delxp(ix,1:it2);
    end
    nmdpai = f0*sum(c0s.*nmdpaxi,1); %primary/meta params
    nmdpbi = nmm(1:it2)*df0i; %primary/non-meta params
    
    nmdsai = wj(1:it2).*nmm(1:it2)*f0.*(nmrat.*lnnmtm(1:it2)*dkm2i+dlnkm1i); %secondary/meta params
    
    % estimate change in number metastasis by summing nmdpai, nmdpbi,
    % nmdsai (effects already calculated)
    sumnmd1 = zeros(1,it2+1);
    for it3 = 1:it2
        sumnmd1(it3+1) = sumnmd1(it3) + nmdpai(it3) + nmdpbi(it3) + nmdsai(it3);
    end
    
    % need to estimate change in number of metastasis due to parameter at
    % each time for further calculations - initial estimate is based on
    % already calculated values
    dnmdpe = 0.5*(sumnmd1(1:end-1)+sumnmd1(2:end));
    
    % estimate change in secondary metastasis based on estimate of change
    % due to primary
    nmdsbi = zeros(1,it2);
    n0a = find(nmtumm(1:it2)>0);
    nmdsbi(n0a) = (1+shed_exp)*wj(n0a).*nmm(n0a)*f0.*nmrat(n0a)./nmtumm(n0a).*dnmdpe(n0a); %secondary/all params
    % add in changes due to variations in N1
    dkg1 = rates.growth_parderiv(1,parname);
    dkr1 = rates.death_parderiv(1,parname);
    dkr2 = rates.death_parderiv(2,parname);
    dks2 = rates.shed_parderiv(2,parname);
    if(f0>0)
        n1dpai = (nmdpai+nmdpbi+nmdsai+nmdsbi)/f0-(dkg1+dkr1)*N1m*dt+(dkr2+dks2)*N2m*dt;
    else
        n1dpai = zeros(1,it2);
    end
    sumen1d = zeros(1,it2+1); % the formula is missing a factor of exp(-t*(kr1+kg1)) that needs to be added later!
    for it3 = 1:it2
        if(key.USING_EXPLICIT_CTC)
            meta_par_it3 = meta_par(it3);
        else
            meta_par_it3 = meta_par;
        end
        sumen1d(it3+1) = sumen1d(it3) + n1dpai(it3).*exp(tm(it3)*(pars.growth1+pars.death1+meta_par_it3));
    end
    if(key.USING_EXPLICIT_CTC)
        meta_par_t = [(1.5*meta_par(1)-0.5*meta_par(2)) 0.5*(meta_par(1:end-1)+meta_par(2:end)) (1.5*meta_par(end)-0.5*meta_par(end-1))];
    else
        meta_par_t = meta_par;
    end
    dn1dpe = exp(-t.*(pars.growth1+pars.death1+meta_par_t)).*sumen1d;
    dn1dpem = 0.5*(dn1dpe(1:end-1)+dn1dpe(2:end));
    nmdsci(not0m) = -f0*wj(not0m).*nmm(not0m).*dn1dpem(not0m)./metm(not0m);
    
    % get contribution to dnmdp from new calculations
    sumnmd2 = zeros(1,it2+1);
    for it3 = 1:it2
        sumnmd2(it3+1) = sumnmd2(it3) + nmdsbi(it3) + nmdsci(it3);
    end
    
    % calculate integrands
    nmdpa{ip} = nmdpai; %primary/meta params
    nmdpax{ip} = nmdpaxi;
    nmdpb{ip} = nmdpbi; %primary/non-meta params
    nmdsa{ip} = nmdsai; %secondary/meta params
    nmdsb{ip} = nmdsbi; %secondary/all params
    nmdsc{ip} = nmdsci; %secondary/change in N1
    dnmpdp{ip} = sumnmd1;
    dnmsdp{ip} = sumnmd2;
end


% make the chains!
for ip = 1:npar
    for it = 1:nt

        xcomp = xcomps{it};
        pxc = pxca{it};
        hxc = hxcs{it};
        
        it2 = round(data.t(it)/dt);        
        
        if(it2 > 0) % treat t=0 differently since init condition never changes
            % find appropriate indices for xcomp in x1 trajectory
            dx1dpi = dx1dp{ip};
            ix1 = zeros(1,length(xcomp));
            for jxc = 1:length(xcomp)
                jxci = find(xp1>xcomp(jxc),1,'first');
                if(isempty(jxci) || jxci>it2+1)
                    ix1(jxc) = it2+1;
                else
                    ix1(jxc) = jxci;
                end
            end
            % spreading factor
            sfac = xp1(ix1)./(Dlen(ix1)+xp1(ix1));
            ix2 = zeros(1,length(xcomp));
            for jxc = 1:length(xcomp)
                jxci = find(xp1+Dlen>xcomp(jxc),1,'first');
                if(isempty(jxci) || jxci>it2+1)
                    ix2(jxc) = it2+1;
                else
                    ix2(jxc) = jxci;
                end
            end
            
            % 1) change due to how distribution is moving
            % 1a) use weights to get dxpdp
            wpi = wp{it}; %xcomp by Np
            wDi = wD{it}; %xcomp by Np
            wmi = wm{it}; %xcomp by 1
            % effect of change in primary position
            dxpdpi = dxdps{ip}(:,it2+1);
            wdxp = wpi*dxpdpi;
            % effect of change in diffusivity
            dlnDdpi = Ddins{ip}(:,it2)./Dins(:,it2);
            wdlnD = wDi*dlnDdpi;
            dlnD1 = Ddin1{ip}(ix1-1)'./Din1(ix1-1)'./hxc;

            dhdp1 = pxc.*((1-wmi).*wdxp+wmi.*(dx1dpi(ix1).*sfac)')...
            -0.5*hxc.*((1-wmi).*wdlnD+wmi.*dlnD1.*sfac');
            % 2) change in metastasis
            dnmpdpi = dnmpdp{ip}';
            dnmsdpi = dnmsdp{ip}';
            
            dhdp2 = dnmpdpi(it2-ix2+2)+dnmsdpi(it2-ix2+2);
            
            % finish it!
            chains{ip,it} = dhdp1+dhdp2;
            if(isempty(chains{ip,it}))
                chains{ip,it} = 0;
            end
%             disp(['ip: ' num2str(ip) ', it: ' num2str(it) ', ' num2str(chains{ip,it}')]);
        else
            chains{ip,it} = zeros(max(1,length(xcomp)),1);
        end
    end
end

end