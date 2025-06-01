
% some details regarding conversions etc 
D2si = 0.01^2/60; % convert creep slopes fitted to SI units (measurements are in cm/min)

%% read data

fnames = dir('*.dat');
numfiles = length(fnames);

lidweightarr = NaN(numfiles,1);
sizearr = NaN(numfiles,1);
addedweightarr = NaN(numfiles,1);
lengarr = NaN(numfiles,1);

timearr = NaN(numfiles, 50); % make enough space
displarr = NaN(numfiles, 50);
errarr = NaN(numfiles, 50);

% for creep curves
slopearr = NaN(numfiles,1);
offsetarr = NaN(numfiles,1);
errorarr = NaN(numfiles,2);

% for lin curves
slopearr_lin = NaN(numfiles,1);
offsetarr_lin = NaN(numfiles,1);
errorarr_lin = NaN(numfiles,2);
conf_lin = NaN(numfiles,2,2);

%settings to extract best fit exponent
bexpy = 0.2; % lower exp
expyr = 0.025; % step exponent
eexpy = 1.4; % max exp
expyarr = NaN(numfiles,1+round((eexpy-bexpy)/expyr)); % save all exp fit rsquares here
stdarr = NaN(numfiles,1+round((eexpy-bexpy)/expyr)); % stats of residuals
kurtarr = NaN(numfiles,1+round((eexpy-bexpy)/expyr)); % kurtosis stats of residuals
skewarr = NaN(numfiles,1+round((eexpy-bexpy)/expyr)); % skewness stats of residuals
distopyarr = NaN(numfiles,1+round((eexpy-bexpy)/expyr)); % distance to optimum

bestexparr = NaN(numfiles,1); %store best exponent (as defined below) here
optexparr = NaN(numfiles,1); %another definition of "optimal" exponent


typearr = []; % save type: cyl or sphere. Not very relevant here as all cases were sphere. Type is used for sth else.
lidonarr = []; % extract whether lid is on or not

% extract prefactors for all nonlin fits (new for these data)
fitexpyarr = bexpy:expyr:eexpy;
slopearr_nonlin = NaN(numfiles,length(fitexpyarr));
offsetarr_nonlin = NaN(numfiles,length(fitexpyarr));
exparr_nonlin = NaN(numfiles,length(fitexpyarr));
conf_nonlin = NaN(numfiles,2,2);

for i=1:numfiles
    
    if ~mod(i,2), fprintf('numfile # %d\n', i); end
    
    % units all in centimeters and minutes!!!!!!!!!!
    readfile = char(fnames(i).name);
    
    lidweightarr(i) = str2num(readfile(2:4));
    typearr = [typearr, readfile(7)];
    lidonarr = [lidonarr, readfile(1)];
    sizearr(i) = str2num(readfile(5:6));
    addedweightarr(i) = str2num(readfile(8 :9));
    
    tmp = importdata(readfile,' ');
    tmpdata = tmp.data;
    tmpstring = tmp.textdata;
    tmpstringerr = tmpstring{3};
    datl = size(tmpdata);
    
    lengarr(i) = datl(1);
    timearr(i,1:datl) = tmpdata(:,1);
    displarr(i,1:datl) = tmpdata(:,2);
    
    if contains(tmpstringerr, 'err') % if there is an error estimate, it will be in the third column.
        fprintf('error estimate in file # %d\n', i)
        errarr(i,1:datl(1)) = tmpdata(:,3);
    end
    
    %%%%%%%%%%%% we first assume square root behavior:
    % delta = sqrt(Dx)
    % delta^2 = Dx
    % x = time
    % delta = displacement

    xtmp = [timearr(i,1:lengarr(i))];
    deltatmp = [displarr(i,1:lengarr(i)).^2]; % square first then do linear fit. Prefactor is D
    
    % fit poly1, offset should be zero hence added point at origin
    [fitobj,gof2] = fit(xtmp',deltatmp','poly1');
    
    tmp = coeffvalues(fitobj);    
    %prefactor is slope on y^2(t) plot
    slopearr(i) = tmp(1);
    offsetarr(i) = tmp(2);
    errorarr(i,1) = gof2.rmse;%/datl(1);
    errorarr(i,2) = gof2.rsquare;
    
    %%%%%%%%%%%%%%%% now fit the linear parts. may not exist for lid
    %%%%%%%%%%%%%%%% experiments
    deltatmp = [displarr(i,1:lengarr(i))];
    
    % fit poly1, offset should be zero hence added point at origin
    [fitobj,gof2] = fit(xtmp',deltatmp','poly1');
    
    %figure
    %plot(fitobj,xtmp',deltatmp')
    
    tmp = coeffvalues(fitobj);    
    slopearr_lin(i) = tmp(1);
    offsetarr_lin(i) = tmp(2);
    errorarr_lin(i,1) = gof2.rmse;%/datl(1);
    errorarr_lin(i,2) = gof2.rsquare;
    conf_lin(i,:,:) = confint(fitobj,0.95); 
    
    %%%%%%%%%%%%%%%% fit a range of different exponents just to see which
    %%%%%%%%%%%%%%%% might work best
    xtmp = timearr(i,2:lengarr(i));
    deltatmp = displarr(i,2:lengarr(i));
    ctr = 0;
    for expy = bexpy:expyr:eexpy
        ctr = ctr + 1;
        opty = fitoptions('power2','Lower', [0,expy,-10], 'Upper', [Inf,expy,10]);
       [fitobj,gof2] = fit(xtmp',deltatmp','power2',opty);
       expyarr(i,ctr) = gof2.rsquare;
       %thought: check stats of residuals
       resi = fitobj(xtmp)-deltatmp';
       stdarr(i,ctr) = nanstd(resi);
       kurtarr(i,ctr) = kurtosis(resi);
       skewarr(i,ctr) = skewness(resi);

       coiffure = coeffvalues(fitobj);
       exparr_nonlin(i,ctr) = expy; 
       slopearr_nonlin(i,ctr) = coiffure(1);
       offsetarr_nonlin(i,ctr) = coiffure(3);
    end
    [~,indexy] = max(expyarr(i,:));
    bestexparr(i) = bexpy+(indexy-1)*expyr;
    
    distopyarr(i,:) = sqrt(skewarr(i,:).^2 + (kurtarr(i,:)-3).^2 + (gradient(expyarr(i,:)).^2));
    [~,indexy] = min(distopyarr(i,:));
    optexparr(i) = bexpy+(indexy-1)*expyr;
end

fprintf('done\n')

close all

%% some admin stuff

%extract locations of sphere of 23mm. Trivially always the case for these experiments
masksp = find(sizearr == 23);

% which weights are used to add stress to sphere. Extracted from filenames
addweiuniq = unique(addedweightarr);

%what range of weights are used to add stress to sphere
minlid = min(unique(lidweightarr));
maxlid = max(unique(lidweightarr));
lidrange = maxlid-minlid;

% calculate stresses per experiment
denssphere = 913;  %polypropylene
denswater = 1006;  
weightrodsp = 0.00224; % rod for spheres is longer than for cylinders
weighttray = 0.002; % newly printed tray is lower in weight?
weightplunger23 = (4*3.1415*0.0115^3)*0.333*(denssphere-denswater); %(4/3)*pi*r^3*Delta-rho
surfplung23 = 3.1415*0.0115^2; % surface area of sphere 23mm


spherestressarr(masksp) = 9.81*(weightplunger23+weightrodsp+weighttray + addedweightarr(masksp)/1000)/surfplung23;

%% compute stress on surface. Identify three cases: no lid vs holland lid (no rod, just added weights) vs oxford lid (+weights + rods + ring)

lidstressarr = NaN(numfiles,1);

% with oxford lid: compute effective constant stress from submersed lid plus
% weight of three rods+ring plus
% weight of edditional weights

% The oxford lid is 123mm in diameter and 8mm thick. It weighs 68.1 grams. 
% It is held by three metal rods which pass through the upper lid of the apparatus i.e the one which holds the bearing. 
% Each rod weighs 2 grams and each is 100mm long and 3.4mm diameter. 
% The upper end of the rods are attached to a perspex ring which is 100mm diameter, 
% 5mm thickness and ring width of 20mm. The ring carries the added mases has been drilled out to lighten it . 
% The ring itself weighs 41.29 grams.

%Hi Joshua, I checked through your list again. There were three small nuts used to attach the rods to the ring. 
%They each weighed 0.35gms. Hence need to add 1.05gms to the total weight of the Oxford lid. Best Tom

radlid = 0.0615; % 6.15cm radius
surflid = 3.1415*radlid^2; % surface area of lid, same for both

densxlid = 1190; % acrylic
thicklidx = 0.00485; % almost 5mm thick
weightlidxabs = surflid*thicklidx*densxlid; % NOT 68.1 grams?????????
weightlidxrel = (3.1415*thicklidx*radlid^2)*(densxlid-denswater); % correct for buoyancy
weightring = .04129; % kg, as given by Tom
weightrod = .002; % kg, times three rods. Buoyancy correction from sinking into liquid is likely negligible
weightnuts = .00035; %kg, nut for each rod, see note above.

lidstressarr = 9.81*(3*weightrod + 3*weightnuts + weightring + weightlidxrel + lidweightarr/1000)/surflid;

% the holland lid case
denshlid = 1240; % PLA
% is effectively weightless when submersed 
weightlidhrel = (3.1415*thicklidx*radlid^2)*(denshlid-denswater); % correct for buoyancy
tmp = strfind(lidonarr,'h');
lidstressarr(tmp) = 9.81*(weightlidhrel+lidweightarr(tmp)/1000)/surflid;

% without lid: provide just a constant stress = surface tension case
surftensoffs = 150; % assume surface tension offset when no lid is present. 
% Note: surface tension is not there when lid is present because particles
% do not poke through surface
tmp = strfind(lidonarr,'l');
lidstressarr(tmp) = surftensoffs;

% data interpretation settings for all figures
crtreshsp = 125; % threshold for seeing dens dep
lintreshsp = 700000000; % threshold for linear regime - does not exist yet for lid experiment

platysp = 1.6E-5; % plateau value for sphere
sloperhosp = 0.045; % slope rho dependence sphere
Dsigmsp = 3.75E-5; % prefacor sigma for sphere


%% fit individual sigma_L data to extract sigma_S(sigma_L): free & fixed

tmplidstress = unique(lidstressarr(masksp));
lengylid = length(tmplidstress);

sigmaS0arrfree = NaN(lengylid,1);
sigmaS0amparrfree = NaN(lengylid,1);
sigmaerrarrfree = NaN(lengylid,2);



for i=1:lengylid

    indytmp = find(lidstressarr == tmplidstress(i));
    tmpx = spherestressarr(indytmp);
    tmpy = slopearr(indytmp);

    %[fitobj,gof2] = fit(tmpx',tmpy,'exp1');
    opty = fitoptions('poly1','Lower', [0,-16], 'Upper', [100,-1]);
    [fitobj,gof2] = fit(tmpx',log(tmpy),'poly1',opty);
    % turns out that the weight of the largest sigm_L skews the fit, so
    % perform on lin-lin scale

    %subplot 131
    %plot(fitobj,tmpx,log(tmpy))
    %legend('off')
    %hold all
    %xlabel('\sigma_S [Pa]')
    %ylabel('log(D)')

    tmp = coeffvalues(fitobj);    
    sigmaS0amparrfree(i) = tmp(1);
    sigmaS0arrfree(i) = tmp(2);
    sigmaerrarrfree(i,1) = gof2.rmse;
    sigmaerrarrfree(i,2) = gof2.rsquare;



end



sigmaS0arrfix = NaN(lengylid,1);
sigmaS0amparrfix = NaN(lengylid,1);
sigmaerrarrfix = NaN(lengylid,2);

for i=1:lengylid

    indytmp = find(lidstressarr == tmplidstress(i));
    tmpx = spherestressarr(indytmp);
    tmpy = slopearr(indytmp);

    %[fitobj,gof2] = fit(tmpx',tmpy,'exp1');
    opty = fitoptions('poly1','Lower', [0,mean(sigmaS0arrfree)], 'Upper', [100,mean(sigmaS0arrfree)]);
    [fitobj,gof2] = fit(tmpx',log(tmpy),'poly1',opty);
    % turns out that the weight of the largest sigm_L skews the fit, so
    % perform on lin-lin scale

    tmp = coeffvalues(fitobj);    
    sigmaS0amparrfix(i) = tmp(1);
    sigmaS0arrfix(i) = tmp(2);
    sigmaerrarrfix(i,1) = gof2.rmse;
    sigmaerrarrfix(i,2) = gof2.rsquare;
end

%% save the data to a mat file for later re-use

save('SM_data.mat','-mat')

%% Figure 1 compression data 

% requires import of a separate data file, with manual readings of the
% height of the packing vs amount of weight added on packing.
compy = importdata('heights.data');
compydata = compy.data;


figure(4)

plot(compydata(:,2),compydata(:,6),'o')
hold on
errorbar(compydata(:,2),compydata(:,6),compydata(:,7),'o',"MarkerSize",1,...
    "MarkerEdgeColor","blue")

plot([0,550],[0,0.20],'-.k')


xlim([0,550])
ylim([0,0.2])
xlabel('Confinement stress [Pa]')
ylabel('\DeltaL/L [-]')
box on
ax=gca;
ax.FontSize = 14;
hold off


% save
print('figcomp', '-depsc','-r600')
print('figcomp', '-dpdf', '-bestfit')

%% Figure 1 interface position

% requires import of a separate data file, with manual readings of the
% height of the packing vs amount of hydrogels added.
tmp = importdata('interface1b.txt');

tmpdata = tmp.data;
tmpstring = tmp.textdata;
tmpstring = tmpstring{1};

figure(5)

plot(tmpdata(:,1),tmpdata(:,2),'o')
hold on
errorbar(tmpdata(:,1),tmpdata(:,2),tmpdata(:,3),'o',"MarkerSize",1,...
    "MarkerEdgeColor","blue")

plot([0,7],[165,165],'-.k')
plot([0,5.2],[0,165],':k')

xlim([0,6.5])
ylim([0,175])
xlabel('\rho [g/L]')
ylabel('L [mm]')
box on
ax=gca;
ax.FontSize = 14;
hold off

% save
print('figintf', '-depsc','-r600')
print('figintf', '-dpdf', '-bestfit')


%% Figure 2 in article: overview of data


figure(1)

% overview of 1/3rd of all data to show creep behavior
subplot 311
cmap = colormap(parula(length(unique(lidstressarr(masksp)))));

for i=1:3:length(masksp)
    
    indexer = find(unique(lidstressarr(masksp))==lidstressarr(masksp(i)));
    
    if lidweightarr(masksp(i)) < lintreshsp
       
        if lidonarr(i) == 'h'
            errorbar(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01,errarr(masksp(i),:)*0.01,'Color',cmap(indexer,:))
            scatter(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, '^k','MarkerEdgeColor',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
            
            hold on
            plot(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, '^k','Color',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
        elseif lidonarr(i) == 'l'
            scatter(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'ok','MarkerEdgeColor',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
            hold on
            plot(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'ok','Color',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
        elseif lidonarr(i) == 'x'
            scatter(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'sk','MarkerEdgeColor',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
            hold on
            plot(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'sk','Color',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
        end
       
    end
    
end

box on
xlabel('t [sec]')
ylabel('\delta [m]')
text(2500,.01,'(a)','FontSize',18)

xlim([0 2750])
ylim([0,0.1])
ax=gca;
ax.FontSize = 14;
hold off


% show same data examples, but then on log-log scale
subplot 312
cmap = colormap(parula(length(unique(lidstressarr(masksp)))));

for i=1:3:length(masksp)

    
    indexer = find(unique(lidstressarr(masksp))==lidstressarr(masksp(i)));
    
    if lidweightarr(masksp(i)) < lintreshsp
       
        if lidonarr(i) == 'h'
            scatter(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, '^k','MarkerEdgeColor',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
            hold on
            plot(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, '^k','Color',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
        elseif lidonarr(i) == 'l'
            scatter(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'ok','MarkerEdgeColor',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
            hold on
            plot(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'ok','Color',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
        elseif lidonarr(i) == 'x'
            scatter(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'sk','MarkerEdgeColor',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
            hold on
            plot(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'sk','Color',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
        end

        
        
    end
    
end

box on
xlabel('t [sec]')
ylabel('\delta [m]')
text(1900,.001,'(b)','FontSize',18)
set(gca,'Yscale','log')
set(gca,'Xscale','log')
xlim([20 3000])
ylim([0.0004,0.12])
ax=gca;
ax.FontSize = 14;
yticks([1E-3,1E-2,1E-1])

plot([1,3000],.0002*[1,3000].^.5,'-.k')
plot([1,3000],.003*[1,3000].^.5,'-.k')

hold off


% show collapse of data by dividing by D
% Note one outlier: # 17: no lid weight, leiden lid, travels only 0.330 cm.
subplot 313

for i=1:3:length(masksp) 
    
    indexer = find(unique(lidstressarr(masksp))==lidstressarr(masksp(i)));
    
    if lidstressarr(masksp(i)) < lintreshsp

        if lidonarr(i) == 'h'
            plot(timearr(masksp(i),:)*60,displarr(masksp(i),:).^2/(10000*slopearr(masksp(i))*D2si), '-^','Color',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
            hold on
        elseif lidonarr(i) == 'l'
            plot(timearr(masksp(i),:)*60,displarr(masksp(i),:).^2/(10000*slopearr(masksp(i))*D2si), '-o','Color',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
        elseif lidonarr(i) == 'x'
            plot(timearr(masksp(i),:)*60,displarr(masksp(i),:).^2/(10000*slopearr(masksp(i))*D2si), '-s','Color',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
        end
             
        hold on
        
    end
    
end
plot([10,3000],[10,3000],'-.k')
plot([10,3000],1.2*[10,3000],':k')
plot([10,3000],0.8*[10,3000],':k')
xlim([0 3000])
ylim([0 3000])

xlabel('t [sec]')
ylabel('\delta^2/D [sec]')
text(2700,500,'(c)','FontSize',18)
ax=gca;
ax.FontSize = 14;
hold off




print('fig1', '-depsc','-r600')
print('fig1', '-dpdf', '-bestfit')

%% Figure 3 then show overview figure with all data. Include lines of D function for one regime.
figure(2)


% set parameters for lines according to old exponential model
D0 = 2E-6;
sigS0 = 35;
sigL0 = 22;
prefy = 1; %irrelevant prefactor; always set to 1

% first plot the D dynamics for both cases

subplot(2,2,[1 2])

tmpspherestress = unique(spherestressarr(masksp));
cmap = colormap(jet(length(tmpspherestress)));


for fixmass  = 0:600
 
    for i=1:length(masksp)
        
        indexer = find(unique(spherestressarr(masksp))==spherestressarr(masksp(i)));


        if addedweightarr(masksp(i)) == fixmass,
            psh1 = scatter([lidstressarr(masksp(i))],[slopearr(masksp(i))*D2si],[30],[spherestressarr(masksp(i))],'Fill',...
            'MarkerEdgeColor',cmap(indexer,:),...
            'MarkerFaceColor',cmap(indexer,:),...
            'LineWidth',1.5);
        
            if lidonarr(i) == 'l'
                psh1.Marker = 'd';
            end

            if lidonarr(i) == 'h'
                psh1.Marker = '^';
            elseif lidonarr(i) == 'l'
                psh1.Marker = 'o';
            elseif lidonarr(i) == 'x'
                psh1.Marker = 'square';
            end

            hold on
            %plot(timearr(maskcyl30(i),:)*60,0.01*(timearr(maskcyl30(i),:)*slopearr(maskcyl30(i))).^0.5,'-.','Color',cmap(indexer,:))
              
            sigL = 100:1:600; % draw a line
            sigS = spherestressarr(masksp(i));
            Dfit = D0*exp(prefy*(sigS./sigS0 - (sigL)./sigL0));
            plot(sigL,Dfit,'Color',cmap(indexer,:))
        end

    end
    
end


xlim([0,600])
ylim([1e-9,1e-4])
set(gca,'Yscale','log')
xlabel('\sigma_L [Pa]')
ylabel('D [m^2/s]')
box on
hold off
c = colorbar('EastOutside','Ticks',[0,300,600,900],...
         'TickLabels',{'0','300','600', '900'});
c.Label.String = '\sigma_S [Pa]';
caxis([min(spherestressarr),max(spherestressarr)]);
ax=gca;
ax.FontSize = 14;
text(30,3E-9,'(a)','FontSize',18)
yticks([1E-9,1E-8,1E-7,1E-6,1E-5,1E-4])


%then plot the rescaled version
subplot(2,2,[3 4])

tmpspherestress = unique(spherestressarr(masksp));
cmap = colormap(jet(length(tmpspherestress)));


for fixmass  = 0:600
 
    for i=1:length(masksp)
        
        indexer = find(unique(spherestressarr(masksp))==spherestressarr(masksp(i)));


        if addedweightarr(masksp(i)) == fixmass

            sigL = lidstressarr(masksp(i));
            sigS = spherestressarr(masksp(i));

            Dcorr = D0*exp(prefy*(sigS./sigS0)); 


            psh1 = scatter([lidstressarr(masksp(i))],[slopearr(masksp(i))*D2si/Dcorr],[30],[spherestressarr(masksp(i))],'Fill',...
            'MarkerEdgeColor',cmap(indexer,:),...
            'MarkerFaceColor',cmap(indexer,:),...
            'LineWidth',1.5);
        
            if lidonarr(i) == 'h'
                psh1.Marker = '^';
            elseif lidonarr(i) == 'l'
                psh1.Marker = 'o';
            elseif lidonarr(i) == 'x'
                psh1.Marker = 'square';
            end

            hold on
             
        end

    end
    
end

sigL = 100:1:600; % draw a line
Dfit = exp(-prefy*((sigL)./sigL0));
plot(sigL,Dfit,'-.k')

ylim([5e-13,1e-2])
set(gca,'Yscale','log')
xlabel('\sigma_L [Pa]')
ylabel('D/D_0exp((\sigma_S/\sigma_0)) [-]')
box on
text(30,7E-12,'(b)','FontSize',18)

hold off

c = colorbar('EastOutside','Ticks',[0,300,600,900],...
         'TickLabels',{'0','300','600', '900'});
c.Label.String = '\sigma_S [Pa]';
caxis([min(spherestressarr),max(spherestressarr)]);
ax=gca;
ax.FontSize = 14;
s=text(100,1E-10,horzcat('\sigma_{L0} = ',num2str(sigL0),' Pa'));
s.FontSize = 20;

yticks([1E-12,1E-10,1E-8,1E-6,1E-4,1E-2])

% save
print('figsurfstress', '-depsc','-r600')
print('figsurfstress', '-dpdf', '-bestfit')



%% Figure 4. Show how D depends on added weights on sphere, for different lid stresses. 


figure(9)

% first, show the D versus sigma_S to highlight new sigma_L dependent
% behavior

subplot(2,2,[1 2]);

tmplidstress = unique(lidstressarr(masksp));
cmap = colormap(jet(length(tmplidstress)));


for i=1:length(masksp)
        
    indexer = find(unique(lidstressarr(masksp))==lidstressarr(masksp(i)));
    
    if lidweightarr(masksp(i)) < lintreshsp, % if statement has no influence; always occurs. Leftover from old code.
        psh1 = scatter([spherestressarr(masksp(i))],[slopearr(masksp(i))*D2si],[90],'Fill',...
            'MarkerEdgeColor',cmap(indexer,:),...
            'MarkerFaceColor',cmap(indexer,:),...
            'LineWidth',1.5');
            if lidonarr(i) == 'h'
                psh1.Marker = '^';
            elseif lidonarr(i) == 'l'
                psh1.Marker = 'o';
            elseif lidonarr(i) == 'x'
                psh1.Marker = 'square';
            end
            hold on
    end
    
    sigL = lidstressarr(masksp(i));
    if sigL >= 0 % above threshold, lid stress dep exists
        %use values from Fig 3
        sigS = 50:1:1000;
        Dfit = D0*exp(sigmaS0arrfree(indexer))*exp((sigS).*sigmaS0amparrfree(indexer));
        plot(sigS,Dfit,'Color',cmap(indexer,:)); 
    else % below threshold, no lid stress dependence
    	% does not occur
    end

end


xlim([0,1000])
ylim([1e-9,1E-4])
set(gca,'Yscale','log')
xlabel('\sigma_S [Pa]')
ylabel('D [m^2/s]')
box on
c = colorbar('EastOutside', 'Ticks',[0,100,200,300,400,500],...
         'TickLabels',{'0','100','200','300','400','500'});
clim([min(tmplidstress),max(tmplidstress)]);
ax=gca;
ax.FontSize = 14;
c.Label.String = '\sigma_L [Pa]';
text(20,3E-5,'(a)','FontSize',18)


% plot free & fixed fit parameters

% first plot 1/A
subplot 223

loglog(tmplidstress,1./sigmaS0amparrfree,'-ok')
hold on
loglog(tmplidstress,1./sigmaS0amparrfix,'-sb')

plot(1:500, 18+0.09*(1:500),'-.r')

xlim([10,600])

xlabel('\sigma_L [Pa]')
ylabel('1/A [Pa]')
ax=gca;
ax.FontSize = 14;
text(12,67,'(b)','FontSize',18)

% then plot B
subplot 224
plot(tmplidstress,(sigmaS0arrfree),'-ok')
hold on
plot([0,600],[mean(sigmaS0arrfree),mean(sigmaS0arrfree)],'-b')

xlabel('\sigma_L [Pa]')
ylabel('B [-]')
ax=gca;
ax.FontSize = 14;
text(500,-8.7,'(c)','FontSize',18)

% save
print('figfit', '-depsc','-r600')
print('figfit', '-dpdf', '-bestfit')


%% Fig 5 Show collapse in separate figure

tmplidstress = unique(lidstressarr(masksp));
cmap = colormap(jet(length(tmplidstress)));

figure(67)

% first the collapse with non-constant B
subplot(2,1,1);

for i=1:length(masksp)

    indexer = find(unique(lidstressarr(masksp))==lidstressarr(masksp(i)));

    sigL = lidstressarr(masksp(i));
    if sigL >= 0 % above threshold, lid stress dep exists
        sigS = 50:1:1000;
        Dfit = D0*exp(sigS./sigS0)*exp(-sigL./sigL0); 

        rescalecoeff = 1./sigmaS0amparrfree(indexer);


        psh1 = scatter([spherestressarr(masksp(i))]./rescalecoeff,[slopearr(masksp(i))*D2si]/(exp(sigmaS0arrfree(indexer))),[90],'Fill',...
        'MarkerEdgeColor',cmap(indexer,:),...
        'MarkerFaceColor',cmap(indexer,:),...
        'LineWidth',1.5');
            if lidonarr(i) == 'h'
                psh1.Marker = '^';
            elseif lidonarr(i) == 'l'
                psh1.Marker = 'o';
            elseif lidonarr(i) == 'x'
                psh1.Marker = 'square';
            end
        hold on
    end

end

sigS = 0:.1:20; % draw a line
Dfit = exp(sigS);
plot(sigS,D0*Dfit,'-.k');
xlim([0,18])%([0,17])
ylim([.9E-5,100])%([1,2.5E7])
set(gca,'Yscale','log')
ylabel('D/exp(B(\sigma_L)) [-]')
xlabel('\sigma_S/A(\sigma_L) [-]')
clim([min(tmplidstress),max(tmplidstress)]);
ax=gca;
ax.FontSize = 14;
box on
text(1,10,'(a)','FontSize',18)

hold off

% subplot with different rescaling choice: const B   %%%%%%%

subplot(2,1,2);

for i=1:length(masksp)

    indexer = find(unique(lidstressarr(masksp))==lidstressarr(masksp(i)));

    sigL = lidstressarr(masksp(i));
    if sigL >= 0 % above threshold, lid stress dep exists
        sigS = 50:1:1000;
        Dfit = D0*exp(sigS./sigS0)*exp(-sigL./sigL0); 

        %rescaling with fixed B
        rescalecoeff = 1./sigmaS0amparrfix(indexer);

        psh2 = scatter([spherestressarr(masksp(i))]./rescalecoeff,[slopearr(masksp(i))*D2si]/(exp(-11.2)),[90],'Fill',...
        'MarkerEdgeColor',cmap(indexer,:),...
        'MarkerFaceColor',cmap(indexer,:),...
        'LineWidth',1.5');
            if lidonarr(i) == 'h'
                psh2.Marker = '^';
            elseif lidonarr(i) == 'l'
                psh2.Marker = 'o';
            elseif lidonarr(i) == 'x'
                psh2.Marker = 'square';
            end
        hold on
    end

end

sigS = 0:.1:20; % show exponential behavior
Dfit = exp(sigS);
plot(sigS,D0*Dfit,'-.k');
xlim([0,18])
ylim([.9E-5,10])
set(gca,'Yscale','log')

box on

c = colorbar('SouthOutside', 'Ticks',[0,100,200,300,400,500],...
         'TickLabels',{'0','100','200','300','400','500'});
clim([min(tmplidstress),max(tmplidstress)]);
ax=gca;
ax.FontSize = 14;
c.Label.String = '\sigma_L [Pa]';
ylabel('D/exp(B) [-]')

xlabel('\sigma_S/A(\sigma_L) [-]')
box on

text(1,1,'(b)','FontSize',18)


hold off


% save
print('figfincoll', '-depsc','-r600')
print('figfincoll', '-dpdf', '-bestfit')

%% Verify nonlinearity with lid 
% (a) Time series for cases with lid plus (b) collapse 
% indexchoice 13 is equivalent to an exponent of 0.5

indexchoice = 13;

lidchoice = 1;
alphy = 1/fitexpyarr(indexchoice);

masklid = 1:numfiles;

figure(10)
subplot(2,2,[1 2])
cmap = colormap(parula(length(unique(spherestressarr(masklid)))));

for i=1:1:length(masklid)
    
    indexer = find(unique(spherestressarr(masklid))==spherestressarr(masklid(i)));
    

    scatter(timearr(masklid(i),:)*60,displarr(masklid(i),:), 'k','MarkerEdgeColor',cmap(indexer,:))
    hold on
    plot(timearr(masklid(i),:)*60,displarr(masklid(i),:), 'k','Color',cmap(indexer,:))

   
end

box on
xlabel('t [sec]')
ylabel('\delta [cm]')
xlim([0 2500])
ylim([0, 10])
ax=gca;
ax.FontSize = 14;
hold off
c = colorbar('EastOutside');
c.Label.String = '\sigma_S [Pa]';

clim([0,250]);
text(3550,90,'(a)','FontSize',18)


subplot 223

for i=1:1:length(masklid) 
    
    indexer = find(unique(spherestressarr(masklid))==spherestressarr(masklid(i)));
  

    scatter(timearr(masklid(i),:)*60,(displarr(masklid(i),:)./slopearr_nonlin(masklid(i),indexchoice)).^alphy, 'k','MarkerEdgeColor',cmap(indexer,:))
    hold on
   
              
end

plot([10,10000],0.012*[10,10000],'-.k')

box on
xlabel('t [sec]')
ylabel(horzcat('(\delta\eta_{eff})^{',num2str(alphy,3),'} [A.U]'))
ax=gca;
ax.FontSize = 14;
hold off

clim([0,250])
xlim([20 3000])
ylim([1E-2, 2E2])
set(gca,'Yscale','log')
set(gca,'Xscale','log')

%save
print('checknl', '-dpdf', '-bestfit')




