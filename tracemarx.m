
% Visualising topo for active margin model
clear all;
clf;
% Approximate number of sampled markers and *.prn files
msample=22000;
fsample=1;
% Create marker arrays
mmm=zeros(msample,1); % marker global index
mty=zeros(msample,fsample); % type
mti=zeros(msample,fsample); % time, yr
mxx=zeros(msample,fsample); % x, m
myy=zeros(msample,fsample); % y, m
mtk=zeros(msample,fsample); % T, K
mpb=zeros(msample,fsample); % P, bar
% File name
nname='cdf78_';
movieaname = 'test1.avi';
moviebname = 'test2.avi';
moviecname = 'test3.avi';
savename = 'test.txt';
% Coordinates of the sampling area, m
xmin=500000;
xmax=1260000;
ymin=17500;
ymax=28500;
% Counters
mfind=0; % Marker counter
fcur=1; % File counter
a=1; % axis counter
q=1; %% movie frame counter
axisconstant = 7; % Speed of axis march
    for nnf=2:1:4%1:1:fsample
        % Data filename
        fileprn=[nname,num2str(nnf*10),'.prn'];
        % Opening data file
        fdata=fopen(fileprn,'rb');
        % Read sizes of variables
        A=fread(fdata,4,'uchar');
        % Read model parameters
        % Grid resolution
        xnumx=fread(fdata,1,'int64');
        ynumy=fread(fdata,1,'int64');
        % Markers per cell
        mnumx=fread(fdata,1,'int64');
        mnumy=fread(fdata,1,'int64');
        % Number of markers
        marknum=fread(fdata,1,'int64');
        % Model sizes
        xsize=fread(fdata,1,'float64');
        ysize=fread(fdata,1,'float64');
        % Pressure value
        pinit=fread(fdata,5,'float64');
        % Gravity
        GXKOEF=fread(fdata,1,'float64');
        GYKOEF=fread(fdata,1,'float64');
        % Number of rocks
        rocknum=fread(fdata,1,'int32');
        % Number of Boundary conditions
        bondnum=fread(fdata,1,'int64');
        % Stage,time
        n1=fread(fdata,1,'int32');
        timesum=fread(fdata,1,'float64');
        % Skip rock properties
        curpos0=4+2*4+16*8+rocknum*(8*24+4);
        fseek(fdata,curpos0,'bof');
        
        % Read nodes information
        for i=1:1:xnumx
            for j=1:1:ynumy
                vbuf=fread(fdata,3,'float32');
                pr(j,i)=vbuf(1);
                vx(j,i)=vbuf(2);
                vy(j,i)=vbuf(3);
                vbuf1=fread(fdata,3,'int64');
                vbuf2=fread(fdata,16,'float32');
                exx(j,i)=vbuf2(1);
                eyy(j,i)=vbuf2(2);
                exy(j,i)=vbuf2(3);
                sxx(j,i)=vbuf2(4);
                syy(j,i)=vbuf2(5);
                sxy(j,i)=vbuf2(6);
                ro(j,i)=vbuf2(7);
                nu(j,i)=vbuf2(8);
                nd(j,i)=vbuf2(9);
                mu(j,i)=vbuf2(10);
                ep(j,i)=vbuf2(11);
                et(j,i)=vbuf2(12);
                pr0(j,i)=vbuf2(13);
                prb(j,i)=vbuf2(14);
                dv(j,i)=vbuf2(15);
                tk(j,i)=vbuf2(16);
                vbuf3=fread(fdata,1,'int64');
                vbuf4=fread(fdata,3,'float32');
                cp(j,i)=vbuf4(1);
                kt(j,i)=vbuf4(2);
                ht(j,i)=vbuf4(3);
            end
        end
        
        % Skip all nodes
        curpos2=curpos0+(4*22+8*4)*xnumx*ynumy;
        fseek(fdata,curpos2,'bof');
        % Read Gridline positions
        gx=fread(fdata,xnumx,'float32');
        gy=fread(fdata,ynumy,'float32');
        eii=ones(ynumy,xnumx)*1e-16;
        sii=ones(ynumy,xnumx)*1e+4;
        for i=1:1:xnumx-2
            for j=1:1:ynumy-2
                eii(j+1,i+1)=(exy(j+1,i+1)^2+((exx(j+1,i+1)+exx(j+2,i+1)+exx(j+1,i+2)+exx(j+2,i+2))/4)^2)^0.5;
                sii(j+1,i+1)=(sxy(j+1,i+1)^2+((sxx(j+1,i+1)+sxx(j+2,i+1)+sxx(j+1,i+2)+sxx(j+2,i+2))/4)^2)^0.5;
            end
        end
        
        % Skip all BC
        curpos3=curpos2+(xnumx+ynumy)*4+(4*4+8*3)*(bondnum-1);
        
        % Pressure points coordinates
        px=gx;px(2:xnumx)=(gx(1:xnumx-1)+gx(2:xnumx))/2;
        py=gy;py(2:ynumy)=(gy(1:ynumy-1)+gy(2:ynumy))/2;
        % Find Markers
        if(mfind==0)
            fseek(fdata,curpos3,'bof');
            for m=1:1:marknum
                mbuf=fread(fdata,9,'float32');
                mt=fread(fdata,1,'uchar');
                mx=mbuf(1);
                my=mbuf(2);
                mk=mbuf(3);
                % Save markers from the interest area
                if(mx>=xmin && mx<=xmax && my>=ymin && my<=ymax && mt>1)
                    mfind=mfind+1;
                    mmm(mfind)=m; % global index
                    mty(mfind,nnf)=mt; % type
                    mti(mfind,nnf)=timesum; % time, yr
                    mxx(mfind,nnf)=mx; % x, m
                    myy(mfind,nnf)=my; % y, m
                    mtk(mfind,nnf)=mk; % T, K
                    % Interpolate pressure
                    % Find indexes for the upper-left pressure node by bisection
                    jmin=2;jmax=ynumy;
                    while(jmax-jmin>1)
                        j=fix((jmax+jmin)/2);
                        if(py(j)>my)
                            jmax=j;
                        else
                            jmin=j;
                        end
                    end
                    j=jmin;
                    if(j>ynumy-1)
                        j=ynumy-1;
                    end
                    imin=2;imax=xnumx;
                    while(imax-imin>1)
                        i=fix((imax+imin)/2);
                        if(px(i)>mx)
                            imax=i;
                        else
                            imin=i;
                        end
                    end
                    i=imin;
                    if(i>xnumx-1)
                        i=xnumx-1;
                    end
                    % Bilinear interpolation
                    dxm=(mx-px(i))/(px(i+1)-px(i));
                    dym=(my-py(j))/(py(j+1)-py(j));
                    % P, bar
                    mpb(mfind,nnf)=1e-5*((1-dxm)*(1-dym)*pr(j,i)+(dxm)*(1-dym)*pr(j,i+1)+...
                        (1-dxm)*(dym)*pr(j+1,i)+(dxm)*(dym)*pr(j+1,i+1));
                end
            end
            if(mfind==0)
                break;
            end
            
            
            % Update PT-paths of selected markers
        else
            for mi=1:1:mfind
                m=mmm(mi);
                curpos4=curpos3+(m-1)*(9*4+1);
                fseek(fdata,curpos4,'bof');
                mbuf=fread(fdata,9,'float32');
                mt=fread(fdata,1,'uchar');
                mx=mbuf(1);
                my=mbuf(2);
                mk=mbuf(3);
                % Save marker data
                mty(mi,nnf)=mt; % type
                mti(mi,nnf)=timesum; % time, yr
                mxx(mi,nnf)=mx; % x, m
                myy(mi,nnf)=my; % y, m
                mtk(mi,nnf)=mk; % T, K
                
                % Interpolate pressure
                % Find indexes for the upper-left pressure node by bisection
                jmin=2;jmax=ynumy;
                while(jmax-jmin>1)
                    j=fix((jmax+jmin)/2);
                    if(py(j)>my)
                        jmax=j;
                    else
                        jmin=j;
                    end
                end
                j=jmin;
                if(j>ynumy-1)
                    j=ynumy-1;
                end
                imin=2;imax=xnumx;
                while(imax-imin>1)
                    i=fix((imax+imin)/2);
                    if(px(i)>mx)
                        imax=i;
                    else
                        imin=i;
                    end
                end
                i=imin;
                if(i>xnumx-1)
                    i=xnumx-1;
                end
                % Bilinear interpolation
                dxm=(mx-px(i))/(px(i+1)-px(i));
                dym=(my-py(j))/(py(j+1)-py(j));
                % P, bar
                mpb(mi,nnf)=1e-5*((1-dxm)*(1-dym)*pr(j,i)+(dxm)*(1-dym)*pr(j,i+1)+...
                    (1-dxm)*(dym)*pr(j+1,i)+(dxm)*(dym)*pr(j+1,i+1));
                if(mtk(mi,nnf)<0 || mxx(mi,nnf)<102000 || myy(mi,nnf)>200000 || mty(mi,nnf)==0 || mty(mi,nnf)==1 || mty(mi,nnf)==5 || mty(mi,nnf)==6 || mty(mi,nnf)==8 || mty(mi,nnf)==9 || mty(mi,nnf)>=10)
                    mxx(mi,nnf) = NaN;
                    myy(mi,nnf) = NaN;
                    mtk(mi,nnf) = NaN;
                    mpb(mi,nnf) = NaN;
                    mty(mi,nnf) = NaN;
                end
            end
            
        end
        fclose(fdata);
        if(mfind==0)
            break;
        end
        
        % Draw
        figure(1);clf;
        colormap('jet');
        p2 = pcolor(gx/1000,gy/1000,log10(nd));
        shading interp;
        axis ij image;
        colorbar;
        title(['Log Viscosity (Pa*s) | Time = ',num2str(timesum*1e-6, '%.2f'), ' Ma']);
        hold on
        [c,h]=contour(gx/1000,gy/1000,(tk-273.15)/300+17,(100:200:2500)/300+17,'w');
        xlabel(gca, 'Distance (km)','FontSize',13, 'FontName', 'Arial','fontweight','bold');
        ylabel(gca, 'Depth (km)','FontSize', 13, 'FontName', 'Arial','fontweight','bold');
        set(gca,'FontSize', 12, 'FontName', 'Arial','fontweight','bold','LineWidth', 2);
        set(gcf,'color','w');
        
        if (nnf>0)
            hold on
            sz = 40;
            for mi=1:1:mfind
                if (mty(mi,(nnf))==0)
                    c = [1 1 1];
                elseif (mty(mi,(nnf))==1)
                    c = [0.50588 0.99608 0.78824];
                elseif (mty(mi,(nnf))==2)
                    c = [1 0.50196 0];
                elseif (mty(mi,(nnf))==3)
                    c = [0.68235 0.34118 0];
                elseif (mty(mi,(nnf))==4)
                    c = [1 0.50196 0];
                elseif (mty(mi,(nnf))==5)
                    c = [0.75294 0.75294 0.75294];
                elseif (mty(mi,(nnf))==6)
                    c = [0.50196 0.50196 0.50196];
                elseif (mty(mi,(nnf))==7)
                    c = [0 0.50196 0];
                elseif (mty(mi,(nnf))==8)
                    c = [0 0.84314 0];
                elseif (mty(mi,(nnf))==9)
                    c = [0 0 0.71765];
                elseif (mty(mi,(nnf))==10)
                    c = [0.3 0.3 0.9];
                elseif (mty(mi,(nnf))==11)
                    c = [0.14118 0.72157 0.99216];
                elseif (mty(mi,(nnf))==12)
                    c = [0 0.50196 1];
                elseif (mty(mi,(nnf))==13)
                    c = [0 0 0.4902];
                elseif (mty(mi,(nnf))==14)
                    c = [0.4 0 0];
                elseif (mty(mi,(nnf))==15)
                    c = [0.8549 0.59608 0.36078];
                else
                    c = [0 0 0];
                end
                p3 = plot(mxx(mi,nnf)/1000,myy(mi,nnf)/1000,'Marker','s','MarkerSize',10,'MarkerFaceColor', c,'MarkerEdgeColor', 'k');%,'Color',c);
                p3.Color(4) = 0;
            end
        end
        
        mova(q) = getframe(gcf);
        
        figure(2);clf;
        colormap('Jet');
        hold on
        for mi=1:1:mfind
            if (mty(mi,nnf)==0)
                c = [1 1 1];
            elseif (mty(mi,nnf)==1)
                c = [0.50588 0.99608 0.78824];
            elseif (mty(mi,nnf)==2)
                c = [1 0.50196 0];
            elseif (mty(mi,nnf)==3)
                c = [0.68235 0.34118 0];
            elseif (mty(mi,nnf)==4)
                c = [1 0.50196 0];
            elseif (mty(mi,nnf)==5)
                c = [0.75294 0.75294 0.75294];
            elseif (mty(mi,nnf)==6)
                c = [0.50196 0.50196 0.50196];
            elseif (mty(mi,nnf)==7)
                c = [0 0.50196 0];
            elseif (mty(mi,nnf)==8)
                c = [0 0.84314 0];
            elseif (mty(mi,nnf)==9)
                c = [0 0 0.71765];
            elseif (mty(mi,nnf)==10)
                c = [0.3 0.3 0.9];
            elseif (mty(mi,nnf)==11)
                c = [0.14118 0.72157 0.99216];
            elseif (mty(mi,nnf)==12)
                c = [0 0.50196 1];
            elseif (mty(mi,nnf)==13)
                c = [0 0 0.4902];
            elseif (mty(mi,nnf)==14)
                c = [0.4 0 0];
            elseif (mty(mi,nnf)==15)
                c = [0.8549 0.59608 0.36078];
            else
                c = [0 0 0];
            end
            
            plot1= plot(mtk(mi,1:nnf)-273.15,mpb(mi,1:nnf)/1000,'o-', 'Color', c,'MarkerFaceColor', c, 'MarkerSize', 8);
            plot1.Color(4) = 0.1;
            box on
            text(825, 2.5, ['t = ',num2str(timesum*1e-6, '%.2f'), ' Ma'],'FontSize', 16, 'FontName', 'Arial','fontweight','normal');
            axis([0 1000 0 40])
            set(gca, 'DataAspectRatio', [20 1 1])
            xlabel(gca, 'Temperature (C)','FontSize', 13, 'FontName', 'Arial','fontweight','bold');
            ylabel(gca, 'Pressure (kbar)','FontSize', 13, 'FontName', 'Arial','fontweight','bold');
            set(gca,'FontSize', 12, 'FontName', 'Arial','fontweight','bold','LineWidth', 2);
            set(gcf,'color','w');
        end
        %title('P-T Trajectory');
        %export_fig(jpgname2,'-png','-r300');%, '-c[NaN,450,NaN,550]');
        movb(q) = getframe(gcf);
        
        figure(3);clf;
        colormap('Jet');
        hold on
        for mi=1:1:mfind
            if (mty(mi,nnf)==0)
                c = [1 1 1];
            elseif (mty(mi,nnf)==1)
                c = [0.50588 0.99608 0.78824];
            elseif (mty(mi,nnf)==2)
                c = [1 0.50196 0];
            elseif (mty(mi,nnf)==3)
                c = [0.68235 0.34118 0];
            elseif (mty(mi,nnf)==4)
                c = [1 0.50196 0];
            elseif (mty(mi,nnf)==5)
                c = [0.75294 0.75294 0.75294];
            elseif (mty(mi,nnf)==6)
                c = [0.50196 0.50196 0.50196];
            elseif (mty(mi,nnf)==7)
                c = [0 0.50196 0];
            elseif (mty(mi,nnf)==8)
                c = [0 0.84314 0];
            elseif (mty(mi,nnf)==9)
                c = [0 0 0.71765];
            elseif (mty(mi,nnf)==10)
                c = [0.3 0.3 0.9];
            elseif (mty(mi,nnf)==11)
                c = [0.14118 0.72157 0.99216];
            elseif (mty(mi,nnf)==12)
                c = [0 0.50196 1];
            elseif (mty(mi,nnf)==13)
                c = [0 0 0.4902];
            elseif (mty(mi,nnf)==14)
                c = [0.4 0 0];
            elseif (mty(mi,nnf)==15)
                c = [0.8549 0.59608 0.36078];
            else
                c = [0 0 0];
            end
            
            if(nnf>1)
                plot3= plot(mtk(mi,1:nnf)-273.15,mpb(mi,1:nnf)/1000,'-', 'Color', c,'MarkerFaceColor', 'none');
                plot3.Color(4) = 0.5;
                box on
            end
            plot2= plot(mtk(mi,nnf)-273.15,mpb(mi,nnf)/1000,'o', 'Color', c,'MarkerFaceColor', c, 'MarkerSize', 8);
            box on
            text(825, 2.5, ['t = ',num2str(timesum*1e-6, '%.2f'), ' Ma'],'FontSize', 16, 'FontName', 'Arial','fontweight','normal');
            axis([0 1000 0 40])
            set(gca, 'DataAspectRatio', [20 1 1])
            xlabel(gca, 'Temperature (C)','FontSize', 13, 'FontName', 'Arial','fontweight','bold');
            ylabel(gca, 'Pressure (kbar)','FontSize', 13, 'FontName', 'Arial','fontweight','bold');
            set(gca,'FontSize', 12, 'FontName', 'Arial','fontweight','bold','LineWidth', 2);
            set(gcf,'color','w');
        end
        %title('P-T Trajectory');
        %export_fig(jpgname2,'-png','-r300');%, '-c[NaN,450,NaN,550]');
        movc(q) = getframe(gcf);
        
        pause(0.1);
        a=a+1;
        q=q+1;
    end

figure(4);clf;
set(gca, 'Visible', 'off','position',[0 0 1 1],'units','normalized');
movie(mova,2,5)
movie(movb,2,5)
movie(movc,2,5)

va = VideoWriter([nname,movieaname]);
va.FrameRate=3;
va.Quality=100;
open(va);
writeVideo(va,mova);
close(va);

vb = VideoWriter([nname,moviebname]);
vb.FrameRate=3;
vb.Quality=100;
open(vb);
writeVideo(vb,movb);
close(vb);

vc = VideoWriter([nname,moviecname]);
vc.FrameRate=3;
vc.Quality=100;
open(vc);
writeVideo(vc,movc);
close(vc);

save([nname,savename],'mty','mti','mxx','myy','mtk','mpb','-ascii','-tabs')