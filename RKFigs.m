%% Plot RK figs
gail.InitializeWorkspaceDisplay

differ = @(t,x) (reshape(t,[size(t,1),1,size(t,2)]) - ...
    reshape(x,[1,size(x,1),size(x,2)]));

normsq = @(t,x,c) squeeze(sum((differ(t,x) .* reshape(c,[1,1,size(c,2)])).^2,3)); 

normof = @(t,x,c) sqrt(normsq(t,x,c));

sqexpker = @(t,x,c) exp(-normsq(t,x,c));


ctrdiscker = @(t,x,c) squeeze(prod(1 + (abs(reshape(t,[size(t,1),1,size(t,2)])-1/2) +...
    abs(reshape(x,[1,size(x,1),size(x,2)])-1/2) - ...
    differ(t,x)) .* reshape(c,[1,1,size(c,2)]),3));

maternkerthreehalfs = @(t,x,c) (1 + normof(t,x,c)).*exp(-normof(t,x,c));

deltaker = @(t,x,c) 1 + c(1).*(differ(t,x)==0);


kernelEx = {sqexpker,ctrdiscker,maternkerthreehalfs};
kernelName = {'sqexpker','ctrdiscker','maternkerthreehalfs'};
cVal = [0.5 2 8]';
c2Val = [8 2];
tc = 0.5;
tc2 = [0.5 0.5];
xplot = (0:0.002:1)';
xx = (0:0.01:1)';
[xxplot2, yyplot2] = meshgrid(xx);
xplot2 = [xxplot2(:) yyplot2(:)];
kplot(length(xplot),length(cVal)) = 0;
kplot2(size(xplot2,1),1) = 0;
legend
for jj = 1:length(kernelEx)
    kernel = kernelEx{jj};
    figure(2*jj-1);
    hold on
    for ii = 1:length(cVal)
        kplot(:,ii) = kernel(xplot,tc,cVal(ii));
    end
    kplot(isnan(kplot)) = 1;
    h = plot(xplot,kplot);
    xlabel('\(x\)')
    if jj < 3
        ylabel(['\(K(x,' num2str(tc) ')\)'])
    else
        ylabel(['\(K_{3/2}(x,' num2str(tc) ')\)'])
    end
    legendCell = cellstr(num2str(cVal, '\\(\\gamma = %-g\\)'));
    legend(h,legendCell,'Box','off','location','best')
    print('-depsc',['RK-' kernelName{jj} '.eps'])

    figure(2*jj);
    kplot2(:) = kernel(xplot2,tc2,c2Val);
    kk = reshape(kplot2,size(xx,1),size(xx,1));
    h = surf(xxplot2,yyplot2,kk,'FaceColor','interp','EdgeColor','interp');
    xlabel('\(x_1\)')
    ylabel('\(x_2\)')
    zlabel(['\(K\bigl((x_1,x_2),(' num2str(tc2(1)) ',' num2str(tc2(1)) ') \bigr)\)'])
    title(['\(\gamma_1 = ' num2str(c2Val(1)) ', \gamma_2 = ' num2str(c2Val(2)) '\)' ])
    print('-depsc',['RK2-' kernelName{jj} '.eps'])

end

%% For delta kernel
figure
kplot2 = 1 + (xplot2(:,1) == xplot2(:,2));
kk = reshape(kplot2,size(xx,1),size(xx,1));
h = surf(xxplot2,yyplot2,kk,'FaceColor','interp','EdgeColor','interp');
xlabel('\(x_1\)')
ylabel('\(x_2\)')
zlabel(['\(K(x_1,x_2)\)'])
title('\(\gamma = 1\)')
print('-depsc','RK-deltaker.eps')

