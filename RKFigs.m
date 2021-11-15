%% Plot RK figs
gail.InitializeWorkspaceDisplay

absdif = @(t,x) (reshape(t,[size(t,1),1,size(t,2)]) - ...
    reshape(x,[1,size(x,1),size(x,2)]));

sqexpker = @(t,x,c) squeeze(exp(-sum((absdif(t,x) .* reshape(c,[1,1,size(c,2)])).^2,3)));


ctrdiscker = @(t,x,c) squeeze(prod(1 + (abs(reshape(t,[size(t,1),1,size(t,2)])-1/2) +...
    abs(reshape(x,[1,size(x,1),size(x,2)])-1/2) - ...
    absdif(t,x)) .* reshape(c,[1,1,size(c,2)]),3));


kernelEx = {sqexpker,ctrdiscker};
kernelName = {'sqexpker','ctrdiscker'};
cVal = [0.5 2 4]';
c2Val = [4 2];
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
    h = plot(xplot,kplot);
    xlabel('\(x\)')
    ylabel(['\(K(x,' num2str(tc) ')\)'])
    legendCell = cellstr(num2str(cVal, 'c = %-g'));
    legend(h,legendCell,'Box','off','location','best')
    print('-depsc',['RK-' kernelName{jj} '.eps'])

    figure(2*jj);
    kplot2(:) = kernel(xplot2,tc2,c2Val);
    kk = reshape(kplot2,size(xx,1),size(xx,1));
    h = surf(xxplot2,yyplot2,kk,'FaceColor','interp','EdgeColor','interp');
    xlabel('\(x_1\)')
    ylabel('\(x_1\)')
    zlabel(['\(K\bigl((x_1,x_2),(' num2str(tc2(1)) ',' num2str(tc2(1)) ') \bigr)\)'])
    print('-depsc',['RK2-' kernelName{jj} '.eps'])
    title(['\(c_1 = ' num2str(c2Val(1)) ', c_2 = ' num2str(c2Val(2)) '\)' ])

end
