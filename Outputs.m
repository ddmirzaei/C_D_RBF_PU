function Outputs(N,InfErr,Orders,SetupTime,SolveTime,nz_percent,...
                 MethodType,WeightType,BoundType,RBFtype,pVec)
% This function prints the outputs and plots the errors 
%   and computational costs in terms of N
% 
if length(pVec) == 3
    if strcmp(MethodType,'pu')
        StrMeth = 'Compact D-RBF-PU';
    else
        StrMeth = 'RBF-HFD';
    end
    if strcmp(WeightType,'Smooth')
        StrWeight = 'Smooth Weight';
    else
        StrWeight = 'ConstGen Weight';
    end
    
    disp('CPUtime ='); disp(SetupTime+SolveTime)
    disp('InfErr ='); disp(InfErr)
    disp('Orders ='); disp(Orders)
    
    figure;
    loglog(sqrt(N),InfErr(:,1),'rs-','MarkerSize',6,'MarkerFaceColor','r','LineWidth',1.3)
    hold on
    loglog(sqrt(N),InfErr(:,2),'bo-','MarkerSize',6,'MarkerFaceColor','b','LineWidth',1.3)
    hold on
    loglog(sqrt(N),InfErr(:,3),'m>-','MarkerSize',6,'MarkerFaceColor','m','LineWidth',1.3)
    xlabel('$\sqrt{N}$','interpreter','latex')
    ylabel('$\|e\|_\infty/\|u\|_\infty$','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    if strcmp(RBFtype,'p')
        leg = legend('PHS5+P2','PHS7+P3', 'PHS9+P4','Location','northeast');
    else % 'tp'
        leg = legend('PHS4+P2','PHS6+P3', 'PHS8+P4','Location','northeast');
    end
    set(leg,'Interpreter','latex');
    set(gcf, 'Position', [300 300 350 400])
    title(sprintf('%s, %s',StrMeth,StrWeight),'Interpreter','latex');
    xlim([10 2*10^2])
    xticks([10 20 50 100 200])
    text(35,0.01,num2str(sprintf('%1.2f\n',Orders(:))),'Interpreter','latex');
    hold off
    
else
    
    if strcmp(BoundType,'Dirichlet')
        StrBound = 'Dirichlet BC';
    else
        StrBound = 'Neumann BC';
    end
    
    % plot sparsity
    figure
    loglog(sqrt(N(3:end)),nz_percent(3:end,1),'rs-','MarkerSize',6,'MarkerFaceColor','r','LineWidth',1.3)
    hold on
    loglog(sqrt(N(3:end)),nz_percent(3:end,2),'go-','MarkerSize',6,'MarkerFaceColor','g','LineWidth',1.3)
    hold on
    loglog(sqrt(N(3:end)),nz_percent(3:end,3),'k>-','MarkerSize',6,'MarkerFaceColor','k','LineWidth',1.3)
    
    ylabel('Non zeros ($\%$)', 'interpreter', 'latex')
    xlabel('$\sqrt N$', 'interpreter', 'latex')
    set(gcf, 'Position', [300 300 400 275])
    set(gca,'TickLabelInterpreter','latex')
    title('$\rho=4h$','Interpreter','latex')
    leg = legend('D-RBF-PU, Smooth','D-RBF-PU, ConstGen', 'RBF-HFD','Location','northeast');
    set(leg,'Interpreter','latex');
    xlim([15 2*10^2])
    ylim([0 35])
    xticks([20 50 100 200])
    yticks([0.5 1 3 10 20])
    
    % plot errors
    figure;
    loglog(sqrt(N),InfErr(:,1),'rs-','MarkerSize',6,'MarkerFaceColor','r','LineWidth',1.3)
    hold on
    loglog(sqrt(N),InfErr(:,2),'go-','MarkerSize',6,'MarkerFaceColor','g','LineWidth',1.3)
    hold on
    loglog(sqrt(N),InfErr(:,3),'k>-','MarkerSize',6,'MarkerFaceColor','k','LineWidth',1.3)
    xlabel('$\sqrt{N}$','interpreter','latex')
    ylabel('$\|e\|_\infty/\|u\|_\infty$','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    leg = legend('D-RBF-PU, Smooth','D-RBF-PU, ConstGen', 'RBF-HFD','Location','northeast');
    set(leg,'Interpreter','latex');
    set(gcf, 'Position', [300 300 350 400])
    title(sprintf('PHS7+P3, %s',StrBound),'Interpreter','latex');
    xlim([10 2*10^2])
    xticks([10 20 50 100 200])
    text(35,10^-3,num2str(sprintf('%1.2f\n',Orders(:))),'Interpreter','latex');
    
    % plot cpu times
    figure;
    loglog(sqrt(N),SetupTime(:,1),'rs-','MarkerSize',6,'MarkerFaceColor','r','LineWidth',1.3)
    hold on
    loglog(sqrt(N),SetupTime(:,2),'go-','MarkerSize',6,'MarkerFaceColor','g','LineWidth',1.3)
    hold on
    loglog(sqrt(N),SetupTime(:,3),'k>-','MarkerSize',6,'MarkerFaceColor','k','LineWidth',1.3)
    xlabel('$\sqrt{N}$','interpreter','latex')
    ylabel('CPU time (sec.)','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    leg = legend('D-RBF-PU, Smooth','D-RBF-PU, ConstGen', 'RBF-HFD','Location','southeast');
    set(leg,'Interpreter','latex');
    set(gcf, 'Position', [300 300 400 370])
    title('Setup times','Interpreter','latex');
    xlim([10 2*10^2])
    xticks([10 20 50 100 200])
    
    figure;
    loglog(sqrt(N),SolveTime(:,1),'rs-','MarkerSize',6,'MarkerFaceColor','r','LineWidth',1.3)
    hold on
    loglog(sqrt(N),SolveTime(:,2),'go-','MarkerSize',6,'MarkerFaceColor','g','LineWidth',1.3)
    hold on
    loglog(sqrt(N),SolveTime(:,3),'k>-','MarkerSize',6,'MarkerFaceColor','k','LineWidth',1.3)
    xlabel('$\sqrt{N}$','interpreter','latex')
    ylabel('CPU time (sec.)','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    leg = legend('D-RBF-PU, Smooth','D-RBF-PU, ConstGen', 'RBF-HFD','Location','southeast');
    set(leg,'Interpreter','latex');
    set(gcf, 'Position', [300 300 400 370])
    title('Solving times','Interpreter','latex');
    xlim([10 2*10^2])
    xticks([10 20 50 100 200])
end
