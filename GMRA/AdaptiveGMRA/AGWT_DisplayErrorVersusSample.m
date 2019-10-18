ErrOpts = 0;


Nsample = DataErrorSample.Nsample;

figure
if ErrOpts == 0
    plot(Nsample,DataErrorSample.UniformAbsError,'b*--',Nsample,DataErrorSample.AdaptiveAbsError,'rp--','MarkerSize',8)
elseif ErrOpts == 1
    plot(Nsample,DataErrorSample.UniformReError,'b*--',Nsample,DataErrorSample.AdaptiveReError,'rp--','MarkerSize',8)
end
title('Error versus sample size','FontSize',15)
legend('Uniform','Adaptive')
set(gca,'FontSize',14,'FontWeight','bold')


%% Least square fitting
figure
if ErrOpts == 0
    UniformFit    = polyfit(log10(Nsample),log10(DataErrorSample.UniformAbsError),1);
    AdaptiveFit   = polyfit(log10(Nsample),log10(DataErrorSample.AdaptiveAbsError),1);
    plot(log10(Nsample),log10(DataErrorSample.UniformAbsError),'b*--',log10(Nsample),log10(DataErrorSample.AdaptiveAbsError),'rp--','MarkerSize',8)
    legend(['Uniform:  y=',poly2str(UniformFit,'x')],['Adaptive: y=',poly2str(AdaptiveFit,'x')])
    ylabel('log10(Absolute error)')
elseif ErrOpts == 1
    UniformFit    = polyfit(log10(Nsample),log10(DataErrorSample.UniformReError),1);
    AdaptiveFit   = polyfit(log10(Nsample),log10(DataErrorSample.AdaptiveReError),1);
    plot(log10(Nsample),log10(DataErrorSample.UniformReError),'b*--',log10(Nsample),log10(DataErrorSample.AdaptiveReError),'rp--','MarkerSize',8)
    legend(['Uniform:  y=',poly2str(UniformFit,'x')],['Adaptive: y=',poly2str(AdaptiveFit,'x')])
    ylabel('log10(Relative error)')
end
title('log10(Error) versus log10(Sample size)','FontSize',15)
xlabel('log10(Sample size)')
set(gca,'FontSize',14,'FontWeight','bold')




