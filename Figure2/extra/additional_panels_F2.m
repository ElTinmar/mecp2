clear all
load("fig2.mat")

[DFF_properties, Binary_properties] = response_properties(...
        neuralActivityAroundStim,...
        resp_perStim,...
        resp_allStim,...
        6);   

f = figure('Position',[0 0 1080 1080]);

ax1 = subplot(3,3,1);
quantity{1} = DFF_properties{1}.perStim.pk;
quantity{2} = DFF_properties{2}.perStim.pk;
plot_curve(ax1,quantity,[0 0.6],'<peak \Delta{}F/F>_{fish}')
ax2 = subplot(3,3,2);
quantity{1} = DFF_properties{1}.perStim.itg;
quantity{2} = DFF_properties{2}.perStim.itg;
plot_curve(ax2,quantity,[0 5],'<\Sigma{}\Delta{}F/F>_{fish}')
ax3 = subplot(3,3,3);
quantity{1} = DFF_properties{1}.perStim.fwhm;
quantity{2} = DFF_properties{2}.perStim.fwhm;
plot_curve(ax3,quantity,[0 2.5],'<duration (s)_{\Delta{}F/F}>_{fish}')

ax4 = subplot(3,3,4);
quantity{1} = Binary_properties{1}.perStim.pk;
quantity{2} = Binary_properties{2}.perStim.pk;
plot_curve(ax4,quantity,[0 0.4],'<peak resp. prob.>_{fish}')
ax5 = subplot(3,3,5);
quantity{1} = Binary_properties{1}.perStim.itg;
quantity{2} = Binary_properties{2}.perStim.itg;
plot_curve(ax5,quantity,[0 3],'<number events>_{fish}')
ax6 = subplot(3,3,6);
quantity{1} = Binary_properties{1}.perStim.fwhm;
quantity{2} = Binary_properties{2}.perStim.fwhm;
plot_curve(ax6,quantity,[0 2.5],'<duration (s)_{binary}>_{fish}')

ax7 = subplot(3,3,7);
quantity{1} = DFF_properties{1}.perStim.prop_resp;
quantity{2} = DFF_properties{2}.perStim.prop_resp;
plot_curve(ax7,quantity,[0 10],'<% resp. cells>_{fish}')

print(f,'additionnal_fig2.pdf','-dpdf','-bestfit');

