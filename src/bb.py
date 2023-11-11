import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
import pickle
pclFileName='RdNMSE_2depth_12bits'
data = [
    {'mode': 'SEAD stat', 'numCntrs': 4, 'Avg': 1.9265724765029656e-06, 'Lo': 1.9265724765029656e-06, 'Hi': 1.9265724765029656e-06},
    {'mode': 'SEAD stat', 'numCntrs': 16, 'Avg': 4.6070930433444664e-07, 'Lo': 4.6070930433444664e-07, 'Hi': 4.6070930433444664e-07},
    {'mode': 'SEAD stat', 'numCntrs': 64, 'Avg': 1.0797237154678404e-07, 'Lo': 1.0797237154678404e-07, 'Hi': 1.0797237154678404e-07},
    {'mode': 'F2P', 'numCntrs': 4, 'Avg': 1.991837331938293e-06, 'Lo': 1.991837331938293e-06, 'Hi': 1.991837331938293e-06},
    {'mode': 'F2P', 'numCntrs': 16, 'Avg': 4.5308619742228003e-07, 'Lo': 4.5308619742228003e-07, 'Hi': 4.5308619742228003e-07},
    {'mode': 'F2P', 'numCntrs': 64, 'Avg': 1.0835588325709969e-07, 'Lo': 1.0835588325709969e-07, 'Hi': 1.0835588325709969e-07},
    {'mode': 'CEDAR', 'numCntrs': 4, 'Avg': 1.9402111004479975e-06, 'Lo': 1.9402111004479975e-06, 'Hi': 1.9402111004479975e-06},
    {'mode': 'CEDAR', 'numCntrs': 16, 'Avg': 4.6007496731622586e-07, 'Lo': 4.6007496731622586e-07, 'Hi': 4.6007496731622586e-07},
    {'mode': 'CEDAR', 'numCntrs': 64, 'Avg': 1.0868835918903677e-07, 'Lo': 1.0868835918903677e-07, 'Hi': 1.0868835918903677e-07},
    {'mode': 'Morris', 'numCntrs': 4, 'Avg': 1.976681623397342e-06, 'Lo': 1.976681623397342e-06, 'Hi': 1.976681623397342e-06},
    {'mode': 'Morris', 'numCntrs': 16, 'Avg': 4.519877163819317e-07, 'Lo': 4.519877163819317e-07, 'Hi': 4.519877163819317e-07},
    {'mode': 'Morris', 'numCntrs': 64, 'Avg': 1.0708445143696941e-07, 'Lo': 1.0708445143696941e-07, 'Hi': 1.0708445143696941e-07},
    {'mode': 'RealCntr', 'numCntrs': 4, 'Avg': 1.9764812167089565e-06, 'Lo': 1.9764812167089565e-06, 'Hi': 1.9764812167089565e-06},
    {'mode': 'RealCntr', 'numCntrs': 16, 'Avg': 4.581019395580588e-07, 'Lo': 4.581019395580588e-07, 'Hi': 4.581019395580588e-07},
    {'mode': 'RealCntr', 'numCntrs': 64, 'Avg': 1.091234207603814e-07, 'Lo': 1.091234207603814e-07, 'Hi': 1.091234207603814e-07}
]
MARKER_SIZE            = 20
MARKER_SIZE_SMALL      = 1
LINE_WIDTH             = 5
LINE_WIDTH_SMALL       = 1
FONT_SIZE              = 20
FONT_SIZE_SMALL        = 5
LEGEND_FONT_SIZE       = 14
LEGEND_FONT_SIZE_SMALL = 5
setPltParams = lambda size='large': matplotlib.rcParams.update({'font.size'      : FONT_SIZE,
                                                                          'legend.fontsize': LEGEND_FONT_SIZE,
                                                                          'xtick.labelsize': FONT_SIZE,
                                                                          'ytick.labelsize': FONT_SIZE,
                                                                          'axes.labelsize' : FONT_SIZE,
                                                                          'axes.titlesize' : FONT_SIZE, }) if (
            size == 'large') else matplotlib.rcParams.update({
        'font.size'      : FONT_SIZE_SMALL,
        'legend.fontsize': LEGEND_FONT_SIZE_SMALL,
        'xtick.labelsize': FONT_SIZE_SMALL,
        'ytick.labelsize': FONT_SIZE_SMALL,
        'axes.labelsize' : FONT_SIZE_SMALL,
        'axes.titlesize' : FONT_SIZE_SMALL
    })
colorOfMode = {'F2P'         : 'green',
                             'RealCntr': 'blue',
                            'CEDAR'       : 'brown',
                            'Morris'      : 'red',
                            'SEAD stat'    : 'Purple'}

markerOfMode = {'F2P'         : 'o',
                              'RealCntr': 'v',
                             'CEDAR'       : '<',
                             'Morris'      : '>',
                             'SEAD stat'    : '*'}
# Extract x and y values for each mode
modes = list(set([entry['mode'] for entry in data]))
setPltParams()  # set the plot's parameters (formats of lines, markers, legends etc.).
_, ax = plt.subplots()
for mode in modes:
    x = [entry['numCntrs'] for entry in data if entry['mode'] == mode]
    y = [entry['Avg'] for entry in data if entry['mode'] == mode]
    y_lo = [entry['Lo'] for entry in data if entry['mode'] == mode]
    y_avg=[entry['Avg'] for entry in data if entry['mode'] == mode]
    y_hi = [entry['Hi'] for entry in data if entry['mode'] == mode]
    numCntrs = [entry['numCntrs'] for entry in data if entry['mode'] == mode]
    ax.plot((numCntrs, numCntrs), (y_lo, y_hi), color=colorOfMode[mode])  # Plot the conf' interval line
    ax.plot(numCntrs, y_avg, color=colorOfMode[mode], marker=markerOfMode[mode],
                     markersize=MARKER_SIZE, linewidth=LINE_WIDTH, label=mode, mfc='none')
    ax.set_xticks(numCntrs)
    ax.set_yticks(y_avg)

plt.xlabel('numCntrs')
plt.ylabel('AvgRdError')
plt.title('AvgRdError vs. numCntrs')
# Set the exact values on the x-axis
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), fontsize=LEGEND_FONT_SIZE)
plt.savefig('../res/{}.pdf'.format(pclFileName), bbox_inches='tight')

