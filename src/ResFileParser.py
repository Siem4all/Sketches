import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker, cms
import pickle

MARKER_SIZE            = 16
MARKER_SIZE_SMALL      = 1
LINE_WIDTH             = 3
LINE_WIDTH_SMALL       = 1
FONT_SIZE              = 20
FONT_SIZE_SMALL        = 5
LEGEND_FONT_SIZE       = 14
LEGEND_FONT_SIZE_SMALL = 5


class ResFileParser(object):
    """
    Parse res files, and generate plots from them.
    """

    # Set the parameters of the plot (sizes of fonts, legend, ticks etc.).
    # mfc='none' makes the markers empty.
    setPltParams = lambda self, size='large': matplotlib.rcParams.update({'font.size'      : FONT_SIZE,
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

    def __init__(self):
        """
        Initialize a pcl_file_parser, used to parse result files, and generate plots.
        """
        # List of algorithms' names, used in the plots' legend, for the dist' case
        self.labelOfMode = {}

        # The colors used for each alg's plot, in the dist' case
        self.colorOfMode = {'F2P'         : 'green',
                             'RealCntr': 'blue',
                            'CEDAR'       : 'brown',
                            'Morris'      : 'red',
                            'SEAD stat'    : 'Purple'}

        # The markers used for each alg', in the dist' case
        self.markerOfMode = {'F2P'         : 'o',
                              'RealCntr': 'v',
                             'CEDAR'       : '<',
                             'Morris'      : '>',
                             'SEAD stat'    : '*'}
        self.points      = []

    def rdRes(self):
        """
        Given a RdRmse.res, read all the data itcontains into self.points
        """
        for index in range(len(cms.data)):
            self.points=cms.data[index]
    def NRMSEVsWidthPlot(self, modes):
        """
        Generate a plot showing the Normalized_RMSE vs. width.
        """

        self.setPltParams()  # set the plot's parameters (formats of lines, markers, legends etc.).
        _, ax = plt.subplots()
        for mode in modes:
            pointsOfThisMode = [point for point in self.points if point['mode'] == mode]
            if pointsOfThisMode == []:
                print(f'No points found for mode {mode}')
                continue
            numCntrs = [point['numCntrs'] for point in pointsOfThisMode]
            y_lo = [point['Lo'] for point in pointsOfThisMode]
            y_avg=[point['Avg'] for point in pointsOfThisMode]
            y_hi = [point['Hi'] for point in pointsOfThisMode]
            ax.plot((numCntrs, numCntrs), (y_lo, y_hi), color=self.colorOfMode[mode])  # Plot the conf' interval line
            ax.plot(numCntrs, y_avg, color=self.colorOfMode[mode], marker=self.markerOfMode[mode],
                     markersize=MARKER_SIZE, linewidth=LINE_WIDTH, label=mode, mfc='none')
            ax.set_xticks(numCntrs)

        plt.xlabel('numCntrs')
        plt.ylabel('AvgRdError')
        plt.title('AvgRdError vs. numCntrs')
        # Set the exact values on the x-axis
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), fontsize=LEGEND_FONT_SIZE)
        plt.show()


if __name__ == '__main__':
    counter_modes = ['SEAD stat', 'F2P', 'CEDAR', 'Morris', 'RealCntr']
    parser = ResFileParser()
    # Read in data from a PCL file
    parser.rdRes()
    # Generate a plot showing normRmseAvg versus number of counter for each counter modes.
    parser.NRMSEVsWidthPlot(modes=counter_modes)
