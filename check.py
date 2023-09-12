import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
import pickle

MARKER_SIZE = 16
MARKER_SIZE_SMALL = 1
LINE_WIDTH = 3
LINE_WIDTH_SMALL = 1
FONT_SIZE = 20
FONT_SIZE_SMALL = 5
LEGEND_FONT_SIZE = 14
LEGEND_FONT_SIZE_SMALL = 5


class ResFileParser(object):
    """
    Parse result files, and generate plots from them.
    """

    # Set the parameters of the plot (sizes of fonts, legend, ticks etc.).
    # mfc='none' makes the markers empty.
    setPltParams = lambda self, size='large': matplotlib.rcParams.update({'font.size': FONT_SIZE,
                                                                          'legend.fontsize': LEGEND_FONT_SIZE,
                                                                          'xtick.labelsize': FONT_SIZE,
                                                                          'ytick.labelsize': FONT_SIZE,
                                                                          'axes.labelsize': FONT_SIZE,
                                                                          'axes.titlesize': FONT_SIZE, }) if (
            size == 'large') else matplotlib.rcParams.update({
        'font.size': FONT_SIZE_SMALL,
        'legend.fontsize': LEGEND_FONT_SIZE_SMALL,
        'xtick.labelsize': FONT_SIZE_SMALL,
        'ytick.labelsize': FONT_SIZE_SMALL,
        'axes.labelsize': FONT_SIZE_SMALL,
        'axes.titlesize': FONT_SIZE_SMALL
    })

    def __init__(self):
        """
        Initialize a Res_file_parser, used to parse result files, and generate plots.
        """
        # List of algorithms' names, used in the plots' legend, for the dist' case
        self.labelOfMode = {}

        # The colors used for each alg's plot, in the dist' case
        self.colorOfMode = {'realCounter': 'green',
                            'CEDAR': 'blue',
                            'Morris': 'red'}

        # The markers used for each alg', in the dist' case
        self.markerOfMode = {'realCounter': 'v',
                             'CEDAR': '<',
                             'Morris': '>'}
        self.points = []

    def rdPcl(self):
        """
        Given a RdRmse.pcl, read all the data it contains into self.points
        """
        pclFile = open('../res/pcl_files/RdRmse.pcl', 'rb')
        while True:
            try:
                self.points.append(pickle.load(pclFile))
            except EOFError:
                break

    def genErVsCntrSizePlot(self, modes=['realCounter', 'CEDAR', 'Morris']):
        """
        Generate a plot showing the error as a function of the counter's size.
        """

        self.setPltParams()  # set the plot's parameters (formats of lines, markers, legends etc.).
        _, ax = plt.subplots()
        all_keys = []
        all_values = {}
        index = 0
        for mode in modes:
            all_keys.append(f'{mode}_width')
            all_keys.append(f'{mode}_NRMSE')
            for index in range(index, len(all_keys)):
                pointsOfThisMode = [point for point in self.points if point['mode'] == mode]
                if pointsOfThisMode == []:
                    print(f'No points found for mode {mode}')
                    continue
                all_values[f'{all_keys[index]}'] = [point['width'] for point in pointsOfThisMode]
                all_values[f'{all_keys[index + 1]}'] = [point['Normalized_RMSE'] for point in pointsOfThisMode]
                print(all_keys[index], all_keys[index + 1])
                ax.plot(all_values[f'{all_keys[index]}'], all_values[f'{all_keys[index + 1]}'],
                        color=self.colorOfMode[mode], marker=self.markerOfMode[mode],
                        markersize=MARKER_SIZE, linewidth=LINE_WIDTH, label=f'{mode}', mfc='none')
                index += 2
                break
        plt.ylabel('Normalized_RMSE')
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), fontsize=LEGEND_FONT_SIZE)
        plt.savefig ('../res/RdRMSE.pdf', bbox_inches='tight')

v = ResFileParser()
v.rdPcl()
v.genErVsCntrSizePlot()
