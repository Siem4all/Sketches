import os
import pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def plot_counters():
    """
    Load the dictionaries from the pcl files and plot the Normalized_Morris_RMSE in Y-axis and the counter Width in X-axis for both the counter types to compare and understand their difference easily.
    """
    WrRmse = []
    with open('../res/pcl_files/WrRmse.pcl', 'rb') as f:
        while True:
            try:
                my_dict = pickle.load(f)
                WrRmse.append(my_dict)
            except EOFError:
                break

    # Separate data for each mode
    morris_data = [entry for entry in WrRmse if entry['mode'] == 'Morris']
    cedar_data = [entry for entry in WrRmse if entry['mode'] == 'CEDAR']
    realcounter_data = [entry for entry in WrRmse if entry['mode'] == 'realCounter']
    # Extract width and Normalized_RMSE values for each mode
    morris_width = [entry['width'] for entry in morris_data]
    morris_rmse = [entry['Normalized_RMSE'] for entry in morris_data]

    cedar_width = [entry['width'] for entry in cedar_data]
    cedar_rmse = [entry['Normalized_RMSE'] for entry in cedar_data]

    realcounter_width = [entry['width'] for entry in realcounter_data]
    realcounter_rmse = [entry['Normalized_RMSE'] for entry in realcounter_data]

    # Plotting with markers and connecting points based on width
    plt.scatter(morris_width, morris_rmse, label='Morris', marker='o')
    plt.scatter(cedar_width, cedar_rmse, label='CEDAR', marker='s')
    plt.scatter(realcounter_width, realcounter_rmse, label='realCounter', marker='^')

    # Connect points based on width
    plt.plot((morris_width, cedar_width, realcounter_width), (morris_rmse, cedar_rmse, realcounter_rmse), color='pink')
    plt.plot(morris_width, morris_rmse, linestyle='-', color='blue')
    plt.plot(cedar_width, cedar_rmse, linestyle='-', color='orange')
    plt.plot(realcounter_width, realcounter_rmse, linestyle='-', color='green')

    plt.xlabel('Width')
    plt.ylabel('normRmseAvg')
    plt.title('normRmseAvg vs. Width')
    plt.legend()

    # Save the plot as a PDF file
    with PdfPages('../res/WrRMSE.pdf') as pdf:
        pdf.savefig(bbox_inches='tight')
        plt.close()
