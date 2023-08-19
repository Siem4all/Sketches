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
    realcounter_data = [entry for entry in WrRmse if entry['mode'] == 'realCounter']
    # Extract width and Normalized_RMSE values for each mode
    morris_width = [entry['width'] for entry in morris_data]
    morris_rmse = [entry['Normalized_RMSE'] for entry in morris_data]

    realcounter_width = [entry['width'] for entry in realcounter_data]
    realcounter_rmse = [entry['Normalized_RMSE'] for entry in realcounter_data]

    # Plotting with markers and connecting points based on width
    morris_plot = plt.scatter(morris_width, morris_rmse, label='Morris', marker='o', color='blue')
    realcounter_plot = plt.scatter(realcounter_width, realcounter_rmse, label='realCounter', marker='s', color='green')

    # Connect points based on width
    plt.plot((morris_width, morris_width), (morris_rmse, realcounter_rmse), color='pink')
    plt.plot(morris_width, morris_rmse, linestyle='-', color='blue')
    plt.plot(realcounter_width, realcounter_rmse, linestyle='-', color='green')

    plt.xlabel('Width')
    plt.ylabel('Normalized_RMSE')
    plt.title('Normalized_RMSE vs. Width')

    # Create the legend with handles and labels
    handles = [morris_plot, realcounter_plot]
    labels = ['Morris', 'realCounter']
    plt.legend(handles, labels)

    # Save the plot as a PDF file
    with PdfPages('../res/WrRMSE.pdf') as pdf:
        pdf.savefig(bbox_inches='tight')
        plt.close()

