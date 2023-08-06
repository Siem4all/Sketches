import os
import pickle
import matplotlib.pyplot as plt

def plot_counters():
    """
    Load the dictionaries from the pcl files and plot the Normalized_Morris_RMSE in Y-axis and the counter Width in X-axis f
    or both the counter types to compare and understand their difference easily.
    """
    Morris_increments = []
    with open('../res/pcl_files/Morris_increments.pcl', 'rb') as f:
        while True:
            try:
                my_dict = pickle.load(f)
                Morris_increments.append(my_dict)
            except EOFError:
                break
    realCounter_increments = []
    with open('../res/pcl_files/realCounter_increments.pcl', 'rb') as f:
        while True:
            try:
                my_dict = pickle.load(f)
                realCounter_increments.append(my_dict)
            except EOFError:
                break
    CEDAR_increments = []
    with open('../res/pcl_files/CEDAR_increments.pcl', 'rb') as f:
        while True:
            try:
                my_dict = pickle.load(f)
                CEDAR_increments.append(my_dict)
            except EOFError:
                break

    # Extract the data from the dictionaries
    Morris_widths = []
    Normalized_Morris_RMSE = []
    realCounter_widths= []
    Normalized_realCounter_RMSE = []
    CEDAR_widths = []
    Normalized_CEDAR_RMSE = []
    for m in Morris_increments:
        Morris_widths.append(m['width'])
        Normalized_Morris_RMSE.append(m['Normalized_RMSE'])
    for r in realCounter_increments:
        realCounter_widths.append(r['width'])
        Normalized_realCounter_RMSE.append(r['Normalized_RMSE'])
    for c in CEDAR_increments:
        CEDAR_widths.append(c['width'])
        Normalized_CEDAR_RMSE.append(c['Normalized_RMSE'])
    if max(Normalized_Morris_RMSE) >= max(Normalized_CEDAR_RMSE):
        maxYValue = max(Normalized_Morris_RMSE) + max(Normalized_Morris_RMSE)/2
    else:
        maxYValue = max(Normalized_CEDAR_RMSE) + max(Normalized_CEDAR_RMSE)/2

    # Create a figure and axis object
    fig, axs = plt.subplots()

    # Plot the depth vs. avg_error data on the first subplot
    axs.plot(CEDAR_widths, Normalized_CEDAR_RMSE, '-o', linewidth=2, markersize=8, color='green', label='CEDAR')
    axs.plot(Morris_widths, Normalized_Morris_RMSE, '-o', linewidth=2, markersize=8, color='blue', label='Morris')
    axs.plot(realCounter_widths, Normalized_realCounter_RMSE, '-o', linewidth=2, markersize=8, color='brown', label='realCounter')
    axs.set_xlabel('Width', fontsize=12)
    axs.set_ylabel('Normalized_RMSE ', fontsize=12)
    axs.set_title('width vs. Normalized_RMSE', fontsize=14)
    axs.grid(True)
    axs.set_ylim([0, max(Normalized_CEDAR_RMSE) + max(Normalized_CEDAR_RMSE)/5])
    axs.set_xlim([0, max(CEDAR_widths) + 1])

    # Add legends with titles to the subplots
    axs.legend()
    # Save the figure as a JPEG image in a specific directory
    save_path = os.path.join('../res/images', 'normalized_RMSE.jpg')

    # Check if the file already exists and remove it if it does
    if os.path.exists(save_path):
        os.remove(save_path)

    plt.savefig(save_path, dpi=300, bbox_inches='tight')

    # Show the plot
    plt.show()
