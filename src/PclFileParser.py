import os
import pickle
import matplotlib.pyplot as plt

def plot_counters():
    # Load the dictionaries from the pcl files
    Morris_increments = []
    with open('../res/pcl_files/Morris_increments.pcl', 'rb') as f:
        while True:
            try:
                my_dict = pickle.load(f)
                Morris_increments.append(my_dict)
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
    widths_Morris = []
    Normalized_Morris_RMSE = []
    widths_CEDAR = []
    Normalized_CEDAR_RMSE = []
    for m in Morris_increments:
        widths_Morris.append(m['width'])
        Normalized_Morris_RMSE.append(m['Normalized_RMSE'])
    for c in CEDAR_increments:
        widths_CEDAR.append(c['width'])
        Normalized_CEDAR_RMSE.append(c['Normalized_RMSE'])
    if max(Normalized_Morris_RMSE) >= max(Normalized_CEDAR_RMSE):
        maxYValue = max(Normalized_Morris_RMSE) + max(Normalized_Morris_RMSE)/2
    else:
        maxYValue = max(Normalized_CEDAR_RMSE) + max(Normalized_CEDAR_RMSE)/2

    # Create a figure and axis object
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

    # Plot the depth vs. avg_error data on the first subplot
    axs[0].plot(widths_Morris, Normalized_Morris_RMSE, '-o', linewidth=2, markersize=8, color='blue', label='Morris')
    axs[0].set_xlabel('Width', fontsize=12)
    axs[0].set_ylabel('Normalized_RMSE ', fontsize=12)
    axs[0].set_title('width vs. Normalized_RMSE', fontsize=14)
    axs[0].grid(True)
    axs[0].set_ylim([0, max(Normalized_Morris_RMSE) + max(Normalized_Morris_RMSE)/2])
    axs[0].set_xlim([0, max(widths_Morris) + 1])

    # Plot the width vs. avg_error data on the second subplot
    axs[1].plot(widths_CEDAR, Normalized_CEDAR_RMSE, '-o', linewidth=2, markersize=8, color='green', label='CEDAR')
    axs[1].set_xlabel('Width', fontsize=12)
    axs[1].set_ylabel('Normalized_RMSE', fontsize=12)
    axs[1].set_title('Width vs. Normalized_RMSE', fontsize=14)
    axs[1].grid(True)
    axs[1].set_ylim([0, max(Normalized_CEDAR_RMSE) + max(Normalized_CEDAR_RMSE)/2])
    axs[1].set_xlim([0, max(widths_CEDAR) + 1])
    axs[1].grid(True)

    # Add legends with titles to the subplots
    axs[0].legend()
    axs[1].legend()

    # Adjust the spacing between the subplots
    plt.subplots_adjust(wspace=0.4)

    # Save the figure as a JPEG image in a specific directory
    save_path = os.path.join('../res/images', 'normalized_RMSE.jpg')

    # Check if the file already exists and remove it if it does
    if os.path.exists(save_path):
        os.remove(save_path)

    plt.savefig(save_path, dpi=300, bbox_inches='tight')

    # Show the plot
    plt.show()
