import os
import pickle
import matplotlib.pyplot as plt

def plot_counters():
    # Load the dictionaries from the pcl files
    first_increments = []
    with open('../res/pcl_files/first_increments.pcl', 'rb') as f:
        while True:
            try:
                my_dict = pickle.load(f)
                first_increments.append(my_dict)
            except EOFError:
                break

    # Extract the data from the dictionaries
    numOfIncrements = []
    Normalized_RMSE = []
    counters = first_increments[1]['counters']
    for f in first_increments:
        numOfIncrements.append(f['numOfIncrements'])
        Normalized_RMSE.append(f['Normalized_RMSE'])

    # Create a figure and axis object
    fig, axs = plt.subplots()

    # Plot the data on the axis
    axs.plot(numOfIncrements, Normalized_RMSE, '-o', linewidth=2, markersize=8, color='blue', label=f'{counters} counters')

    # Add labels and a title to the plot
    axs.set_xlabel('num Of Increments', fontsize=12)
    axs.set_ylabel('Normalized RMSE', fontsize=12)
    axs.set_title('num Of Increments vs. Normalized RMSE', fontsize=14)

    # Set the axis limits and add a grid
    axs.set_ylim([0, max(Normalized_RMSE) +max(Normalized_RMSE)/2])
    axs.set_xlim([0, max(numOfIncrements) + max(numOfIncrements)//10])
    axs.grid(True)

    # Add a legend with a title
    axs.legend()

    # Adjust the spacing between the subplots

    # Save the figure as a JPEG image in a specific directory
    save_path = os.path.join('../res/images', 'normalized_RMSE.jpg')

    # Check if the file already exists and remove it if it does
    if os.path.exists(save_path):
        os.remove(save_path)

    plt.savefig(save_path, dpi=300, bbox_inches='tight')

    # Show the plot
    plt.show()
