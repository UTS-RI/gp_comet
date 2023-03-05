# Visualise the given track (track number given as argument of script)

import sys
import os
import numpy as np
import matplotlib.pyplot as plt


def main(path, track_number):
    
    trajectory_path = path + "track_" + str(track_number) + ".txt"
    event_undist_path = path + "track_" + str(track_number) + "_undistorted_events.txt"
    event_path = path + "track_" + str(track_number) + "_raw_events.txt"

    # Check if the files exist
    if not os.path.isfile(trajectory_path):
        print("Trajectory file " + trajectory_path + " does not exist")
        return
    if not os.path.isfile(event_undist_path):
        print("Event file " + event_undist_path + " does not exist")
        return
    if not os.path.isfile(event_path):
        print("Event file " + event_path + " does not exist")
        return

    # Load the files
    trajectory = np.loadtxt(trajectory_path)
    if(trajectory.size == 0):
        print("Trajectory file " + trajectory_path + " is empty")
        return
    event_undist = np.loadtxt(event_undist_path)
    event = np.loadtxt(event_path)

    # Create a simgle image to plot in order, the raw event, the undistorted event and the trajectory
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(event[:, 1], event[:, 2], s=0.1, c='r')
    ax.scatter(event_undist[:, 1], event_undist[:, 2], s=0.1, c='b', alpha=0.3)
    ax.plot(trajectory[:, 1], trajectory[:, 2], c='g')

    ax.set_title("Track " + str(track_number) + " of length " + str(round(trajectory[-1, 0] - trajectory[0,0],3)) + "s has " + str(np.size(event_undist,0)) + " events")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend(["Raw events", "Undistorted events", "Trajectory"], markerscale=10)
    # Set aspect ratio to 1
    ax.set_aspect('equal')
    plt.show()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])

