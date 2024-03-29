# GP COMET


**G**aussian **P**rocess for **CO**mpentation of **M**otion for **E**vent-based **T**tracking

This is the implementation of the IEEE ICRA2023 paper __Continuous-Time Gaussian Process Motion-Compensation for Event-vision Pattern Tracking with Distance Fields__ available [here](https://arxiv.org/abs/2303.02672).

Given seeds (x,y,t), this piece of software performs motion compensation and pattern tracking using solely event data.


__This code is not optimised. If you are interrested to make it run faster or to improve the performance, feel free to send me an email, I'll be happy to help and provide further insight.__

## How to compile

### Dependencies

```
sudo apt-get install libeigen3-dev
sudo apt-get install libboost-all-dev
sudo apt-get install libceres-dev
sudo apt-get install libyaml-cpp-dev
sudo apt-get install libopencv-dev
```

Known working configuration

- Ubuntu 20.04
- Eigen 3.3.7
- Boost 1.71.0
- Ceres 1.14.0
- OpenCV 4.6.0
- OpenMP 4.5

It probably works well with other versions but we never tested it.


### Compiling
```
mkdir build
cd build
cmake ../
make -j4
```


## How to run

From the `build` folder
```
./app/gp_comet -c ../config.yaml
```

The file `config.yaml` contains the parameters and the path to both the input data and the output folder.
Explanations are provided as comments in the file itself.

The input data must follow the format of [The Event-Camera Dataset and Simulator](https://rpg.ifi.uzh.ch/davis_data.html).
That is:

- `events.txt`: One event per line (timestamp x y polarity)
- `calib.txt`: Camera parameters (fx fy cx cy k1 k2 p1 p2 k3)

The track seeds can be provided as a `.txt` file with one track per line (id timestamp x y),

The tracks generated by our code are written in the `results_path/` folder:

- `track_[id].txt`: Trajectory of the track in the image
- `track_[id]_raw_events.txt` (if `write_events: true`): Events used for the track generation
- `track_[id]_undistorted_events.txt` (if `write_events: true`): Events of the track after motion compensation

Each line of these three files is (timestamp x y polarity).
Note that polarity is not used in this work.

We provided sample data in `data/sample_data.zip`.
It is a short sample from the `poster_6dof` dataset in [The Event-Camera Dataset and Simulator](https://rpg.ifi.uzh.ch/davis_data.html).
To run this software, you just need to unzip the sample data archive in place, and run the command `./app/gp_comet -c ../config.yaml` with the provided configuration file.

In its current state, the configuration file has the same parameters as used in the paper except for `downsample` (comment that line and you have the same parameters as in the evaluation section of the paper)

EDIT 29/05/2023: This directory now includes a "preprocessing" step based on the lumped matrix approximation for the GP of the SE2 motion compensation. It seems to help convergences a bit. This approximation can be used on its own (as the main SE2 compensation mechanism) by switching the `approx_only` to `true` for fast computation (I did not evaluate the performance of the `approx_only` version other than qualitatively but it is definitely faster).



## Other utils

In the folder `other_utils`, we provide a python script for the visualisation of a track from the results. 
It can be used as follows
```
python visualise_track.py ../results/ id
```
With `id` the number of the track you want to display.
The file `requierements.txt` contains the dependencies for that script.


