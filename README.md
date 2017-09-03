# Writeup for P2 of Self-Driving Car Nanodegree Term2: Unscented Kalman Filter
[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)

### Matt, Min-gyu, Kim
---

[image1]: ./result/UKF_Process.png "UKF_PROCESS"
[image2]: ./result/result_unscented_kalman_filter.png "SCREENSHOT"
[image3]: ./result/NIS.png "NIS"

### Brief Description

I built the [Unscented Kalman Filter](http://www.cs.unc.edu/~welch/kalman/media/pdf/Julier1997_SPIE_KF.pdf) which takes inputs from sensor fusion - LiDAR and RADAR). The following figure is the processign flow the filter.

![alt_text][image1]

All of the process are implemented in `ukf.cpp`

---

### Simulation

The compiled program can communicate with the simulator by [uWebSocket](https://github.com/uNetworking/uWebSockets). Here is the captured image after the simulation is done for `Dataset 1`.

![alt_text][image2]

### Parameter Tunning : Process Noise

In this project, two process noises - standard deviations of `longitudinal acceleration` and `yaw acceleration` experimentally. To check the consistency of the parameters, I calculated the NIS(Noirmalized Innovation Squared) between the predicted state and the ground truth. The following is the plot of the NIS values earned from final parameters.


![alt_text][image3]
