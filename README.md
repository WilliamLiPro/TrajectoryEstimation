# TrajectoryEstimation
Fast trajectory MAP and its corresponding application in path planning and visual tracking.  
Tips: part of the tests need path configuration. Please replace the relative path into your absolute path (this program) if necessary.  
For more information, please contact author by: williamli_pro@163.com

## Citation
Welcome to use this code in your project, this paper is under review.

## Abstract
Trajectory estimation under regional correlations is applied in numerous tasks like dynamic path planning, navigation and tracking. Many previous works get impressive results in the strictly controlled condition with accurate prior statistics and specific dynamic model for certain object. While in a more challenging situation without specific dynamic model and accurate prior statistics, the performance of these methods significantly declines. To estimate the trajectory without specific dynamic model, we propose a general model called the power-limited steering model (PLS), which is a natural combination of instantaneous power and instantaneous angular velocity. As a nonlinear model it needs less parameter to describe the switching of states compared with linear Markov process model. Then we derive the corresponding form in discrete time and compensate the nonlinear effect of perturbation with renormalization group. To overcome the biased prior statistics, the observations drift and linear growing computation in trajectory estimation, we propose the adaptive trajectory estimation (AdaTE) where the online updated statistics and the adaptive truncation time is applied. The experiment of typical trajectory estimation demonstrates that compared with EKF and UKF, AdaTE is more robust to the biased prior statistics and observation drift. Another task shows that with slight modification, AdaTE can be used in path planning that needs obstacle avoidance. Finally, we exhibit its potential application on trajectory optimization in visual tracking.

## Results
![Trajectory estimation from 3D observations (under noise)](https://github.com/WilliamLiPro/TrajectoryEstimation/tree/master/present-result/trajectory-estimation.png)

![Local navigation and obstacle avoidance](https://github.com/WilliamLiPro/TrajectoryEstimation/tree/master/present-result/online-path-planning.png)

![visual tracking](https://github.com/WilliamLiPro/TrajectoryEstimation/tree/master/present-result/visual-tracking-filtering.png)
