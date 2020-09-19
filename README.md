# TrajectoryEstimation
Fast trajectory MAP and its corresponding application in path planning and visual tracking.  
Tips: part of the tests need path configuration. Please replace the relative path into your absolute path (this program) if necessary.  
For more information, please contact author by: williamli_pro@163.com

## Citation
Welcome to use this code in your project, this paper is under review.

## Abstract
Trajectory estimation of maneuvering objects is applied in numerous tasks like navigation, path planning and visual tracking. Many previous works get impressive results in the strictly controlled condition with accurate prior statistics and dedicated dynamic model for certain object. But in challenging conditions without dedicated dynamic model and precise prior statistics, the performance of these methods significantly declines. To solve the problem, a stochastic nonlinear model called the power-limited steering model (PLS) is proposed to describe the motion of non-cooperative object. It is a natural combination of instantaneous power and instantaneous angular velocity, relying on the nonlinearity to achieve the change of states. And the renormalization group is introduced to compensate the nonlinear effect of perturbation in PLS model. For robust and efficient trajectory estimation, an adaptive trajectory estimation (AdaTE) algorithm is proposed. By updating the statistics and truncation time online, it corrects the estimation error caused by biased prior statistics and observation drift, while reducing the computational complexity lower than O(n). The experiment of trajectory estimation demonstrates the convergence and robust of AdaTE. Other experiments demonstrate through slight modification, AdaTE can also be applied to local navigation in random obstacle environment, and trajectory optimization in visual tracking.

## Results
![Trajectory estimation from 3D observations (under noise)](https://github.com/WilliamLiPro/TrajectoryEstimation/tree/master/present-result/trajectory-estimation.png)

![Local navigation and obstacle avoidance](https://github.com/WilliamLiPro/TrajectoryEstimation/tree/master/present-result/online-path-planning.png)

![visual tracking](https://github.com/WilliamLiPro/TrajectoryEstimation/tree/master/present-result/visual-tracking-filtering.png)
