# TrajectoryEstimation
Fast trajectory MAP and its corresponding application in path planning and visual tracking.  
For more information, please contact author by: williamli_pro@163.com

## Citation
Welcome to use this code in your project, the preprint can be found from: <br>
[Weipeng Li, Xiaogang Yang, Ruitao Lu, Chuan He, "Adaptive Trajectory Estimation with Power Limited Steering Model under Perturbation Compensation", arXiv:1912.03661, 2019.](https://arxiv.org/abs/1912.03661)

Tips: part of the tests need path configuration. Please replace the relative path into your absolute path (this program) if necessary.

## Abstract
Trajectory estimation under regional correlations is applied in numerous tasks like dynamic path planning, navigation and tracking. Many previous works get impressive results in the strictly controlled condition with accurate prior statistics and specific dynamic model for certain object. While in a more challenging situation without specific dynamic model and accurate prior statistics, the performance of these methods significantly declines. To estimate the trajectory without specific dynamic model, we propose a general model called the power-limited steering model (PLS), which is a natural combination of instantaneous power and instantaneous angular velocity. As a nonlinear model it needs less parameter to describe the switching of states compared with linear Markov process model. Then we derive the corresponding form in discrete time and compensate the nonlinear effect of perturbation with renormalization group. To overcome the biased prior statistics, the observations drift and linear growing computation in trajectory estimation, we propose the adaptive trajectory estimation (AdaTE) where the online updated statistics and the adaptive truncation time is applied. The experiment of typical trajectory estimation demonstrates that compared with EKF and UKF, AdaTE is more robust to the biased prior statistics and observation drift. Another task shows that with slight modification, AdaTE can be used in path planning that needs obstacle avoidance. Finally, we exhibit its potential application on trajectory optimization in visual tracking.

## Results
![Trajectory estimation from 3D observations (under noise)](https://github.com/WilliamLiPro/TrajectoryEstimation/tree/master/present-result/trajectory-estimation.png)

![Local navigation and obstacle avoidance](https://github.com/WilliamLiPro/TrajectoryEstimation/tree/master/present-result/online-path-planning.png)

![visual tracking](https://github.com/WilliamLiPro/TrajectoryEstimation/tree/master/present-result/visual-tracking-filtering.png)

## References
>[1]	E.A. Wan, and R. Van Der Merwe. “The unscented Kalman filter for nonlinear estimation”. IEEE Adaptive Systems for Signal Processing, Communications, and Control Symposium. Alberta, Canada, 2000, pp. 153–158.  
>[2]	Merwe, R. Van Der, and E. A. Wan. “The square-root unscented Kalman filter for state and parameter-estimation.” IEEE International Conference on Acoustics, Speech, and Signal Processing, 2001. Proceedings IEEE, 2002:3461-3464 vol.6.  
>[3]	Grisetti, G., Kummerle, R., Stachniss, C., &Burgard, W. (2011). “A tutorial on graph-based slam”. IEEE Intelligent Transportation Systems Magazine, 2(4), 31-43.  
>[4]	Kummerle R, Grisetti G, Strasdat H, et al. “G2o: A general framework for graph optimization” IEEE International Conference on Robotics and Automation. IEEE, 2011:3607-3613.  
>[5]	Li, X. Rong, and V. P. Jilkov. "Survey of maneuvering target tracking. Part I. Dynamic models." IEEE Transactions on Aerospace &Electronic Systems 39.4(2004):1333-1364.  
>[6]	Fan, X., Fan G., and J. P. Havlicek. "Generative Model for Maneuvering Target Tracking." IEEE Transactions on Aerospace &Electronic Systems 46.2(2010):635-655.    
>[7]	Maybeck, P. S., Worsley, W. H., and Flynn, P. M. (1982) Investigation of constant turn-rate dynamics for airborne vehicle tracking. In Proceedings of the National Aerospace and Electronics Conference (NAECON), 1982, 896–903.  
>[8]	Bryan, R. S. (1980)Cooperative estimation of targets by multple aircraft. M.S. thesis, Air Force Institute of Technology,Wright-Patterson AFB, Ohio, June, 1980.  
>[9]	RAUCH, et al. "Maximum likelihood estimates of linear dynamic systems." Aiaa Student Journal American Institute of Aeronautics & Astronautics 3 i8.3 i8(2012):1445-1450.  
>[10]	W. Li, X. Yang, C. He, S. Ren and Y. Zhang, "Online maximum likelihood filtering for aircraft tracking under low accuracy observations," 2017 Chinese Automation Congress (CAC), Jinan, 2017, pp. 1998-2005.  
>[11]	Xinyan Yan, VadimIndelman, and Byron Boots. "Incremental Sparse GP Regression for Continuous-Time Trajectory Estimation and Mapping." Robotics and Autonomous Systems 87(2017):120-132.  
>[12]	Mercieca, J, P. Aram, and V. Kadirkamanathan. "Unscented Rauch-Tung-Striebel smoothing for nonlinear descriptor systems." Control Conference IEEE, 2015:491-496.  
>[13]	KatjaNummiaro, Esther Koller-Meier, and Luc Van Gool. "An adaptive color-based particle filter." Image and Vision Computing 21.1(2003):99-110.  
>[14]	Zhang, Tianzhu, C. Xu, and M. H. Yang. "Multi-task Correlation Particle Filter for Robust Object Tracking." IEEE Conference on Computer Vision and Pattern Recognition IEEE Computer Society, 2017:4819-4827.  
>[15]	Hol, Jeroen D., T. B. Schon, and F. Gustafsson. "On Resampling Algorithms for Particle Filters." IEEE Nonlinear Statistical Signal Processing Workshop IEEE, 2007:79-82.  
>[16]	Moral, P. Del. "Nonlinear filtering: Interacting particle solution." Markov Processes & Related Fields 2.4(1996):555-580.  
>[17]	Danelljan, Martin, et al. "Discriminative Scale Space Tracking." IEEE Transactions on Pattern Analysis & Machine Intelligence 39.8(2017):1561-1575.  
>[18]	Li Y, Zhu J. “A Scale Adaptive Kernel Correlation Filter Tracker with Feature Integration”. IEEE European Conference on Computer Vision 2014, 8926:254-265.  
>[19]	Arasaratnam, Ienkaran, and S. Haykin. "Cubature Kalman Filters." IEEE Transactions on Automatic Control 54.6(2009):1254-1269.  
>[20]	Nguyen, Trung, et al. "Developing a Cubature Multi-state Constraint Kalman Filter for Visual-Inertial Navigation System." Conference on Computer and Robot Vision IEEE Computer Society, 2017:321-328.  
>[21]	Johansen, Tor A., and T. I. Fossen. "The eXogenous Kalman Filter (XKF)." International Journal of Control 90.2(2017):161-167.  
>[22]	Chen, Badong, et al. "Maximum correntropy Kalman filter." Automatica 76(2017):70-77.  
>[23]	M. Zorzi, "Robust Kalman Filtering Under Model Perturbations," in IEEE Transactions on Automatic Control, vol. 62, no. 6, pp. 2902-2907, June 2017.  
>[24]	Carlone, Luca, et al. "Planar Pose Graph Optimization: Duality, Optimal Solutions, and Verification." IEEE Transactions on Robotics 32.3(2016):545-565.  
>[25]	F. Lu, and E. Milios. "Globally Consistent Range Scan Alignment for Environment Mapping." Autonomous Robots 4.4(1997):333-349.  
>[26]	Kaess, Michael. "Incremental smoothing and mapping." IEEE Transactions on Robotics 24.6(2008):1365-1378.  
>[27]	Hwang, Joo Young, et al. "A fast path planning by path graph optimization." Systems Man & Cybernetics Part A Systems & Humans IEEE Transactions on 33.1(2003):121-129.  
>[28]	Davison, Andrew J, et al. "MonoSLAM: Real-Time Single Camera SLAM." IEEE Transactions on Pattern Analysis and Machine Intelligence 29.6(2007):1052-1067.  
>[29]	Engel, Jakob, T. Schöps, and D. Cremers. "LSD-SLAM: Large-Scale Direct Monocular SLAM." European Conference on Computer Vision8690(2014):834-849.  
>[30]	Mur-Artal, Raúl, J. M. M. Montiel, and J. D. Tardós. "ORB-SLAM: A Versatile and Accurate Monocular SLAM System." IEEE Transactions on Robotics 31.5(2017):1147-1163.  
>[31]	Liu, D. C, and J. Nocedal. "On the limited memory BFGS method for large scale optimization." Mathematical Programming 45.1-3(1989):503-528.  
>[32]	Mandic, D. P. "A generalized normalized gradient descent algorithm." Signal Processing Letters IEEE 11.2(2004):115-118.  
