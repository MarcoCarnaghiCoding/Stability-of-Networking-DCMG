# Stability-of-Networking-DCMG

## Stability of Networking DC Microgrids with Active Loads

This repository contains MATLAB source codes related to the research article titled "Stability of Networking DC Microgrids with Active Loads." The work focuses on the stability challenges posed by active loads and communication delays in NMGs.

### Abstract

Direct current networking microgrids (NMGs) enhance the integration of renewable energy sources and promote system resilience through the interconnection of nearby microgrids. However, they face stability challenges influenced by communication delays, fluctuating renewable energy sources, and uncertain loads. This study presents a delay-dependent stability analysis for NMGs, emphasizing robust operation amidst active loads and varying conditions. Circuit simulations are performed to validate insights derived from the stability analysis.

### Repository Structure

```plaintext
├── MatLab_Analysis/
│   ├── analysis_1_Ptot_distvsTau.m
│   ├── analysis_2_aij_PtotvsTau.m
│   ├── analysis_3_Kii_KpivsTau.m
│   └── Analysis_4_AlphavsTau.m
├── Simulink/
│   └── Simulation_alpha_vs_tau.m
│   └── Simulation_Base_4_mg.m
│   └── Simulink_mg_1234_with_DistributedControl_ControlConmutado.slxc
│   └── Simulink_mg_1234_with_DistributedControl_fijo.slxc
│   └── Simulink_Simulation_AlPHA_vs_Tau.slxc
│   └── Simulink_Simulation_resistiveLoad.slxc
└── README.md
```

### Contents

- **MatLab_Analysis/**: Contains MATLAB scripts for the four analyses conducted in the study.
  - `analysis_1.m`: Examines the impact of total load power distribution among clusters.
  - `analysis_2.m`: Analyzes the effects of consensus gains on stability.
  - `analysis_3.m`: Explores tuning of controller gains.
  - `analysis_4.m`: Studies the impact of active and passive load proportions on stability margins.

- **Simulink/**: Contains a Simulink model used for circuit simulations of the NMG scenario.
  

### Getting Started

To run the MATLAB scripts and Simulink model:

1. **Clone the repository**:
   ```bash
   git clone https://github.com/MarcoCarnaghiCoding/Stability-of-Networking-DCMG.git
   cd Stability-of-Networking-DCMG
   ```

2. **Open MATLAB** and navigate to the `MatLab_Analysis` folder to run the analysis scripts. You can execute each script individually to observe the results.

3. **Open Simulink** and load the `NMG_Simulation.slx` file to visualize the circuit simulation.

### Requirements

- MATLAB with the necessary toolboxes (e.g., Control System Toolbox, Simulink)
- CVX Research library for MATLAB for convex formulation analyses

### References

+ [1] M. Shafiullah, A. Refat, M. E. Haque, D. M. H. Chowdhury, S. Hossain, A. Alharbi, S. Alam, A. Ali, and S. Hossain, “Review of recent developments in microgrid energy management strategies,”
Sustainability, vol. 14, pp. 1–30, 11 2022.
+ [2] Q. Zhou, M. Shahidehpour, A. Paaso, S. Bahramirad, A. Alabdulwahab, and A. Abusorrah, “Distributed control and communication strategies in networked microgrids,” IEEE Communications Surveys Tutorials, vol. 22, no. 4, pp. 2586–2633, 2020.
+ [3] J. Liu, W. Zhang, and G. Rizzoni, “Robust stability analysis of dc microgrids with constant power loads,” IEEE Trans. Power Syst., vol. 33, no. 1, pp. 851–860, 2018.
+ [4] M. Jeeninga, C. De Persis, and A. van der Schaft, “Dc power grids with constant-power loads—part i: A full characterization of power flow feasibility, long-term voltage stability, and their correspondence,” IEEE Transactions on Automatic Control, vol. 68, no. 1, pp. 2–17, 2023.
+ [5] Z. Liu, M. Su, Y. Sun, X. Zhang, X. Liang, and M. Zheng, “A comprehensive study on the existence and stability of equilibria of dc-distribution networks with constant power loads,” IEEE Transactions
on Automatic Control, vol. 67, no. 4, pp. 1988–1995, 2022.
+ [6] J. Liu, Y. Zhang, A. J. Conejo, and F. Qiu, “Ensuring transient stability with guaranteed region of attraction in dc microgrids,” IEEE Transactions on Power Systems, vol. 38, no. 1, pp. 681–691, 2023.
+ [7] Z. Liu, M. Su, Y. Sun, W. Yuan, H. Han, and J. Feng, “Existence and stability of equilibrium of dc microgrid with constant power loads,” IEEE Trans. Power Syst., vol. 33, no. 6, pp. 6999–7010, 2018.
+ [8] M. Carnaghi, P. Cervellini, M. Judewicz, R. Garcia Retegui, and M. Funes, “Stability analysis of a networking dc microgrid with distributed droop control and cpls,” IEEE Latin America Transactions,
vol. 21, p. 966–975, Sep. 2023.
+ [9] M. Dong, L. Li, Y. Nie, D. Song, and J. Yang, “Stability analysis of a novel distributed secondary control considering communication delay in dc microgrids,” IEEE Transactions on Smart Grid, vol. 10, no. 6,
pp. 6690–6700, 2019.
+ [10] A. B. Shyam, S. Anand, and S. R. Sahoo, “Effect of communication delay on consensus-based secondary controllers in dc microgrid,” IEEE Trans. Ind. Electron., vol. 68, no. 4, pp. 3202–3212, 2021.
+ [11] W. Yao, Y. Wang, Y. Xu, and C. Dong, “Small-signal stability analysis and lead-lag compensation control for dc networked-microgrid under multiple time delays,” IEEE Trans. Power Syst., vol. 38, no. 1,
pp. 921–933, 2023.
+ [12] Z. Liu, M. Su, Y. Sun, H. Han, X. Hou, and J. M. Guerrero, “Stability analysis of dc microgrids with constant power load under distributed control methods,” Automatica, vol. 90, p. 62–72, Apr. 2018.
+ [13] M. Carnaghi, P. Cervellini, R. G. Retegui, M. Judewicz, and M. Funes, “Analysis of constant power loads impact on dc microgrid with distributed control,” in 2022 IEEE Biennial Congress of Argentina
(ARGENCON), pp. 1–7, 2022.
+ [14] D. Liu, K. Jiang, L. Yan, X. Ji, K. Cao, and P. Xiong, “A fully distributed economic dispatch method in dc microgrid based on consensus algorithm,” IEEE Access, vol. 10, pp. 119345–119356, 2022.
+ [15] W. Xie, M. Han, W. Cao, J. M. Guerrero, and J. C. Vasquez, “System-level large-signal stability analysis of droop-controlled dc microgrids,” IEEE Trans. Power Electron., vol. 36, no. 4, pp. 4224–4236 2021.
+ [16] S. Liu, X. Li, M. Xia, Q. Qin, and X. Liu, “Takagi-sugeno multimodeling-based large signal stability analysis of dc microgrid clusters,” IEEE Trans. Power Electron., vol. 36, no. 11, pp. 12670–12684, 2021.
+ [17] R. Han, M. Tucci, A. Martinelli, J. M. Guerrero, and G. Ferrari-Trecate, “Stability analysis of primary plug-and-play and secondary leader-based controllers for dc microgrid clusters,” IEEE Trans. Power Syst., vol. 34, no. 3, pp. 1780–1800, 2019.
+ [18] Y. Yu, G.-P. Liu, and W. Hu, “Coordinated distributed predictive control for voltage regulation of dc microgrids with communication delays and data loss,” IEEE Transactions on Smart Grid, vol. 14, no. 3,
pp. 1708–1722, 2023.
+ [19] Y. He, Q.-G. Wang, L. Xie, and C. Lin, “Further improvement of free-weighting matrices technique for systems with time-varying delay,” IEEE Trans. Autom. Control, vol. 52, no. 2, pp. 293–299, 2007.
+ [20] C. Li, F. de Bosio, S. K. Chaudhary, M. Graells, J. C. Vasquez, and J. M. Guerrero, “Operation cost minimization of droop-controlled dc microgrids based on real-time pricing and optimal power flow,” in IECON 2015 - 41st Annual Conference of the IEEE Industrial Electronics Society, pp. 003905–003909, 2015.
+ [21] M. Grant and S. Boyd in Graph implementations for nonsmooth convex programs, pp. 95–110, Springer-Verlag Limited, 2008.
+ [22] M. Zaery and M. A. Abido, “Distributed optimal power dispatch for islanded dc microgrids with time delays,” IEEE Access, vol. 12, pp. 12533–12544, 2024.


### Acknowledgments

This work was supported in part by the Universidad Nacional de Mar del Plata (UNMDP), the Consejo Nacional de Investigaciones Científicas y Tecnológicas (CONICET), Argentina.


### Contact Authors

{mcarnaghi,paulacervellini,marcosj,rgarcia,mfunes}@fi.mdp.edu.ar}


