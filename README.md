**Mathematical models of the VTA-NAc-mPFC neural circuit under different dopamine concentration inputs**  
**Please cite these references:**  
1. Li Y, Zhang B, Pan X, et al. Dopamine-mediated major depressive disorder in the neural circuit of ventral tegmental area-nucleus accumbens-medial prefrontal cortex: from biological evidence to computational models[J]. Frontiers in Cellular Neuroscience, 2022, 16: 923039.  
2. Li Y, Zhang B, Liu Z, et al. Neural energy computations based on Hodgkin-Huxley models bridge abnormal neuronal activities and energy consumption patterns of major depressive disorder[J]. Computers in Biology and Medicine, 2023, 166: 107500.
3. Li Y, Zhang B, Liu J, Wang R. Neural energy coding patterns of dopaminergic neural microcircuit and its impairment in major depressive disorder: A computational study. PLoS Comput Biol. 2025 Apr 7;21(4):e1012961.
  
**Matlab was used for modeling.  
Matlab: R2021b**  
  
**1. Neurodynamical models**  
Run test_FiringSequence_Dopamineratio_Change_2.m to calculate the differential equation models for Low, Medium, High, and Full dopamine concentrations.  
The results of the differential equation models will be a matrix with 150001 rows and 1527 columns.  
  
**2. Neural energy models**  
Run CalculateAllCurrents.mlx to get the ion channel currents, synaptic currents, compartment currents, and stimulus currents.  
Run CalculatePower_Expereiment_1.mlx to get the neural power.  
Run CalculateSurfaceArea.mlx to get the surface area of each neuron.  
  
  
**Additional notes:**  
1. We have adjusted the parameter settings for different studies. Please adjust the parameters of this DEMO according to the references if you want to reproduce our results.  
2. Since the results files are too large, we are unable to upload the complete files to GITHUB. Please follow the steps above, and run the scripts based on the parameters described in our papers to reproduce our results.
3. Please contact Dr. Yuanxi Li at dr.yuanxili@gmail.com if you have any questions.  
