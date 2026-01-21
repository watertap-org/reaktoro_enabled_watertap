This contains analysis code that reproduces results in "Optimization of Desalination Systems with Detailed Water Chemistry Through Integration of Reaktoro in WaterTAP", currently under review (DOI will be added once published). To generate data presented in the paper execute analysis code in data_generation folder. All analysis sweeps utilize the softening_acid_ro flowsheet. 

1. treatment_sweep.py will run the four waters across a recovery range of 50-90% and use both HCl and H2SO4 as acidification agent and CaO and Na2CO3 as softening agents. Refer to treat treatment_lime_soda_ash_hcl_h2so4_sweep.yaml for configuration details. Data will be saved in data_generation/outputs/treatment_lime_soda_ash_hcl_h2so4_sweep_analysisType_treatment_sweep.h5

2. validation_sweep.py will run the four waters across a recovery range of 50-90% and use either HCl or H2SO4 as acidification agent and only Na2CO3 as softening agent. Refer to treat validation_soda_ash_hcl_h2so4_sweep.yaml for configuration details. Data will be saved in data_generation/outputs/validation_soda_ash_hcl_h2so4_sweep_analysisType_validation_sweep.h5

3. stability_sweep.py will run the four waters across, 13 different hessian configurations, and  a recovery range of 50-90% in 4% steps. The flowsheet will use both HCl and H2SO4 as acidification agent and CaO and Na2CO3 as softening agents. Refer to treat stability_sweep.yaml for configuration details. Data will be saved in data_generation/outputs/stability_sweep_analysisType_stability_sweep.h5

Figure generaton code process each one of the data files generated above:

1. To generate optimal results for base analysis from treatment_sweep, execute:

    a. cost_breakdown_bgw.py -> generates cost breakdown for three brackish water cases. Figures saved in figure_generation/cost_figures

    b. cost_breakdown_hpro.py -> generates cost breakdown for seawater case. Figures saved in figure_generation/cost_figures

    c. plotting_performance_regions.py -> plots cost optimal operation and output metrics for the four waters and their regions. Figures saved in figure_generation/treatment_figures

2. To generate optimal results for validation analysis against Amusat et al. 2024, execute:
    
    a. validation_plotting.py -> will generate comparison plots of costs, chemical dosing, and scaling potentials against Amusat et al. 2024. Figures saved in figure_generation/validation_figures

3. To generate stability analysis results from stability_sweep.py, execute:
    
    a. stability_plotting.py -> will generate figures for IPOPT performance using different hessian options. Figures saved in figure_generation/stability_figures