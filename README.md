# EnKF-and-PF-DA-for-UEB-SACSMA-rutpix7 

Data assimilation to three coupled (snowmelt, soil moisture accounting, and surface and streamflow routing) models. 

- The Utah Energy Balance (UEB) snowmelt model with snow data asssimilation using the Ensemble Kalman Filter (EnKF)  

- The Sacramento Soil Moisture Accounting Model (SACSMA) and  

- rutpix7 model for surface and streafflow routing with watershed outlet discharge data assimilation using the Particle Filter (PF) code  

These are based on National Weather Service's Office of Hydrologic Development (OHD) HL-RDHM (Hydrology Laboratory--Research Distributed Hydrologic Model, https://www.cbrfc.noaa.gov/present/rdhm/RDHM_User_Manual.pdf.   

The code here first implements the UEB model as component of RDHM and SACSMA and rutpix7 are reconfigured to run on GPUs.   

Programming languages: C++, CUDA   Requires: MPI, GPU (CUDA), NetCDF  

For more see the paper: Ensemble Streamflow Forecasting Using an Energy Balance Snowmelt Model Coupled to a Distributed Hydrologic Model with Assimilation of Snow and Streamflow Observations (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019WR025472)
