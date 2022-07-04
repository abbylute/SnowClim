# SnowClim

Seasonal snowpack dynamics shape the biophysical and societal characteristics of many global regions. However, snowpack accumulation and duration have generally declined in recent decades largely due to anthropogenic climate change. Mechanistic understanding of snowpack spatiotemporal heterogeneity and climate change impacts will benefit from snow data products that are based on physical principles, that are simulated at high spatial resolution, and that cover large geographic domains. Most existing datasets do not meet these requirements, hindering our ability to understand both contemporary and changing snow regimes and to develop adaptation strategies in regions where snowpack patterns and processes are important components of Earth systems.
We developed a computationally efficient process-based snow model, SnowClim, that can be run in the cloud. The model was evaluated and calibrated at Snowpack Telemetry sites across the western United States (US), achieving a site-median root mean square error for daily snow water equivalent of 64 mm, bias in peak snow water equivalent of -2.6 mm, and bias in snow duration of -4.5 days when run hourly. Positive biases were found at sites with mean winter temperature above freezing where the estimation of precipitation phase is prone to errors. The model was applied to the western US (a domain covering 3.1 million km2) using newly developed forcing data created by statistically downscaling pre-industrial, historical, and pseudo-global warming climate data from the Weather Research and Forecasting (WRF) model. The resulting product is the SnowClim dataset, a suite of summary climate and snow metrics, including monthly snow water equivalent (SWE) and snow depth, as well as annual maximum SWE and snow cover duration, for the western US at 210 m spatial resolution (Lute et al., 2021). The physical basis, large extent, and high spatial resolution of this dataset enable novel analyses of changing hydroclimate and its implications for natural and human systems. This repository contains the code for the SnowClim model.

Additional details regarding the SnowClim model physics, model calibration, climate data downscaling, model application to the western US, and model performance are available in:
Lute, A. C., Abatzoglou, J., and Link, T.: SnowClim v1.0: high-resolution snow model and data for the western United States, Geosci. Model Dev., 15, 5045–5071, https://doi.org/10.5194/gmd-15-5045-2022, 2022.

The western US snow and climate datasets are available here: https://www.hydroshare.org/resource/acc4f39ad6924a78811750043d59e5d0/

The SnowClim Model code is also available here:

[![DOI](https://zenodo.org/badge/475599638.svg)](https://zenodo.org/badge/latestdoi/475599638)

