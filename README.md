# CCDetectonAttribution_Esbjerg

**Climate change detection and atrribution analysis of the extreme rainfall event on the 27th of September 2024 in Esbjerg, Denmark, with a focus on:**

- Bias correction of climate model data
- Performing sampling method (Rx1day) and linking to global mean surface temperature (GMST) anomaly
- **Employing the standardized protocol of the World Weather Attribution as detailed by Philip et al. (2020) and Otto et al. (2024)**

This repository contains **jupyter notebooks (for the preparation of observations and model data)** and **R scripts (for the detection and attribution of the extreme rainfall event)**.

---

# 🌍 Data Overview

The observation and model datasets used in this analysis are publicly available. A brief overview of the datasets used in the analysis (and where to access them) is provided here. 

# Observational Datasets

- **Era5-Land Reanalysis:** Precipitation data from Era5-Land are publicly available from the Copernicus Climate Data Store. Available at: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land

  Copernicus Climate Change Service (C3S), 2022a. ERA5-Land hourly data from 1950 to present. [Dataset]. Copernicus Climate Change Service (C3S) Climate Data Store (CDS).      https://doi.org/10.24381/cds.e2161bac 

  Muñoz-Sabater, J., Dutra, E., Agustí-Panareda, A., Albergel, C., Arduini, G., Balsamo, G., Boussetta, S., Choulga, M., Harrigan, S., Hersbach, H., Martens, B., Miralles,     D.G., Piles, M., Rodríguez-Fernández, N.J., Zsoter, E., Buontempo, C., Thépaut, J.-N., 2021. ERA5-Land: A state-of-the-art global reanalysis dataset for land                 applications.   Earth Syst Sci Data 13, 4349–4383. https://doi.org/https://doi.org/10.5194/essd-13-4349-2021

- **Copernicus European Regional Reanalysis (CERRA-) Reanalysis:** Precipitation data from CERRA are publicly available from the Copernicus Climate Data Store. Available at: https://cds.climate.copernicus.eu/datasets/reanalysis-cerra-single-levels

  Copernicus Climate Change Service (C3S), 2022b. CERRA sub-daily regional reanalysis data for Europe on single levels from 1984 to present. [Dataset]. Copernicus Climate      Change Service (C3S) Climate Data Store (CDS). https://doi.org/10.24381/cds.622a565a 

  Ridal, M., Bazile, E., Le Moigne, P., Randriamampianina, R., Schimanke, S., Andrae, U., Berggren, L., Brousseau, P., Dahlgren, P., Edvinsson, L., El-Said, A., Glinton, M.,   Hagelin, S., Hopsch, S., Isaksson, L., Medeiros, P., Olsson, E., Unden, P., Wang, Z.Q., 2024. CERRA, the Copernicus European Regional Reanalysis system. Quarterly Journal    of the Royal Meteorological Society 150, 3385–3411. https://doi.org/10.1002/QJ.4764 

- **Station Data:** Precipitation data from the Danish Meteorological Institute are publicly availale from a technical report and DMI's API. Precipitation data used in the analysis are sourced from two stations. The data for Station #6088 (1872-2020) is available at: https://www.dmi.dk/fileadmin/Rapporter/2021/DMIRep21-02.zip. The data from Station #5340 (2021-2025) can be accessed via DMI's API at: https://opendatadocs.dmi.govcloud.dk/en/Data/Meteorological_Observation_Data.

- **Global Temperature Data:** Data on the historical global mean surface temperature are publicly available from Berkeley Earth. Available at: https://berkeleyearth.org/data/

  Rohde, R.A., Hausfather, Z., 2020. The Berkeley Earth Land/Ocean Temperature Record. Earth Syst Sci Data 12, 3469–3479. https://doi.org/10.5194/ESSD-12-3469-2020 


