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

- **Historical Global Temperature Data:** Data on the historical global mean surface temperature are publicly available from Berkeley Earth. Available at: https://berkeleyearth.org/data/

  Rohde, R.A., Hausfather, Z., 2020. The Berkeley Earth Land/Ocean Temperature Record. Earth Syst Sci Data 12, 3469–3479. https://doi.org/10.5194/ESSD-12-3469-2020 

# Model Datasets

- **CORDEX Ensemble:** The individual members that comprise the CORDEX ensemble in this analysis are publicly available via the Earth System Grid Federation federated data nodes, e.g., https://esgf-metagrid.cloud.dkrz.de/search.

  Jacob, D., Petersen, J., Eggert, B., Alias, A., Christensen, O.B., Bouwer, L.M., Braun, A., Colette, A., Déqué, M., Georgievski, G., Georgopoulou, E., Gobiet, A., Menut,     L., Nikulin, G., Haensler, A., Hempelmann, N., Jones, C., Keuler, K., Kovats, S., Kröner, N., Kotlarski, S., Kriegsmann, A., Martin, E., van Meijgaard, E., Moseley, C.,      Pfeifer, S., Preuschmann, S., Radermacher, C., Radtke, K., Rechid, D., Rounsevell, M., Samuelsson, P., Somot, S., Soussana, J.F., Teichmann, C., Valentini, R., Vautard,      R.,   Weber, B., Yiou, P., 2014. EURO-CORDEX: New high-resolution climate change projections for European impact research. Reg Environ Change 14, 563–578.
  https://doi.org/10.1007/S10113-013-0499-2/ 

- **ClimEx ensemble:** The 50 realizations (members) that comprise the ClimEx ensemble are publicly available via the Globus file transfer service. Available at: https://www.climex-project.org/data-access/

  Leduc, M., Mailhot, A., Frigon, A., Martel, J.L., Ludwig, R., Brietzke, G.B., Giguère, M., Brissette, F., Turcotte, R., Braun, M., Scinocca, J., 2019. The ClimEx project:   A 50-member ensemble of climate change projections at 12-km resolution over Europe and northeastern North America with the Canadian Regional Climate Model (CRCM5). J Appl    Meteorol Climatol 58, 663–693. https://doi.org/10.1175/JAMC-D-18-0021.1 

- **Global Climate Model (GCM) Global Temperature Data:** The global surface temperature data for the driving global climate models can be accessed via multiple data repositories. Available at (but not exclusively): https://cds.climate.copernicus.eu/datasets/projections-cmip5-monthly-single-levels; https://esgf-metagrid.cloud.dkrz.de/search; https://pcmdi.llnl.gov/mips/cmip5/availability.html
