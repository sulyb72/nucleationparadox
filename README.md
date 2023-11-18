# Nucleation Paradox: Investigating Earth’s Inner Core Formation

Unravel the mystery of Earth's inner core formation with this GitHub repository, dedicated to investigating the detected presence of the inner core. Our research challenges established cooling rates, seismic data, and previous studies, delving into the impact of silicon in the nucleation process.

## Overview

This project utilizes molecular dynamics simulations with LAMMPS to explore iron–silicon mixtures, analogous to prior studies on iron–oxygen alloys. Despite encountering an initial critical undercooling value of 69 ± 69 K, uncertain data impedes conclusive insights. The study proposes future research combining LAMMPS simulations and pre-freezing nuclei methods, spanning various systems like pure iron and iron alloys (oxygen, sulfur, silicon) to decode the inner core nucleation puzzle.

## Methodology

Python is employed for data processing and statistical analysis, enhancing the repository's capability to unravel crucial insights from complex simulations.

## Project Files

### Simulation Data

- **001rdf2.txt and 001rdf3.txt**: Data used in Radial Distribution Functions.py to plot radial distribution functions of atoms.
- **001step2posfinal.txt**: Data used in Atomic Density.py to plot atomic density of atoms at the beginning of a simulation run.
- **3800K-waitingtimes.txt, 4000K-waitingtimes.txt, 4200K-waitingtimes.txt**: Waiting time data used in Waiting Times Histogram.py for simulations at different temperatures.
- **3800K_EOS_data.txt, 4000K_EOS_data.txt, 4200K_EOS_data.txt**: Thermodynamic data used in Equations of State.py for simulations at 3800K, 4000K, and 4200K, respectively.
- **gruneisen_originals.txt**: Data used in Initial Supercooling Estimate.py to plot an existing estimate of the required supercooling.

### Analysis and Visualization Scripts

- **Atomic Density.py**: Script to plot the atomic density of atoms at a specific time position within the simulation.
- **Data Processor.py**: Primary script to process simulated data, extracting key information and generating data plots.
- **Equations of State.py**: Script to calculate and plot thermodynamic equations of state for a group of data simulations at a given temperature.
- **Gibbs Free Energies.py**: Script to plot the relationship between the radius of nuclei and their Gibbs free energy.
- **Initial Supercooling Estimate.py**: Script to estimate the supercooling of the Earth necessary to match geophysical observations.
- **Radial Distribution Functions.py**: Script to plot radial distribution functions for a simulation run of atoms.
- **Waiting Times Histogram.py and Waiting Times.py**: Scripts to plot histograms and average waiting times, respectively, to observe frozen nuclei in the simulation runs.

## Project Report

For a comprehensive overview of the research process and results, refer to [Dissertation.pdf](Dissertation.pdf).

## Supplementary Calculations

- **Main Spreadsheet.xlsx**: Spreadsheet used to perform various geophysical and thermodynamic calculations used in this research.

Feel free to explore, contribute, and join the quest to unravel the nucleation paradox of Earth's inner core!
