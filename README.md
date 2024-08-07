# OOI RCA 2024 DAS Experiment

This repository provides information on how to understand the organization and access the data from the Ocean Observatories Initiative (OOI) Regional Cabled Array's Distributed Acoustic Sensing (DAS) and Distributed Temperature Sensing (DTS) tests conducted in 2024.

## Overview

The data for this experiment was acquired using the OptoDAS system from Alcatel Subsea Networks. The focus was on the South Cable of the OOI RCA.

## Data Access

The OptoDAS data can be accessed via the following link:

[OptoDAS Data Repository](http://piweb.ooirsn.uw.edu/das24/data/)

The data is organized into directories named by date. Here’s a breakdown:

- **200G**: 2024-05-06
- **617G**: 2024-05-07
- **617G**: 2024-05-08
- **617G**: 2024-05-09
- **382G**: 2024-05-10

For details on how to read and manipulate the HDF5 format data, refer to the [example code in this repository](https://github.com/uwfiberlab/OOI_DAS_2024) prepared by Qibin Shi and Ethan Williams from the University of Washington.

## Cable Geometry

The following map shows the OOI RCA cables (black lines) and the interrogated segment (red line) between the coast and the repeater (red circles). Bathymetric contours are shown in meters and are unevenly spaced. The inset map provides additional context including cable lengths, seafloor nodes, coastal lines, and plate boundaries.

![Cable Map](docs/OOIcable.jpg)

#### Reproducing Channel Locations

You can reproduce the channel locations using the [example code](https://github.com/uwfiberlab/OOI_DAS_2024/blob/main/ooi24_data_example_and_cable_info.ipynb) available in this repository. This Jupyter notebook provides a step-by-step guide to visualize and analyze the cable geometry.
