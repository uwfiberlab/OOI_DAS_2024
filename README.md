# OOI RCA 2024 DAS experiment
Information for understanding the organization and accessing the data from the Ocean Observatories Initiative (OOI) Regional Cabled Array test of distributed acoustic sensing (DAS) and distributed temperature sensing (DTS).

## Data acquired by OptoDAS
The Alcatel Subsea Networks OptoDAS system was used to acquire data on the South Cable.

The OptoDAS data can be found at http://piweb.ooirsn.uw.edu/das24/data/. This includes 5 directories named by date:
- 200G	20240506
- 617G	20240507
- 617G	20240508
- 617G	20240509
- 382G	20240510

This [repository](https://github.com/uwfiberlab/OOI_DAS_2024) prepared by Qibin Shi and Ethan Williams at University of Washington illustrates how to read and manipulate the HDF5 format data.


## Cable Geometry

![map](docs/OOIcable.jpg)
Map of the OOI RCA cables (black lines), the 15 selected channels for earthquake detection (black squares) and the interrogated segment (red line) between the coast and repeater (red circles). Bathymetric contours are denoted in meters and are unevenly spaced. The full length of cables, seafloor nodes, coastal lines and plate boundaries are shown in the inset map.

The channel locations can be reproduced using [this](https://github.com/uwfiberlab/OOI_DAS_2024/blob/main/ooi24_data_example_and_cable_info.ipynb) code.
