# py-cloud-parcel-model
Cloud parcel model written in python designed to be used to simulate experiments in a cloud chamber. 
This model is based on the underlying science of the Aerosol-Cloud-Precipitation-Interactions-Model (http://personalpages.manchester.ac.uk/staff/paul.connolly/research/acpim01.html).

Simple parcel model that calculates pressure, temperature, RH and the growth of aerosol particles into cloud drops and ice crystals. 

How to use:
1. Install anaconda 3
2. Create an ananconda environment using the pyacpim_env.yml file (eg $ conda env create -f pyacpim_env.yml)
3. Enter the new pyacpim environment and run pyACPIM.py (eg $ source activate pyacpim, $ python pyACPIM.py)
