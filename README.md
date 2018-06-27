# py-cloud-parcel-model
Cloud parcel model written in python designed to be used to simulate experiments in a cloud chamber. 
This model is based on the underlying science of the Aerosol-Cloud-Precipitation-Interactions-Model (http://personalpages.manchester.ac.uk/staff/paul.connolly/research/acpim01.html).

Simple parcel model that calculates pressure, temperature, RH and the growth of aerosol particles into cloud drops and ice crystals. 

How to use:
1. edit namelist.py specifying your initial conditions of aerosol, temperature, pressure and RH
2. run main.py
3. results are plotted in 'output.pdf'

How to use gui:
1. from command line run $ python pyacpim_gui_v0.oy
2. fill in gui form with initial conditions to run model
3. plot output from model run using 'plot' button in gui
![Alt text](./gui_screenshot.png?raw=true "Title")

