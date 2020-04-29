import numpy as np
##############################################################################

# Enter False to process the spectra one-by-one, or enter True to process
# multiple spectra simultaneously over multiple cores. If True, enter number
# of process to run simultaneously. See README for more details. 
multiprocess = False
process_no = 2

# Enter True to display a progress bar for the mcmc sampling. Doesnt work when
# multiprocessing, so enter False.
progress_bar = True


###############################################################################


# The following parameters are suitable for the default CaT region case using 
# the templates provided and LR WEAVE spectra. The following can be tweaked if
# desired, but generally they are satisfactory.

### MCMC SAMPLE PARAMETERS ### 

nwalkers = 200		# number of walkers
burn = 250			# number of burn-in steps
runs = 2000			# number of steps
ntemps = 5			# number of temperatures

### PARAMETER BOUNDARIES ###

# Teff, logg, vsini must be same as template grid
Teff = np.arange(5000, 18500, 500)
logg = np.array([2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
vsini = np.array([0, 50, 100, 150, 200, 250, 300])

drv = np.arange(-500., 500., 0.1)
slopes = np.arange(-0.3e-5, 2.0e-5, 1.0e-10)
intercepts = np.arange(-1.0e-1, 1.0e-1, 1.0e-6)

### TARGET SPECTRUM INFORMATION ###

# Specify wavelength range to analyse
min_wav = 8470
max_wav = 8940

# Enter True to exclude a region within the min max wavelength range, or False
# to not.
exclude_region = False

# ...if True, specify the desired INCLUDED wavelength ranges. More 'regions'
# can be added to further split the spectrum.
region1 = [min_wav, 8610]
region2 = [8630, max_wav]

lines = [region1, region2] # a list of included regions, in order of wavelength 

# pixels-per-resolution element value (PPRE = resolution/wavelength sampling)
PPRE = 5.2

###############################################################################

####################### PROCESSING TEMPLATES (OPTIONAL) #######################

# If you wish to process (crop, rotationally broaden, rebin and smooth) your
# own template set to use with ptmcmc, enter True below (or False if not).
# Templates should be text files including wavelength and unnormalised 
# flux/counts data. The templates are processed and then the flux data put into
# a grid. Processing templates takes a while!

# Note: the parameter boundaries stated above must reflect the grid of your
# templates to be used, including appropriate vsini values if templates are to
# be rotationally broadened. 

process_templates = False

# list of templates to be processed. Listed objects should include paths
template_list_file = 'your_own_templates/template_list'

# index of wavelength and unnormalised flux/counts columns
wav_ind = 0
flux_ind = 2

# directory to write processed templates and template grid to
template_write_directory = 'your_own_templates/processed/'

# wavelengths to crop template spectrum to. Should be sufficiently larger than
# ptmcmc wavelength region to allow for Doppler shifting.
template_crop_min_wav = 6000.	
template_crop_max_wav = 9000.	

# The vsini values if rotational broadening is desired. In this case, the
# provided templates should be unbroadened. The vsini values should be given in
# array form i.e. '[50.0, 100.0, 200.0]'. If rotational broadening is not
# desired, enter 'None'. Remember to enter appropriate vsini values in the
# PARAMETER BOUNDARIES section above.
rotbroads = None#[50.0, 100.0, 150.0, 200.0, 250.0, 300.0]

# The width (i.e., standard deviation) of the Gaussian profile used to
# broaden/smooth. Enter 'None' for no smoothing.
sigma = 0.48

# The desired wavelength sampling. Enter 'None' for no rebinning.  
sampling = 0.25

# The starting and ending indexed of temperature, logg, vsini information in
# the listed templates. E.g. for a listed template called 
# 'C:/Users/me/Documents/templates/t05000g3.00-vsini000.dat', the indices would
# be (33, 38), (39,43), (49, 52) respectively.
temp_ind = (30, 36)
logg_ind = (37, 41)
vsini_ind = (46, 49) 	# if templates are not being rotationally broadened




