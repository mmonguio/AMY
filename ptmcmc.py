import numpy as np
from scipy import interpolate
import matplotlib
from datetime import datetime
from scipy.stats.mstats import mquantiles
import os
from astropy.io import fits
import glob
from parameters import *
import sys
import argparse
from astropy.table import Table
##sys.path.insert(0, 'python')
# matplotlib.rc('text', usetex=True)
# matplotlib.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}",
#                                    r"\usepackage{color}"]
# matplotlib.use('Agg')
from ptemcee import Sampler as PTSampler
import matplotlib.pyplot as plt
##import corner
startTime = datetime.now()


parser = argparse.ArgumentParser(description="Calculate radial velocities and stellar parameters of WEAVE target spectra.")

parser.add_argument("--infile", type=str, required=True, help="input file", nargs=1)
parser.add_argument("--outdir", type=str, required=True, help="output directory", nargs=1)
parser.add_argument("--targlist", type=str, required=True, help="path to list of FIBREIDs or TARGIDs to be analysed. To analyse all BA stars, enter 'all'.", nargs=1)
parser.add_argument("--params", type=str, required=False, default='parameters.py', help="path to parameter file", nargs=1)

#parser.add_argument("--setups", type=str, default=None, required=True, help="input setups", nargs='*')

args = parser.parse_args()
write_directory = args.outdir[0]
target_list = args.targlist[0]
data_file = args.infile[0]


# processing (cropping, smoothing, rebinning, rotational broadening) templates if required 

if process_templates==True:
	print('Beginning processing of templates')
	from PyAstronomy import pyasl
	import scipy.stats
	templatelist = np.genfromtxt(template_list_file, dtype=None, encoding=None)

	# Ensuring the number of templates to be processed matches the number of points on the grid. 
	notemps = len(templatelist)
	d_size = len(Teff) * len(logg) * len(vsini)
	if notemps < d_size:
		print('ERROR: The amount of templates to be processed is less than the amount of points in the grid. Add more templates to the list or reduce the grid as necessary (see PARAMETER BOUNDARIES in parameters file).')
		sys.exit()
	if notemps > d_size: 
		print('ERROR: The amount of templates to be processed is greater than the amount of points in the grid. Remove templates from the list or expand the grid as necessary (see PARAMETER BOUNDARIES in parameters file).')
		sys.exit()

	def restrict_range(w, f):		# to crop the templates
		our_range = (w > template_crop_min_wav) & (w < template_crop_max_wav)
		w, f = w[our_range], f[our_range]
		return w, f

	def smooth(w, f, sig):		# to smooth/broaden templates to match resolution of observed spectrum
		f = pyasl.broadGaussFast(w, f, sig)
		return f

	def rotbroad(w, f, vsini):		# to rotationally broaden the templates
		f = pyasl.rotBroad(w, f, 0.6, float(vsini))
		return f

	def rebin(w, f, samp):		# to rebin the templates to match the sampling of the observed spectrum
		f, bin_edges, binnumber = scipy.stats.binned_statistic(w, f, statistic = 'mean', bins = (w[-1] - w[0]) / samp)
		bin_width = bin_edges[1] - bin_edges[0]
		w = bin_edges[1:] - bin_width/2
		return w, f

	def process_template(w, f, sig=sigma, samp=sampling):
		if sigma != None:
			f = smooth(w, f, sig)
		if sampling != None:
			w, f = rebin(w, f, samp)
		file = open(template_write_directory + '/'+ os.path.splitext(os.path.basename(i))[0] + '_processed', "w")
		for index in range(len(w)):
			file.write(str(w[index]) + " " + str(f[index]) + "\n")
		file.close()
		return w, f

	def process_template_vsini(w, f, vsini, sig=sigma, samp=sampling):
		f = rotbroad(w, f, vsini)
		if sigma != None:
			f = smooth(w, f, sig)
		if sampling != None:
			w, f = rebin(w, f, samp)
		file = open(template_write_directory + '/' + os.path.splitext(os.path.basename(i))[0] + '_processed_vsini' + str(int(vsini)), "w")
		for index in range(len(w)):
			file.write(str(w[index]) + " " + str(f[index]) + "\n")
		file.close()
		return w, f

	def write_to_grid(templatename, f, t_ind=temp_ind, l_ind=logg_ind, v_ind=vsini_ind, vsini0=None):
		# finding the teff and logg information from the template name, to write to corresponding grid point
		t1 = float(templatename[temp_ind[0]:temp_ind[1]])
		l1 = float(templatename[logg_ind[0]:logg_ind[1]])
		if vsini0==None:
			v1 = float(templatename[vsini_ind[0]:vsini_ind[1]])
		else: 
			v1 = float(vsini0)

		print('Teff: ' + str(t1) + ', logg: ' + str(l1) + ', vsini: ' + str(v1))

		d[np.where(Teff == t1)[0][0], np.where(logg == l1)[0][0], np.where(vsini == v1)[0][0]] = f

	t_wavelength_0, t_flux_0 = np.genfromtxt(templatelist[0], unpack=True, usecols=(wav_ind, flux_ind), dtype=None, encoding=None)
	t_wavelength_0, t_flux_0 = restrict_range(t_wavelength_0, t_flux_0)
	if sampling != None:
		t_wavelength_0, t_flux_0 = rebin(t_wavelength_0, t_flux_0, sampling)

	d = np.zeros((len(Teff), len(logg), len(vsini), len(t_wavelength_0)))
	d_filled = 0

	for i in templatelist:
		print('Processing: ' + i)
		t_wavelength, t_flux = np.genfromtxt(i, unpack=True, usecols=(wav_ind, flux_ind), dtype=None, encoding=None)
		t_wavelength, t_flux = restrict_range(t_wavelength, t_flux)
		if rotbroads == None:
			_, t_flux = process_template(t_wavelength, t_flux, sigma, sampling)
			print('template processed, writing to grid at:')
			write_to_grid(i, t_flux)
			d_filled += 1
			print('grid is ' + str("{0:.1f}".format(float(d_filled)/float(d_size) * 100.)) + ' percent full')



		else:
			print('processing for vsini: 0')
			_, t_flux_v0 = process_template(t_wavelength, t_flux, sig=sigma, samp=sampling)
			print('template processed, writing to grid at:')
			write_to_grid(i, t_flux_v0, vsini0=0.)
			d_filled += 1
			print('grid is ' + str("{0:.1f}".format(float(d_filled)/float(d_size) * 100.)) + ' percent full')

			for j in rotbroads:
				if j == 0:
					continue
				print('processing for vsini: ' + str(j))
				_, t_flux_vj = process_template_vsini(t_wavelength, t_flux, j, sig=sigma, samp=sampling)
				print('template processed, writing to grid at:')
				write_to_grid(i, t_flux_vj, vsini0=j)
				d_filled += 1
				print('grid is ' + str("{0:.1f}".format(float(d_filled)/float(d_size) * 100.)) + ' percent full')


	np.save(template_write_directory + 'template_grid.npy', d)


# creating /plots and /results folders in write_directory if they dont already exist 
if not os.path.exists(write_directory + '/plots'):
    os.mkdir(write_directory + '/plots')
    print(write_directory + 'plots folder created')
if not os.path.exists(write_directory + '/results'):
    os.mkdir(write_directory + '/results')
    print(write_directory + 'results folder created')

if process_templates==False:
	# location of folder containing this script
	templatedirectory = os.path.dirname(os.path.abspath(__file__)) + '/templates/'
	# loading the template flux grid
	flux_data_all = np.load(templatedirectory + 'template_grid.npy')
	#loading the wavelength data
	templatewavelength = np.loadtxt(templatedirectory + 'wavelength_data.dat')
else:
	# the first template for the template wavelength data
	flux_data_all = np.load(template_write_directory + 'template_grid.npy')
	templatewavelength = t_wavelength_0

template_min_wav = min_wav*(1.0 + (min(drv)/299792.458)) - 10.
template_max_wav = max_wav*(1.0 + (max(drv)/299792.458)) + 10.

templatemask = (templatewavelength > template_min_wav) & (templatewavelength < template_max_wav)
flux_data = np.zeros((len(Teff), len(logg), len(vsini), len(templatewavelength[templatemask])))
templatewavelength = templatewavelength[templatemask]

for ii in range(len(Teff)):
	for jj in range(len(logg)):
		for kk in range(len(vsini)):
			flux_data[ii,jj,kk] = flux_data_all[ii,jj,kk][templatemask]


ipo = interpolate.RegularGridInterpolator((Teff, logg, vsini), flux_data, method='linear')

# grid of parameters for walker initial positions
ini_grid = [Teff, logg, vsini, drv, slopes, intercepts]
ndim=6

# define edges of parameter space
teffmin, teffmax, loggmin, loggmax, vsinimin, vsinimax, rvmin, rvmax, slopemin, slopemax, interceptmin, interceptmax = min(Teff), max(Teff), min(logg), max(logg), min(vsini), max(vsini), min(drv), max(drv), min(slopes), max(slopes), min(intercepts), max(intercepts)

def model(X, wavelength):
	i, j, k, l, m, n = X 
	# interpolating template grid with teff (i), logg (j), vsini (k) trial parameter
	templateflux = ipo([i, j, k])[0]
	# interpolating on wavelength axis for trial RV (l)
	fi = interpolate.interp1d(templatewavelength*(1.0 + l/299792.458), templateflux)
	return fi(wavelength)

def lnprior(X):
	i, j, k, l, m, n = X
	# flat prior, edges should corespond to template grid
	if (teffmin <= i <= teffmax) & (loggmin <= j <= loggmax) & (vsinimin <= k <= vsinimax) & (rvmin <= l <= rvmax) & (slopemin <= m <= slopemax) & (interceptmin <= n <= interceptmax):
		return 0.0
	else:
		return -np.inf

def lnlike(X, wavelength, flux, noisespec, mask):
	i, j, k, l, m, n = X
	z = m, n 
	f = np.poly1d(z)
	if exclude_region == False:
		return -(np.sum((flux - (model(X, wavelength)*f(wavelength)))**2/(2*PPRE*noisespec**2)))
	else:
		return -(np.sum((flux[mask] - (model(X, wavelength)*f(wavelength))[mask])**2/(2*PPRE*noisespec[mask]**2)))

def mcmc_one(t):

	print("Processing: " + t)  
	targ_start = datetime.now()

	# checking if spectrum file exists, and if result already written
	if os.path.exists(write_directory + '/results/' + os.path.splitext(os.path.basename(t))[0] + '_results'):
		print('WARNING: result for '+t+' already exists, skipping')
		return

	# specify unnormalised calibrated flux
	with fits.open(data_file) as ALLDATA:
		final_spectra = ALLDATA[1].data[info['TARGID'] == t][0]
		spectra_before_sky_subtraction = ALLDATA[3].data[info['TARGID'] == t][0]
		Calibration_function = ALLDATA[5].data[info['TARGID'] == t][0]
		try: 
			idx=list(ALLDATA[6].data['FIBREID']).index(t)
		except:
			idx=list(ALLDATA[6].data['TARGID']).index(t)
		NSPEC = ALLDATA[6].data['Nspec'][idx]
		FIBREID = ALLDATA[6].data['FIBREID'][idx]
		CNAME = ALLDATA[6].data['CNAME'][idx]
	flux = final_spectra * Calibration_function * 1.0e18

	# restrict to desired wavelength range
	flux = flux[targetmask]

	# create corresponsing noise spectrum
	noisespec = np.sqrt((2.*spectra_before_sky_subtraction - final_spectra)*Calibration_function*1.0e18)
	noisespec = noisespec[targetmask]

	# mask for excluding a wavelength region
	if exclude_region == True:
		mask = np.zeros(len(wavelength), dtype=bool)
		for line in lines:
			mask |= (wavelength >= line[0]) & (wavelength <= line[1])
	else:
		mask = np.ones(len(wavelength), dtype=bool)

	# choose initial walker positions
	pos = [[[np.random.choice(i) for i in ini_grid] for i in range(nwalkers)] for i in range(ntemps)]

	# initialise MCMC sampler
	sampler = PTSampler(ntemps=ntemps, nwalkers=nwalkers, dim=ndim, logl=lnlike, logp=lnprior, Tmax=np.inf, loglargs=[wavelength, flux, noisespec, mask])

	# run MCMC sampler for burn period
	if progress_bar==True:
		print("running burn")
	pos, prob, state = sampler.run_mcmc(pos, burn, adapt=True, progress_bar=progress_bar)
	# reset sampler, run MCMC sampler for run period with walkers starting at their positions at the end of burn
	sampler.reset()
	if progress_bar==True:
		print("running runs")
	sampler.run_mcmc(pos, runs, adapt=True, progress_bar=progress_bar)
	samples=sampler.chain[0, :, :, :].reshape((-1, ndim))

	# plot walker paths
	ylabels = ['Teff', 'logg', 'vsini', 'RV', 'slope', 'intercept']
	for m in range(ndim):
		plt.subplot(ndim,1,m+1)
		plt.plot(sampler.chain[0,:,:,m].transpose(), alpha=0.2)
		plt.ylabel(ylabels[m])
	plt.xlabel('Step')
	plt.savefig(write_directory + '/plots/' + t + '_walkers.png', bbox_inches='tight')
	plt.close()

	# calculate 16th, 50th, 84th quantiles of the parameter samples
	quantiles = mquantiles(samples, prob=[0.16, 0.50, 0.84], axis=0)

	targetname = os.path.splitext(os.path.basename(t))[0]
	acceptance_r = np.mean(sampler.acceptance_fraction)

	# print acceptance fraction, should be between 0.2-0.5 for efficient sampling
	print("Mean acceptance fraction: {0:.3f}"
                .format(acceptance_r))

	# The parameter results
	Teff_r = quantiles[1][0]
	Teffminus_r = quantiles[1][0] - quantiles[0][0]
	Teffplus_r = quantiles[2][0] - quantiles[1][0]
	logg_r = quantiles[1][1]
	loggminus_r = quantiles[1][1] - quantiles[0][1]
	loggplus_r = quantiles[2][1] - quantiles[1][1]
	vsini_r = quantiles[1][2]
	vsiniminus_r = quantiles[1][2] - quantiles[0][2]
	vsiniplus_r = quantiles[2][2] - quantiles[1][2]
	RV_r = quantiles[1][3]
	RVminus_r = quantiles[1][3] - quantiles[0][3]
	RVplus_r = quantiles[2][3] - quantiles[1][3]
	slope_r = quantiles[1][4]
	slopeminus_r = quantiles[1][4] - quantiles[0][4]
	slopeplus_r = quantiles[2][4] - quantiles[1][4]
	intercept_r = quantiles[1][5]
	interceptminus_r = quantiles[1][5] - quantiles[0][5]
	interceptplus_r = quantiles[2][5] - quantiles[1][5]

	#fig = corner.corner(samples, quantiles=[0.16, 0.50, 0.84], labels=['Teff', 'log(g)', 'vsini', 'RV', 'slope', 'intercept'], show_titles=True, title_kwargs={"fontsize": 10}, plot_datapoints=True, plot_contours=True, auto_bars=True, data_kwargs={"alpha": 0.005})
	# fig.savefig(write_directory + t + '_cornerplot.png', bbox_inches='tight')
	# plt.close()

	# parameters of best fit (for plotting)
	Xp = Teff_r, logg_r, vsini_r, RV_r, slope_r, intercept_r
	fp = np.poly1d(Xp[-2:])
	fitp = fp(wavelength)

	# plot spectrum and best-fit
	plt.plot(wavelength, flux, wavelength, model(Xp, wavelength) * fitp)
	if exclude_region == True:
		for n,line in enumerate(lines):
			if n == 0:
				plt.vlines([line[1]], flux.min()*0.8, flux.max()*1.2, colors='r', linestyles='dashed')
			elif n == (len(lines)-1):
				plt.vlines([line[0]], flux.min()*0.8, flux.max()*1.2, colors='r', linestyles='dashed')
			else:
				plt.vlines([line[0], line[1]], flux.min()*0.8, flux.max()*1.2, colors='r', linestyles='dashed')

	plt.xlim(min_wav, max_wav)
	plt.ylim(flux.min()*0.8, flux.max()*1.2)
	plt.xlabel(r'Wavelength ($\AA$)')
	plt.ylabel('Calibrated counts')
	plt.savefig(write_directory + '/plots/' + t + '_spectrum.png', bbox_inches='tight')
	plt.close()

	# plot mapping function
	plt.plot(wavelength, fitp, wavelength, flux / model(Xp, wavelength))
	plt.xlim(min_wav, max_wav)
	plt.xlabel(r'Wavelength ($\AA$)')
	plt.ylabel('Calibrated counts')
	plt.savefig(write_directory + '/plots/' + t + '_mapping_function.png', bbox_inches='tight')
	plt.close()

	# write results
	tab = open(write_directory + '/results/' + t +  '_results', "w")
	tab.write(np.str(NSPEC) + " " + np.str(FIBREID) + " " + np.str(CNAME) + " " + t + " " + np.str(acceptance_r) + " " + np.str(Teff_r) + " " + np.str(Teffminus_r) + " " + np.str(Teffplus_r) + " " + np.str(logg_r) + " " + np.str(loggminus_r) + " " + np.str(loggplus_r) + " " + np.str(vsini_r) + " " + np.str(vsiniminus_r) + " " + np.str(vsiniplus_r) + " " + np.str(RV_r) + " " + np.str(RVminus_r) + " " + np.str(RVplus_r) + " " + np.str(slope_r) + " " + np.str(slopeminus_r) + " " + np.str(slopeplus_r) + " " + np.str(intercept_r) + " " + np.str(interceptminus_r) + " " + np.str(interceptplus_r) + "\n")
	tab.close()

	targ_end = datetime.now() - targ_start
	print(targ_end)

def make_output_fits():
	output_files = glob.glob(write_directory + 'results/*_results')
	all_res = []
	for i in output_files:
		with open(i) as outf:
			all_res.append(outf.read().split())
	t = Table(rows=all_res, names=('NSPEC', 'FIBREID', 'CNAME', 'TARGID', 'Acceptance', 'Teff', 'Teff_minus', 'Teff_plus', 'logg', 'logg_minus', 'logg_plus', 'vsini', 'vsini_minus', 'vsini_plus', 'RV', 'RV_minus', 'RV_plus', 'slope', 'slope_minus', 'slope_plus', 'intercept', 'intercept_minus', 'intercept_plus'))
	t.write(write_directory+'results/'+os.path.splitext(os.path.basename(args.infile[0]))[0]+'_ptmcmc.fits', format='fits', overwrite=False)




if __name__ ==  '__main__':
	# checking if results table already exists
	if os.path.exists(write_directory+'results/'+os.path.splitext(os.path.basename(args.infile[0]))[0]+'_ptmcmc.fits'):
		print('WARNING: result table already exists, ending process.')
		sys.exit()
	if target_list == 'all':
		print("Processing all BA stars in fits file")
		BA = []
		with fits.open(data_file) as ALLDATA:
			for n,i in enumerate(ALLDATA[6].data['TARGID']):
				if 'LR-BA' in i:
					BA.append(i)
	else:
		print("Processing BA stars specified in target list")
		targlist = np.genfromtxt(target_list, dtype=None, encoding=None)
		BA = []
		with fits.open(data_file) as ALLDATA:
			for targ in targlist:
				try: 
					idx=list(ALLDATA[6].data['FIBREID']).index(targ)
					BA.append(ALLDATA[6].data['TARGID'][idx])
				except:
					try:
						idx=list(ALLDATA[6].data['TARGID']).index(targ)
						BA.append(ALLDATA[6].data['TARGID'][idx])
					except:	
						print(str(targ)+": Cant find either FIBREID or TARGID in input table.")

ALLDATA = None
info = None
wavelength = None
targetmask = None

def set_globals():
    global ALLDATA
    global info 
    global wavelength
    global targetmask

    with fits.open(data_file) as ALLDATA:
	    head0 = ALLDATA[0].header
	    info = ALLDATA[6].data
	    data1 = ALLDATA[1].data
	    head1 = ALLDATA[1].header
	    wave0 = head1['CRVAL1']  
	    increm = head1['CD1_1']
	    wavelength = np.array([wave0+increm*a for a in range(len(data1[1]))])
	    targetmask = (wavelength > min_wav) & (wavelength < max_wav)
	    wavelength = wavelength[targetmask]

if multiprocess==True:
	from multiprocessing import Pool
	if __name__ ==  '__main__':
		pool = Pool(processes=process_no, initializer=set_globals)
		it = pool.imap_unordered(mcmc_one, BA)
		for nn,i in enumerate(range(len(BA))):
			it.next()
		make_output_fits()
		print(datetime.now() - startTime)

else:
	set_globals()
	for nn,i in enumerate(BA):
		mcmc_one(i)
	make_output_fits()
	print(datetime.now() - startTime)