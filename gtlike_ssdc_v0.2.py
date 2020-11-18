import numpy as np
from fermipy.gtanalysis import GTAnalysis
from astropy.io import fits
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u
from glob import glob
import os
import sys

####################################################################################################################
# Script for online analysis of Fermi-LAT data using Fermipy - written by R. Angioni (roberto.angioni@ssdc.asi.it) #
####################################################################################################################


########################################
# READ INPUTS AND DEFINE CONFIGURATION #
########################################

argv = sys.argv

NAME = argv[1]
RA = argv[2]
DEC = argv[3]
TSTART = argv[4]
TSTOP = argv[5]
EMIN = argv[6]
EMAX = argv[7]
ID = argv[8]
DATAPATH = argv[9]
RAD = argv[10]

#SC = argv[8]
#ROIu = argv[9]
#xml = argv[10]

#defining analysis configuration
roiSize = RAD
binsz = '0.1'
binsperdec = '10'
evclass = '128'
zmax = '90'
tstart = float(TSTART)
tstop = float(TSTOP)
emin = EMIN
emax = EMAX
srcRA = RA
srcDec = DEC
evtype = '3'
fil = '(DATA_QUAL>0)&&(LAT_CONFIG==1)'
edisp = 'True'
irfs = "'P8R3_SOURCE_V2'"
catalog = '4FGL'


#create working directory
wpath = DATAPATH
os.mkdir(wpath+ID)
os.chdir(wpath+ID)

#define data files
ft1_list = glob(wpath+ID+'*PH*')
ft1 = 'events'+ID+'.txt'
#ascii.write(ft1_list, ft1)
f=open(ft1,'w')
ft1_list=map(lambda x:x+'\n', ft1_list)
f.writelines(ft1_list)
f.close()

ft2 = wpath+ID+'_SC00.fits'
    
#writing configuration file

cFile_name = 'config'+ID+'.yaml'
cFile = open(cFile_name, 'w')
cFile.write("data:\n  evfile : "+ft1+"\n  scfile : "+ft2+"\n")
cFile.write("\nbinning:\n  roiwidth : "+str(roiSize)+"\n  binsz : "+binsz+"\n  binsperdec : "+binsperdec+"\n")
cFile.write("\nselection:\n  tmin : "+str(tstart)+"\n  tmax : "+str(tstop)+"\n  evclass : "+evclass+"\n  ra : "+str(srcRA)+"\n  dec : "+str(srcDec)+"\n  zmax : "+zmax+"\n  emin : "+str(emin)+"\n  emax : "+str(emax)+"\n  evtype : "+evtype+"\n  filter : "+fil+"\n")
cFile.write("\ngtlike :\n  edisp : "+edisp+"\n  irfs : "+irfs+"\n  edisp_disable : ['isodiff']\n")
cFile.write("\nltcube: \n  use_local_ltcube : True\n")
cFile.write("\nmodel:\n  src_roiwidth : "+str(RAD)+"\n  galdiff  : '$FERMI_DIFFUSE_DIR/gll_iem_v07.fits'\n  isodiff  : '$FERMI_DIFFUSE_DIR/iso_P8R3_SOURCE_V2_v1.txt'\n  catalogs : [4FGL]")
cFile.close()

#####################
# START OF ANALYSIS #
#####################

#initializing analysis object
gta = GTAnalysis(cFile_name,logging={'verbosity': 3})

gta.setup(overwrite = True)

#first optimization run with output
fit_res = gta.optimize()

gta.write_roi('fit_optimize_'+ID)

#free parameters for full likelihood fit
gta.free_sources(pars='norm')
gta.free_sources(distance = 3.0)
gta.free_source('galdiff')
gta.free_source('isodiff')


#do the likelihood fit
fit_results = gta.fit()
if fit_results['fit_success']!=True:
    gta.load_roi('fit_optimize_'+ID+'.npy')
    gta.free_sources(free=False)
    gta.free_sources(pars='norm', distance = 3.0)
    gta.free_sources(distance = 1.)
    gta.free_source('galdiff')
    gta.free_source('isodiff')
    fit_res2 = gta.fit()
    if fit_results['fit_success']!=True:
        gta.load_roi('fit_optimize_'+ID+'.npy')
        gta.free_sources(free=False)
        gta.free_sources(pars='norm', distance = 1.5)
        gta.free_sources(distance = 0.5)
        gta.free_source('galdiff')
        gta.free_source('isodiff')
        fit_res2 = gta.fit()
    
fit_res = gta.optimize()

#save results with 4FGL-only model
gta.write_roi('fit_gtlike_'+ID)

cat_results = open('cat_results.txt', 'w')
for s in gta.roi.sources:
    cat_results.write(str(s)+'\n')

cat_results.close()

#look for additional sources based on TS map peaks (TS>25)

#average spectrum
model = {'Index' : 2.2, 'SpatialModel' : 'PointSource'}
find = gta.find_sources(model=model, sqrt_ts_threshold=5, min_separation=0.5, multithread = True)

#hard spectrum
model = {'Index' : 1.8, 'SpatialModel' : 'PointSource'}
find = gta.find_sources(model=model, sqrt_ts_threshold=5, min_separation=0.5, multithread = True)

#soft spectrum
model = {'Index' : 2.4, 'SpatialModel' : 'PointSource'}
find = gta.find_sources(model=model, sqrt_ts_threshold=5, min_separation=0.5, multithread = True)

fit_res = gta.optimize()
fit_res = gta.fit()

#write output
gta.write_roi('fit_srcfind_'+ID)



#localize new sources and add to ds9 region file with error circles

for s in gta.roi.sources:
        if s.name != 'isodiff' and s.name != 'galdiff' and s['ts']>25 and s['offset']<3.:
                soi = s.name
                loc=gta.localize(soi, update=True, free_radius=1)

fit_res = gta.optimize()

gta.write_roi('fit_srcfind_loc_'+ID)



#check if there is a source close to the target position, if not put a test source and fit it

srcname = NAME
if gta.roi.sources[0]._data['offset']>gta.roi.sources[0]._data['pos_r99']:
    print '# No source consistent with target position after localization, creating and fitting test source... #'
    gta.add_source(srcname,{ 'ra' : srcRA, 'dec' : srcDec ,'SpectrumType' : 'PowerLaw', 'Index' : 2.2,'Scale' : 1000, 'Prefactor' : 5.0E-12,'SpatialModel' : 'PointSource' })
    gta.free_source(srcname)
    gta.free_sources(pars='norm')
    gta.free_source('galdiff')
    gta.free_source('isodiff')
    fit_res = gta.optimize()
    fit_res = gta.fit()
    gta.write_roi('fit_testsrc_'+ID)
    print '# Test source fitted successfully... #'
    if gta.roi[srcname]._data['ts']>25:
        print '# ...and we have a significant detection! (TS='+str(gta.roi[srcname]._data['ts'])+') Localizing... #'
        loc=gta.localize(srcname, update=True)
        if loc['fit_success']==True:
            print '# Localization succeeded! #'
            print '#Final optimization run...#'
            fit_res = gta.optimize()
            gta.write_roi('fit_detected_localized_optimized_'+ID)
        else:
            print '# Localization failed! #'
            print '# Running SED...#'
            sed = gta.sed(srcname,make_plots=True)
    else:
        print '# ...but the source is not significantly detected (TS='+str(gta.roi[srcname]._data['ts'])+') #'
        print '# Writing source properties to outfile.dat #'
        outfile=open(ID+'_gtlike.txt','w')
        outfile.write('ROI data for source close to target:\n')
        outfile.write(str(gta.roi.sources[0]))
        outfile.write('\nROI center offset\t:\t'+str(gta.roi.sources[0]._data['offset'])+' degrees')
        if gta.roi[srcname]._data['ts']<25:
            outfile.write('\n95% confidence flux upper limit\t:\t'+str(gta.roi.sources[0]._data['flux_ul95'])+' ph/cm**2/s')
            outfile.write('\n95% confidence energy flux upper limit\t:\t'+str(1.6021766e-6*gta.roi.sources[0]._data['eflux_ul95'])+' erg/cm**2/s')
        else:
            outfile.write('\nEnergy flux in cgs units\t:\t'+str(1.6021766e-6*gta.roi.sources[0]._data['eflux'])+' +\- '+str(1.6021766e-6*gta.roi.sources[0]._data['eflux_err'])+' erg/cm**2/s')
        outfile.close()
        gta.write_roi('fit_final_'+ID, make_plots = True)
else:
    print '# Source',gta.roi.sources[0].name,'is positionally consistent with target position! Optimizing and writing output...#'
    gta.free_source(gta.roi.sources[0].name)
    fit_res = gta.optimize()
    fit_res = gta.fit()
    gta.write_roi('fit_final_'+ID)
    outfile=open(ID+'_gtlike.txt','w')
    outfile.write('ROI data for source close to target:\n')
    outfile.write(str(gta.roi.sources[0]))
    if gta.roi.sources[0]._data['ts']<25:
        outfile.write('\n95% confidence flux upper limit\t:\t'+str(gta.roi.sources[0]._data['flux_ul95'])+' ph/cm**2/s')
        outfile.write('\n95% confidence energy flux upper limit\t:\t'+str(1.6021766e-6*gta.roi.sources[0]._data['eflux_ul95'])+' erg/cm**2/s')
    else:
        outfile.write('\n68% error circle\t:\t'+str(gta.roi.sources[0]._data['pos_r68'])+' degrees')
        outfile.write('\n95% error circle\t:\t'+str(gta.roi.sources[0]._data['pos_r95'])+' degrees')
        outfile.write('\n99% error circle\t:\t'+str(gta.roi.sources[0]._data['pos_r99'])+' degrees')
    outfile.write('\nROI center offset\t:\t'+str(gta.roi.sources[0]._data['offset'])+' degrees')
    outfile.close()
    print gta.roi.sources[0]
    e_min = gta.config['selection']['emin']
    e_max = gta.config['selection']['emax']
    srcname = gta.roi.sources[0].name
    if gta.roi.sources[0]._data['ts'] > 10.:
        if gta.roi.sources[0]._data['ts'] > 100.:
            n_sed_bins = 10
        else:
            n_sed_bins = gta.roi.sources[0]._data['ts']/10.
        binsz_sed = (np.log10(e_max)-np.log10(e_min))/n_sed_bins
        sed_bin_edges = np.asarray([0.]*(int(n_sed_bins)+1))
        sed_bin_edges[0]=np.log10(e_min)
        for i in xrange(1,int(n_sed_bins)+1):
            sed_bin_edges[i] = sed_bin_edges[i-1]+binsz_sed
        sed_bin_edges[-1]=np.log10(e_max)
        sed = gta.sed(srcname,loge_bins=sed_bin_edges, make_plots=True)


print '###########################'
print '# Computing HE photons... #'
print '###########################'

print '# Running gtdiffrsp... #'
os.system('gtdiffrsp evfile = ft1_00.fits scfile = %s srcmdl = fit_final_%s_00.xml irfs = P8R3_SOURCE_V2'%(ft2,ID))

print '# Running gtsrcprob... #'
os.system('gtsrcprob evfile = ft1_00.fits scfile = %s outfile = srcprobs_%s.fits srcmdl= fit_final_%s_00.xml irfs = P8R3_SOURCE_V2'%(ft2,ID,ID))

probs_file = fits.open('srcprobs_'+ID+'.fits')
probs = probs_file[1].data
probs_e = probs['ENERGY']/1.e3
probs_t = probs['TIME']
probs_p = probs[gta.roi.sources[0].name]

probs_e2 = []
probs_t2 = []
probs_p2 = []

for i in xrange(0,len(probs_e)):
    if probs_p[i] > 0.80 and probs_e[i]> 10.:
        probs_e2.append(probs_e[i])
        probs_t2.append(Time(51910.0+probs_t[i]/86400., format = 'mjd'))
        probs_p2.append(probs_p[i])


if probs_e2 != []:
    outfile=open(ID+'_gtlike.txt','a')
    outfile.write('\nE > 10 GeV photons with p > 0.8:')
    for i in xrange(0,len(probs_e2)):
        outfile.write('\nE = %.2f GeV; time : %s; prob = %f'%(probs_e2[i], probs_t2[i].iso, probs_p2[i]))
    outfile.close()


    
#cleaning up
for cleanup in glob('*.*'):
    if 'fit_final' not in cleanup and 'gtlike' not in cleanup:
        os.remove(cleanup)




