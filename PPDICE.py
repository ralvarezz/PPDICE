from scipy.ndimage import binary_propagation as reco
from scipy.ndimage import binary_erosion as ero
from scipy.signal import find_peaks as fp
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
from astropy.table import Table
from datetime import datetime
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import multiprocessing
import numpy as np
import matplotlib
import warnings
import random
import glob
import time
import os

warnings.filterwarnings('ignore', category=UserWarning, append=True)
np.seterr(all=None, divide=None, over=None, under=None, invalid=None)


#fix the seed to generate the same spectrum plots in case for comparison between different tries
np.random.seed(10)

#################################################################
#                                                               #
#                   User defined parameters                     #
#                                                               #
#################################################################



#define cube filename 
data_location='demo/'
c = data_location+'G353_feathered_n2hp_pbcorr.fits'

name = 'G353'#region name to be used in the written files

#variable to distinctively separate different runs
date = datetime.today().strftime('%d%m%Y_%H:%M:%S')

#use the total number of threads for the parallelization
#or set it to the number of threads you want to use
nprocess = os.cpu_count()

#diagnostic plots, created at random positions for diagnostics
nplots = 300

#number of points on each side of a given element for the rolling average
#implemented in the function sm
smooth_half_width=2
        
#make matplotlib non interactive to avoid using more resources than needed
matplotlib.use('Agg')

SNR = True #create snr map of the region

#create mean spectra maps
mean_maps = True

#create a raw moment 0 map in order to get a wcs frame
frames=True

#create plots comparing the spectra of the full cube and the cube with the isolated component only
comparison = True

#channel limits to where to estimate the rms
#0->rms_chan_i; rms_chan_u->last chanell 
rms_chan_i = 44
rms_chan_u = 170

#For spectra with SNR<snr_threshold the isolated component will be extracted based on a mean measurement. A refined approach is implemented for SNR>snr_threshold
snr_threshold = 4

#values for the colorscales
std_scale = 5
lw_contour = 0.7

#rough vlsr of the cloud
vlsr = -18


#################################################################
#                                                               #
#    The rest of the code runs based on the parameters above    #
#                                                               #
#################################################################



outfolder = data_location+name+'/'+date#'_width_'+str(width_vel)+

os.makedirs(outfolder,exist_ok=True)


def sm(array,npoints):
    """
    Simple function to smooth a spectrum previous to the usage of find_peaks. This approach reduces the chance of identifying a peak or valley due to noise.
    array: array to smooth
    npoints: half width of the rolling average
    Returns the smoothed array
    """    
    lz = len(array)
    smoothed_array = np.copy(array)
    for z in range(lz):
        if (z<npoints) or (z>lz-npoints) :
            array[z]=np.nan
        else:
            mean_window_value = np.nanmean(array[z-npoints:z+npoints+1])
            smoothed_array[z] = mean_window_value
    return smoothed_array

def sm3d(cube,npoints):
    """
    Simple function to smooth a datacube 'array' by averaging the values ± npoints. The same as the 'sm' function but for a cube.
    Returns the smoothed cube
    """
    lx,ly,lz = cube.shape[::-1]
    data_footprint = np.isfinite(np.nanmax(cube,axis=0))
    smoothed_cube = np.copy(cube)
    zero_plane = np.zeros(lx*ly).reshape(ly,lx)
    zero_plane[~data_footprint]=np.nan
    
    for z in range(lz):
        if (z<npoints) or (z>lz-npoints) :
            smoothed_cube[z]=zero_plane
        else:
            mean_window_value = np.nanmean(cube[z-npoints:z+npoints+1],axis=0)
            smoothed_cube[z] = mean_window_value
    return smoothed_cube


def create_hdu(array):
    """
    Function to create a hdu with celestial coordinates. This is used for the auxiliary fits file.
    """
    
    hdu = fits.PrimaryHDU()
    hdu.data = 1*array#works for the mask to convert from bool to int, doesn't affect floats

    hdu.header['BITPIX'] = -64
    
    shape = len(array.shape)
    if shape==3:
        lz,ly,lx = array.shape
    else:
        ly,lx = array.shape

    
    hdu.header['NAXIS'] = shape
    hdu.header['NAXIS1'] = lx#mask_snr.shape[1]
    hdu.header['NAXIS2'] = ly#mask_snr.shape[0]
    hdu.header['NAXIS3'] = lz#mask_snr.shape[0]

    hdu.header['DATATYPE'] = 'Auxiliary Information'
    hdu.header['WCSAXES'] = (               shape , 'Number of coordinate axes'                       )
    hdu.header['CRPIX1']  = (cube_header['CRPIX1'], 'Pixel coordinate of reference point'             )
    hdu.header['CRPIX2']  = (cube_header['CRPIX2'], 'Pixel coordinate of reference point'             )
    hdu.header['CDELT1']  = (cube_header['CDELT1'], '[deg] Coordinate increment at reference point'   )
    hdu.header['CDELT2']  = (cube_header['CDELT2'], '[deg] Coordinate increment at reference point'   )
    hdu.header['CUNIT1']  =  ('deg'                , 'Units of coordinate increment and value'         )
    hdu.header['CUNIT2']  = ('deg'                , 'Units of coordinate increment and value'         )
    hdu.header['CTYPE1']  = ('RA---SIN'           , 'Right ascension, orthographic/synthesis project' )
    hdu.header['CTYPE2']  = ('DEC--SIN'           , 'Declination, orthographic/synthesis projection'  )
    hdu.header['CRVAL1']  =  (cube_header['CRVAL1'], '[deg] Coordinate value at reference point'       )
    hdu.header['CRVAL2']  =  (cube_header['CDELT2'], '[deg] Coordinate value at reference point'       )
    hdu.header['LONPOLE'] =  (cube_header['LONPOLE'] , '[deg] Native longitude of celestial pole'        )
    hdu.header['LATPOLE'] =  (cube_header['LATPOLE'], '[deg] Native latitude of celestial pole'         )
    hdu.header['RESTFRQ'] =  (cube_header['RESTFRQ'], '[Hz] Line rest frequency'                        )
    
    for l in range(len(planes)):
        hdu.header['PLANE '+str(l)] = planes[l]

        
    fits.writeto(outfolder+'/'+name+'_auxiliary_cube.fits',hdu.data,hdu.header,overwrite=True)


def extract_isolated_components(x,y,saveplots):
    
    """
    Function that does the isolated extraction procedure for a given spectrum with position x,y in a cube.
    """
    
    
    #define placeholder arrays for the intensity, velocity offsets, and weights
    I = np.array((0,0,0))
    V = np.array((1e1,1e1,1e1))
    dV = np.array((1e1,1e1,1e1))
    weights = np.array((0,0,0)) 
    
    saveplots=bool(saveplots)
    
    spec = data_raw[:,y,x]
    
    
    rms_spec = np.copy(spec)
    
    #remove range of values where there is cloud emission to only get
    #the channels for the rms estimation
    rms_spec[rms_chan_i:rms_chan_u]=np.nan

    rms=np.nanstd(rms_spec,axis=0)  

    #length/2 of the box used to smooth
    smooth_half_width=2

    #grab only the velocities where there is cloud emission
    sm_spec = sm(spec,smooth_half_width)
    sm_spec[:rms_chan_i]=0
    sm_spec[rms_chan_u:]=0            

    if saveplots:
        #only create plots if saveplots==True,
        #so I don't waste time creating plots for every pixel and only 
        #saving some of them
        
        plt.clf()
        plt.plot(vel,spec,label='Data')#,c='k',linestyle='dashed')
        plt.plot(vel,sm_spec,'r--',label='Smoothed data')                                                    
        plt.axhline(rms,c='blue',linestyle='dashed',label='rms')
    else:
        pass
    
    #for the data, only implement the dip finder method if the pixel
    #is inside the snr mask, so data[:,y,x].all()==np.nan
    #if np.isnan(spec).all()==True:#data SNR<snr_threshold
    width = mean_width#2.5 #/res#km/s / km/s /pixel
    if snr_map[y,x]<snr_threshold:   
        #for spectra with SNR < threshold take a velocity cut based on the mean dip.
        finder=False
        
        
        #For this case, some variables for the title don't exist so when False, I use a different title
        # create placeholder for the peaks  in order to fill up that array
        #and only to avoid the code failing later insert the raw spectra until the first dip into the new array for the isolated component only, so it takes data from channel 1 up to vel_dip_idx (location of the mean dip)
        peaks=np.array([[1,1],[1,1]])
        
        
        cube_isolated[:vel_dip_idx+1,y,x] = data_raw[:vel_dip_idx+1,y,x]
            
        
        #include noisy channels at the end of the cube, after the line emission in case it's useful for better rms estimation
        
        cube_isolated[rms_chan_u:,y,x] = data_raw[rms_chan_u:,y,x]
        
        if saveplots:
            plt.axvline(vel_dip,c='k',linestyle='dashed',label='mean dip at ' + str(np.round(vel_dip,2))+'km/s')                    
            plt.axvline(vel_dip,c='orange',linestyle='dotted',alpha=0.7,label='final upper limit' + str(np.round(vel_dip,2))+'km/s')                
 
        else:
            pass
    else:
        #dip finder procedure
        finder=True
        
        #find peaks then define a height relative close to the isolated component

        #pre processing of the data
        #we only care for the data around the intensity dip between the isolated and main component. So to avoid missidentifications we remove values after the vlsr
        data_removed_after_vlsr = np.copy(sm_spec)                    
        data_removed_after_vlsr[mean_peak_iso+int(width/res):]=0
        
        #the height conditions says that the peak of the isolated 
        #component has to be between 2rms and the 90% the intensity of the main line for continuum holes
        n_height_peak_finder = 2 
        #could be number larger to avoid the noise, after the smoothing a value of 2 worked well with the tested data
        
        peaks = fp(data_removed_after_vlsr, height=n_height_peak_finder*rms)
        

        #find minima on the right of the isolated component
        # spec = spec*-1

        if len(peaks[0])==0:
            #this can happen when the line has high snr but the emission seems flat line with a top-hat-like shape 
            plt.clf()
            pass
        else:
            idx_peak = peaks[0][0]

            #distance of 2km/s between peaks
            minima = fp(sm_spec*-1,width=2,distance=int(2/res))#,distance=6/res)#,width=0.5/res)#,width=width)

            #array of velocities of possible dip for this specific spectrum
            vel_min = vel[minima[0]]
            
            if len(vel_min)==0:
                pass
            else:
                
                #first guess
                window_lim = vel[idx_peak]+width
                
                #index of vel_min where the dip is located
                #take dip closer to the mean dip
                
                #2nd guess
                idx = abs(vel_min-vel_dip).argmin()
                found_dip_idx = abs(vel_min[idx]-vel).argmin()

                #3rd guess,
                
                idx_3rd_guess = abs(vel-(vel_dip+width)).argmin()
                
                dist_window_to_mean_dip = abs(window_lim-vel_dip)
                #distance from dip for this specific spectra to the mean dip
                dist_found_dip_to_mean_dip = abs(vel[found_dip_idx]-vel_dip)
                #distance from dip to 3rd guess
                distance_from_dip_to_third_guess = abs(vel[idx_3rd_guess]-vel_dip)

                #velocities of possible limits of the isolated component
                limit_value = vel[abs(vel-(window_lim)).argmin()],vel[found_dip_idx],vel[idx_3rd_guess]
                
                #idx of limit_value elements relative to velocity array

                limit_value_idx = abs(window_lim-vel).argmin(), abs(vel[found_dip_idx]-vel).argmin(),idx_3rd_guess
                
                #define an array with the distances between the mean dip to the window limit and to the found dip 
                print(x,y,saveplots)
                #print('#####')
                #print(window_lim,vel_min[idx],vel_dip+mean_width)
                V = np.array((window_lim,vel_min[idx], vel_dip+mean_width))
                
                dV = np.array((dist_window_to_mean_dip,dist_found_dip_to_mean_dip,distance_from_dip_to_third_guess))
                
                
                
                
                I_window = abs(sm_spec[limit_value_idx[0]])
                I_dip = abs(sm_spec[limit_value_idx[1]])
                I_third = abs(sm_spec[limit_value_idx[2]])

                
                I = np.array((I_window,I_dip,I_third))
                
                #by normalizing the parameters we get rid of the units and compare each element of the intensity and distances
                I_normed = I/min(I)
                
                #I have to normalize by the minimum value, given it'll leave the minimum guesses at 1
                
                #when normed, if the closest guess is has a distance of ~0, then the other guess has always a multiple of the channel resolution (steps of 0.23km/s for the N2H+ ALMA-IMF data I tested)
                
                dV_normed = dV/min(dV)
                
                #these constants indicate how relevant is the variation in velocity and intensity. If the guesses are a lot different in velocity, but similar in intensity, then the variation in velocity should play a more important role in the 'weights'
                #These were defined by inspection of the diagnostic plots when saveplots = True
                weights = 1/(dV_normed*0.8+I_normed*0.2)
                
                idx_guess_with_higher_weight = limit_value_idx[weights.argmax()]
                
                choosen_guess = vel[idx_guess_with_higher_weight]

                
                #Select the final guess as the one with higher weight
                I_final = sm(spec,smooth_half_width)[idx_guess_with_higher_weight]
                vel_final = choosen_guess
                
                ##make plots of the weights insert into 2D arrays the intensity of the guesses and the final guess 
                #values in order to plot them
                
                #closest_limit_to_mean_dip_idx = limit_value_idx[distances.argmin()]

                #higher_weight = weights.argmax()
                
                #limit_velocity_idx = abs(vel-value_upper_limit).argmin()
                
                #closest_limit_to_mean_dip_idx=limit_velocity_idx


                '''modification on august 17th                    
                
                #idx of the limit with less distance to the mean dip
                closest_limit_to_mean_dip_idx = limit_value_idx[distances.argmin()]
                value_upper_limit = vel[closest_limit_to_mean_dip_idx]
                
                #final choosen limit:
                final_limit = vel[closest_limit_to_mean_dip_idx]

                #sometimes the window limit of the isolated component
                #is closer to the mean dip than the dip found for that specific
                #spectrum.
                #If that is the case, then choose the velocity with the lowest 
                #intensity
                
                                    
                #find velocity value of the smaller distance to the mean dip
                #so obtain the array "plus" or "vel_min[idx]"
                
                value_upper_limit = vel[closest_limit_to_mean_dip_idx]
                
                #find element index for this velocity
                #this is the orange dashed line
                limit_velocity_idx = abs(vel-value_upper_limit).argmin()
                

                #check if selected guess has the lowest associated
                #intensity
                
                #limit_velocity_idx = 
                
                if I[distances.argmin()]>I[distances.argmax()]:
                    value_upper_limit = limit_value[distances.argmax()]
                    limit_velocity_idx = abs(vel-value_upper_limit).argmin()
                else:#do nothing
                    pass
                #if both distances from the guesses to the mean vel dip
                #are the same, then use the one with lower intensity
                
                if distances[0]==distances[1]:
                    value_upper_limit = limit_value[I.argmin()]
                    limit_velocity_idx = abs(vel-value_upper_limit).argmin()
                else:
                    pass
                '''
                
                
                
                
                cube_isolated[:idx_guess_with_higher_weight+1,y,x] = data_raw[:idx_guess_with_higher_weight+1,y,x]
                    
                cube_isolated[rms_chan_u:,y,x] = data_raw[rms_chan_u:,y,x]
                
                #print('inserting isolated component with snr>snr_threshold')      
                #print('isolated_component for (x,y)=',cube_isolated[:,y,x])
                if saveplots:
                    plt.axvline(vel[idx_peak],c='g',alpha=0.8,label='Isolated component') #isolated peak
                    print(weights)
                    plt.axvline(vel_dip,c='k',linestyle='dashed',label='mean dip at ' + str(np.round(vel_dip,2))+'km/s')
                     
                    plt.axvline(vel[idx_3rd_guess],c='darkred',linestyle='solid',label='3rd guess' + str(np.round(vel[idx_3rd_guess],2))+'km/s')
                                              
                     
                    #find minima
                    #find peaks of inverted spectra (minima)
                    #then find the peak closes to 22km/s (dip in the spectra of the isolated
                    #component
                    #pass
                    plt.axvline(vel_min[idx],c='k',linestyle='solid',alpha=0.7)

                    #choose between black and green vertical line
                    plt.axvline(vel[found_dip_idx],c='gray',linestyle='solid',alpha=1,label='2nd guess, '+str(np.round(vel[found_dip_idx],2))+'km/s')

                    #plt.axvline(minus,c='g',alpha=0.7,linestyle='dashed')
                    plt.axvline(vel[limit_value_idx[0]],c='g',alpha=0.7,linestyle='dashed', label='limit fixed velocity window' +str(np.round(vel[limit_value_idx[0]],2))+'km/s')
                    #plt.axhline(rms,c='blue',linestyle='dashed',label='rms')
                    
                    plt.axvline(vel_final,c='orange',linestyle='dotted',alpha=0.7,label='final upper limit' + str(np.round(vel_final,2))+'km/s')
                else:
                    pass
    if saveplots:  
        #print('#####',peaks,x,y,saveplots)

        try:
            len(peaks[0])
        except:
            print(peaks,x,y,saveplots)

        if len(peaks[0])==0:
            plt.clf()
        else:
            if finder:                    
                plt.title('(X,Y) = ('+str(x+1)+','+str(1+y)+ '), $|$I(window)$|$=' +str((np.round(I_window,2)))+'K , $|$I(dip)$|$=' +str((np.round(I_dip,2)))+'K, $|$I(3rd)$|$=' +str((np.round(I[2],2)))+'K',size='medium')
            else:
                plt.title('(X,Y) = ('+str(x+1)+','+str(1+y)+')')
                
            plt.xlabel('Vel [km/s]')
            plt.ylabel('Intensity [K]')
            legend = plt.legend(loc=2,fontsize='small')
            legend.get_frame().set_alpha(0.3)
    
            plt.savefig(outfolder+'/spectra_plots/x'+str(x+1)+'_y'+str(y+1)+'.png',dpi=350)
            print(name,'('+str(x+1)+','+str(y+1)+')',(yy*x +y)/(xx*yy),' snr=',np.round(snr[y,x],3))#,'I=',I,', I_normed=',I_normed,', dv=',distances,'  dv_normed=',dV_normed,'weights=',weight)
    
            #create comparison plot
            plt.clf()
            
            spectra_raw = plt.plot(vel, data_raw[:,y,x], c='k',label='Original spectrum')
            
            plt.plot(vel,cube_isolated[:,y,x],'r--', label='Isolated component')
            
            plt.xlabel('Velocity [km/s]')
            plt.ylabel('Intensity [K]')
            if finder:
                i1 = str(np.round(I[0]      ,2))
                i2 = str(np.round(I[1]      ,2))
                i3 = str(np.round(I[2]      ,2))
                dV1 = str(np.round(dV[0]      ,2))
                dV2 = str(np.round(dV[1]      ,2))
                dV3 = str(np.round(dV[2]      ,2))
                
                w1 = str(np.round(weights[0],2))
                w2 = str(np.round(weights[1],2))
                w3 = str(np.round(weights[2],2))
                
                #if I use the finder then determine the snr of the found peak
                snr_peak =np.round(spec[idx_peak]/rms,3)
                plt.title('(X,Y)=('+str(x+1)+','+str(y+1)+r'),SNR$_{peak}$='+str(snr_peak)+r', SNR$_{isolated part}$='+str(np.round(snr_map[y,x],2))+ r', $I_{g1}=$'+i1+'K'+r', $I_{g2}=$'+i2+'K'+r', $I_{g3}=$'+i3+'K'+ '\n'+r'$\Delta V_{g1}$='+dV1+'km/s, '+r'$\Delta V_{g2}$='+dV2+r'km/s, $\Delta V_{g3}$='+dV3+r'km/s, $W_{g1}=$'+ w1+ r',$W_{g2}=$'+w2+r', $W_{g3}=$'+w3,size='medium')
            else:
                #otherwise use the snr of the data with v < vel_dip
                snr_peak = snr[y,x]
                plt.title('(X,Y) = ('+str(x+1)+','+str(y+1)+r'), SNR for $V<V_{lsr}$='+str(snr_peak))
            
            plt.legend(loc=1)
            
            
            plt.savefig(outfolder+'/comparison/x'+str(x+1)+'_y'+str(y+1)+'.png',dpi=350)
            plt.close()
            
    else:
        pass
    return [y,x,cube_isolated[:,y,x],(I[weights.argmax()],V[weights.argmax()],dV[weights.argmax()],weights[weights.argmax()]),I,V,dV,weights]




noisee = np.array(())
namee = np.array(())
vell = np.array(())

#names of the planes to save in the auxiliary fits file
planes = 'SNR mask','Final intensity','Final velocity','Final dV','Final weight', 'Guess 1 Intensity','Guess 2 Intensity','Guess 3 Intensity', 'Guess 1 V ','Guess 2 V','Guess 3 V','Guess 2 dV','Guess 3 dV','Weight guess 1','Weight guess 2','Weight guess 3'


plots_done = []
start_time = time.time()

print(name,'snr cut = ',snr_threshold )

hdu = SpectralCube.read(c)
res = hdu.hdulist[0].header['CDELT3']/1000 #spectral resolution in km/s

cube_header = hdu.hdulist[0].header




data = hdu.hdulist[0].data    
lz,ly,lx = data.shape[0],data.shape[1],data.shape[2]

# define placeholder arrays to fill up as the code runs

#create auxiliary cube full of nans to later fill in with the parameters
auxiliary_cube = np.zeros(lx*ly*len(planes)).reshape(len(planes),ly,lx)
auxiliary_cube[True]=np.nan

cube_isolated = np.zeros(lx*ly*lz).reshape(lz,ly,lx)

guess_intensity = np.zeros(lx*ly*3).reshape(3,ly,lx)
guess_cut_velocity = np.copy(guess_intensity)
guess_distance = np.copy(guess_intensity)
guess_weight = np.copy(guess_intensity)
final_guesses = np.copy(guess_intensity) 


nsize = lx


data_raw = np.copy(data)


lx,ly = data.shape[2],data.data.shape[1]
#generate positions
lx = np.arange(lx)
ly = np.arange(ly)

vel = hdu.spectral_axis.value/1000


vlsr_idx = abs(vlsr-vel).argmin()

yy,xx = data.shape[1:]

noise_chans = np.copy(data)
noise_chans[rms_chan_i:rms_chan_u]=np.nan
noise_map = np.nanstd(noise_chans,axis=0)

print(name,', noise = ',np.nanmean(noise_map,axis=(0,1)))


spectra = np.nanmean(data,axis=(1,2))

#data array has nans outside mask_snr
chan_number = np.arange(len(spectra))+1

#to identify a minima invert the mean spectrum
flipped_spectra = -1*sm(spectra,smooth_half_width)

#make noise channels equal to 0 to avoid false positive peaks
flipped_spectra[:rms_chan_i]=0
flipped_spectra[rms_chan_u:]=0

peaks = fp(flipped_spectra,width=5)#distance=7/0.23)#,height=rms_map*3)

vel_dip = vel[peaks[0][0]+1] #velocity of the mean dip
vel_dip_idx = abs(vel-vel_dip).argmin()


rms_mean_spectrum = np.copy(spectra)
rms_mean_spectrum[rms_chan_i:rms_chan_u]=np.nan
rms_mean_spectrum = np.nanstd(rms_mean_spectrum)
mean_peak = fp(sm(spectra,smooth_half_width),height=snr_threshold*rms_mean_spectrum)#rms_mean_value_map)   

mean_peak_iso = mean_peak[0][0]

mean_width = np.round(abs(vel_dip-vel[mean_peak_iso]),3)

isolated_chan = peaks[0][0]+1#+1 because is used in slicing arrays

if mean_maps:
    plt.figure()
    plt.plot(vel,spectra,c='blue',label='Mean spectrum')
    plt.plot(vel,sm(spectra,smooth_half_width),c='r',label='smoothed spectrum',alpha=0.8)

    plt.axvline(vel[mean_peak_iso],c='green',label='mean isolated peak at '+str(np.round(vel[mean_peak_iso],2))+'km/s')

    plt.axvline(vel_dip,label='mean dip at chan = '+str(peaks[0][0])+  ',' + str(np.round(vel_dip,2))+' km/s',c='k',linestyle='dashed')
    
    plt.axvline(vel[rms_chan_i],c='gray',label='Emission limits')
    plt.axvline(vel[rms_chan_u],c='gray')

    plt.ylabel('Intensity [K]')
    plt.xlabel('Velocity [km/s]')
    
    legend = plt.legend(loc=2,fontsize='small')

    plt.title(name+', width = '+str(mean_width)+'km/s')

    os.makedirs(outfolder+'/mean_spectra',exist_ok=True)
    plt.savefig(outfolder+'/mean_spectra/'+name+'.png',dpi=350)


#estimate the SNR map for the mean isolated component
snr_map = np.nanmax(data[:isolated_chan],axis=0)/np.nanstd(noise_chans,axis=0)        


snr = np.copy(snr_map)#in order to use snr and leave snr_map intact
peak_cube =  np.nanmax(data,axis=0)

mask_snr = snr>snr_threshold

#erode the mask, to remove single pixel masks, generally on the edges of the images
eroded = ero(mask_snr,iterations=1)
#use the original mask to obtain a less agressive erotion
mask_snr = reco(eroded,mask=mask_snr)


#create plots comparing the full spectrum and the extracted component at random to avoid plotting the spectrum in every pixel
#every spectrum should have a 1% change of being plotted
prob_to_plot = np.zeros(100)
prob_to_plot[0]=1


#fill data for the table
namee = np.append(namee,name)
noisee = np.append(noisee,np.nanmean(noise_map))
vell = np.append(vell, str(np.round(vel[0])) +' - ' + str(np.round(vel[rms_chan_i]))+', '+str(np.round(vel[rms_chan_u])) +' - ' + str(np.round(vel[len(vel)-1])))

for z in range(data.shape[0]):
    #remove spectra with snr < snr_threshold
    data[z][~mask_snr]=np.nan

if SNR:
    #create and save SNR map
    plt.close('all')
    plt.figure()
    
    mean = np.nanmean(noise_map)
    std = np.nanstd(noise_map)

    plot = plt.imshow(noise_map,origin='lower',cmap='jet',vmin=mean-2*std, vmax=mean+3*std)
    plt.colorbar(plot,label='RMS [K]',fraction=0.046)

    plt.xlabel('X pixel')
    plt.ylabel('Y pixel')


    plt.title(name)

    os.makedirs(outfolder+'/noise_maps',exist_ok=True)

    plt.savefig(outfolder+'/noise_maps/'+name+'_snr_map_isolated.png',dpi=350)
   
    plt.clf()    

    a = snr[np.isfinite(snr)]
    mean = np.nanmean(a)
    std = np.nanstd(a)

    plot = plt.imshow(snr,origin='lower',cmap='jet',vmin=0,vmax=mean+std_scale*std)
    plt.colorbar(plot,label='SNR',fraction=0.046)

    plt.contour(mask_snr,levels=[1,1000], colors=['black'], linewidths=lw_contour)

    plt.xlabel('Y pixel')
    plt.ylabel('X pixel')

    plt.title(name+' - SNR contour at '+str(snr_threshold))
    os.makedirs(outfolder+'/snr_maps',exist_ok=True)

    plt.savefig(outfolder+'/snr_maps/'+name+'_snr_map_isolated.png',dpi=350)



mask = np.isfinite(mask_snr) #detect where the data are located

#backup this code in the outfolder to track changes if needed
os.system('cp multiproc_peaks.py '+outfolder)
os.makedirs(outfolder+'/comparison',exist_ok=True)
os.makedirs(outfolder+'/spectra_plots',exist_ok=True)

#create array of all the pixels combinations
lx_n=np.arange(len(lx))
ly_n=np.arange(len(ly))


#create a list for the 3 arguments
#x,y,saveplots

parm_list = []
#this list should cover all the pixels in the region and randomly parsing saveplots values
pool = multiprocessing.get_context('fork').Pool(nprocess) 
#by default nprocess is the number of cpu threads

# fill out the parm_list that then each element will be parsed to the 'extract_isolated_components' isolated extraction function
for x in range(len(lx_n)):
    for y in range(len(ly_n)):
        spec = data_raw[:,y,x]
        if np.isnan(spec).all():
            #ignore pixels with no data
            pass            
        else:
            if len(plots_done)<=nplots:           
                saveplots = np.random.choice(prob_to_plot)
            else:
                saveplots=False
            parm_list.append((x,y,saveplots))  
            if saveplots==1:
                plots_done.append(1)
        
#async to not wait for processes to finish before going to another pixel, as each extraction is independent from its neighbours
result = pool.starmap_async(extract_isolated_components,iterable=parm_list)


pool.close()#wait for all the processes to finish

print('Closing the pool')
print('Waiting for all the processes to finish')
pool.join()

output = result.get()


#assign obtained parameters to 2d arrays
for parameters in output:
    y = parameters[0] 
    x = parameters[1] 
    cube_isolated[:,y,x] = parameters[2]
    
    for counter in range(3):
        final_guesses[counter][y,x] = parameters[3][counter]
        guess_intensity[counter][y,x] = parameters[4][counter]
        guess_cut_velocity[counter][y,x] = parameters[5][counter]
        guess_distance[counter][y,x] = parameters[6][counter]
        guess_weight[counter][y,x] = parameters[7][counter]


#create cube that includes only data>SNR_threshold
cube_isolated_snr = np.copy(cube_isolated)
#replacing zeroes with with nans outside the region for the output cubes
for ll in range(cube_isolated.shape[0]):
    cube_isolated[ll][~mask]=np.nan#put nans outside the data
    cube_isolated_snr[ll][~mask_snr]=np.nan

for ll in range(3):#the dip finder procedure applies only to data>snr so I clean
    #the cubes below
    guess_intensity[ll][~mask_snr]=np.nan
    guess_distance[ll][~mask_snr]=np.nan
    guess_weight[ll][~mask_snr]=np.nan
    final_guesses[ll][~mask_snr]=np.nan


auxiliary_cube[0]=mask_snr
auxiliary_cube[1:4]=final_guesses    
auxiliary_cube[4:7]=guess_intensity
auxiliary_cube[7:10]=guess_cut_velocity    
auxiliary_cube[10:13]=guess_distance    
auxiliary_cube[13:]=guess_weight
create_hdu(auxiliary_cube)#,'SNR Mask')

#readme file:
#describing the cube with all the parameters, file name = G353+_auxiliary_info
#and suggest to use the mask

fits.writeto(outfolder+'/'+name+'_isolated_component.fits', data = cube_isolated , header=hdu.hdulist[0].header, overwrite=True)

os.makedirs(outfolder+'/moments_isolated',exist_ok=True)


hdu = SpectralCube.read(outfolder+'/'+name+'_isolated_component.fits')
#obtain their moment maps
for l in range(3):
    moment = hdu.spectral_slab(vel[rms_chan_i]*u.km/u.s, vlsr*u.km/u.s).moment(order=l)
    if l==2:
        moment = moment**0.5
    moment = moment/1000
    moment.write(outfolder+'/moments_isolated/'+name+'_moment'+str(l)+'.fits',overwrite=True)
    print('moment ',l,' done')    

#save into a fits file the isolated cube but only the data>SNR_threshold
#header is ok
fits.writeto(outfolder+'/'+name+'_isolated_component_cleaned_snr.fits', data = cube_isolated_snr , header=hdu.hdulist[0].header, overwrite=True)  

os.makedirs(outfolder+'/moments_isolated_snr',exist_ok=True)

hdu = SpectralCube.read(outfolder+'/'+name+'_isolated_component_cleaned_snr.fits')
#create moment maps
#headers are ok too
for l in range(3):
    moment = hdu.spectral_slab(vel[rms_chan_i]*u.km/u.s, vlsr*u.km/u.s).moment(order=l)
    if l==2:
        moment = moment**0.5
    moment = moment/1000
    moment.write(outfolder+'/moments_isolated_snr/'+name+'_moment'+str(l)+'.fits',overwrite=True)
    print('moment ',l,' done')    


folder = data_location+name+'/'+date
table = Table([namee,noisee,vell],names=['Region','RMS [K]','RMS velocity range [km/s]'])
table.write(folder+'/'+name+'_table_rms.dat',format='ascii.fixed_width',overwrite=True)

print(name,' --- ',np.round((time.time() - start_time)/60,3),' minutes --- ')
