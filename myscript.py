# Imports
from astropy.io import fits
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import scipy.optimize
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.colors as colors

def fileintaker(filterstring1, filterstring2):
  '''
  Function extracts all FITS files from folder structure using glob, organises them by given filter name

  Input parameters:
  filestring1 : dtype = string
    Number used to ID filter
    Used to match subdirectory names within folder
  filestring2 : dtype = string
    Number used to ID filter
    Used to match subdirectory names within folder

  Returns:
  filelistlist : dtype = list of lists of strings
    A list containing two sublists
    - filelist1 = List of file paths for FITS files matching filestring1
    - filelist2 = List of file paths for FITS files matching filestring2

  Notes:
  - Function expects folder structure: data/{filter}/*.fits
  '''

  filelist1 = sorted(glob(f"data/{filterstring1}/*.fits"))
  filelist2 = sorted(glob(f"data/{filterstring2}/*.fits"))

  filelistlist = [filelist1, filelist2]

  return filelistlist

#####

def rayremovingfluxconverter(filelist):
  '''
  Function that takes list of FITS files, converts FITS file data to flux data, clips negative flux values, and combines images using median stacking to remove cosmic rays.
  Function parses data headers of FITS file to extract relevant parameters, namely exposure time and photflam, for use in treatment of flux data

  Input parameters:
  filelist : dtype = list of strings
    List of file paths to relevant FITS files. Each file should contain the header and data in ext 1

  Returns:
  median_flux : dtype = array
    Median stack of flux image data with same dimensions as input FITs images, with negative values clipped to 0
  median photflam : dtype = float
    median photflam value across all input FITS files, for later use in magnitude calculation

  Notes:
  - If the FITS header contains BUNIT='COUNTS', the data is divided by EXPTIME
      to convert from counts to count rate (counts/s). Otherwise, the data is used as-is
  - Default exposure time of 700 seconds is used if EXPTIME header keyword is missing
  '''

  fluxstack = []
  photflams = []
  for file in filelist:
    with fits.open(file) as hdul:
        header = hdul[1].header
        data = hdul[1].data
        exptime = header.get('EXPTIME', 700) #s
        photflam = header.get('PHOTFLAM')
        photzpt = header.get('PHOTZPT')

        if 'BUNIT' in header and header['BUNIT'].strip().upper() == 'COUNTS':
          data_conv = data/exptime
        else:
          data_conv = data

        flux_data = data_conv
        flux_data = np.clip(flux_data, 0, None)
        fluxstack.append(flux_data)
        photflams.append(photflam)

  return np.median(fluxstack, axis = 0), np.median(photflams), photzpt

#####

def flux_to_apparent(value, zeropoint):
  '''
  Function takes given flux values and converts to  apparent magnitude,
  using a specified zeropoint value

  Input parameters:
  value : dtype = float
    The flux value to be conveted. Must be positive, and in same units
    as zeropoint value used
  zeropoint : dtype = float
    photometric zeropoint for the magnitude system

  Returns:
  app : dtype = float
    Apparent magnitude value, calculated from input flux and zeropoint

  Notes:
  - If the flux value given is negative, the function will return a magnitude of 99 (a non detection)
  '''
  if value > 0:
    app = -2.5 * np.log10(value) + zeropoint

  else:
    app = 99

  return app

#####

def ellipticitycheck(p1, p2):
  '''
  Function calculates the ellipticity from two axis measurements

  Input parameters:
  p1 : dtype = float
    First axis value
  p2 : dtype = float
    Second axis value

  Returns:
  f : dtype = float
    Value of ellipticity, which will be between 0 and 1
    - f = 0 indicates a perfect circle (p1 = p2)
    - f = 1 indicates an ellipse (one axis > other axis)

  Notes:
  Function automatically determines major and minor axis' by taking minima
  and maxima of both points
  '''
  f = 1 - min(p1, p2) / max(p1, p2)
  return f

#####

def square_aperture(data, center, box_size):
    """
    Function performs a square cutout of data

    Input parameters:

    data : dtype = array
      input data, to be apertured
    center : dtype = list, or tuple
      the center of square aperture, in y-x order
    box_size : dtype = int
      box size in pixels

    Returns:
    res : dtype = array
      square cutout of data
    """
    y, x = map(int, center)

    if box_size % 2 == 0:
        box_size += 1

    half = box_size // 2

    sly = slice(y - half, y + half + 1)
    slx = slice(x - half, x + half + 1)

    return data[sly, slx]

#####

def circular_aperture(data, center, radius):
    """
    Function performs a circular cutout of data

    Input parameters:

    data : dtype = array
      input data, to be apertured
    center : dtype = list, or tuple
      the center of square aperture, in y-x order
    radius : dtype = int
      circle radius in pixels

    Returns:
    cutout : dtype = array
      circular cutout of data
    mask : dtype = array
      circular mask used for cutout
    """
    box_size = int(np.ceil(2 * radius))
    cutout = square_aperture(data, center, box_size)

    ys, xs = np.indices(cutout.shape)
    xc = cutout.shape[1] // 2
    yc = cutout.shape[0] // 2

    y0 = yc - box_size // 2
    x0 = xc - box_size // 2

    x_coords = xs + x0
    y_coords = ys + y0

    distance = np.sqrt((x_coords - xc) ** 2 + (y_coords - yc) ** 2)

    mask = distance <= radius

    return cutout, mask

#####

def multicoord(yc, xc, data, photflam, photzpt):
  '''
  Function extracts flux and magnitude for sources using circular aperture
  photometry with local background subtraction, analysed over an 8-pixel radius.

  Input parameters:
  yc : dtype = list
    list of y-coordinates for source center
  xc : dtype = list
    list of x-coordinates for source center
    must be same length as yc
  data : dtype = array
    2d image array containing flux data
  photflam : dtype = float
    calibration value, to be used in conversion of flux data to physical data
  photzpt : dtype = float
    calibration value, to be used in conversion of flux data to apparent magnitude

  Returns:
  srcmagvals : dtype = list of floats
    Apparent magnitudes for each source after background subtraction and
    flux calibration.
  srcfluxvals : dtype = list of floats
    Net instrumental flux values (counts/s) for each source after background
    subtraction. Negative values are replaced with 1e-20 (this is a redundant check, but I added it for completeness).
  cutouts : dtype = list of arrays
    Image cutouts centered on each source, extracted by the circular aperture.
  masks : list of arrays
    Boolean masks indicating pixels within the aperture (True) for each source.
  background : dtype = list of floats
    Background flux per pixel for each source, estimated from the region
    outside the mask.

  Notes:
  Photometry procedure for each source:
    1. Extract 8-pixel radius circular mask, center = (y, x)
    2. Estimate background from pixels outside the mask
    3. Calculate background per pixel and subtract from total mask flux
    4. Convert net flux to physical units using PHOTFLAM
    5. Convert calibrated flux to apparent magnitude using PHOTZPT
  '''
  cutouts = []
  masks = []
  srcmagvals = []
  background = []
  srcfluxvals = []

  for y, x in zip(yc, xc):
    cutout, mask = circular_aperture(data, center = (y, x), radius = 8)
    invmask = ~mask

    backpix = np.sum(invmask)
    backflux = np.sum(cutout[invmask])

    if backflux != 0:
      backperpix = backflux/backpix
    else:
      backperpix = 0

    sourcepix = np.sum(mask)
    totalflux = cutout[mask].sum()
    totalback = backperpix * sourcepix


    netflux = totalflux - totalback
    if netflux < 0:
      netflux = 1e-20

    srcfluxvals.append(netflux)

    fluxphysical = netflux * photflam
    final = flux_to_apparent(fluxphysical, photzpt)
    srcmagvals.append(final)

    cutouts.append(cutout)
    masks.append(mask)
    background.append(backperpix)
  return srcmagvals, srcfluxvals, cutouts, masks, background

#####

def remove_close(xcoord, ycoord, mags, minsep = 1.5):
  '''
  Function removes sources which are too close together in an image,
  by calculating distance between given coords, and seeing if said distance is above
  or below the given threhsold

  Input parameters:
  xcoord : dtype = list
    List of x coordinates of detected sources
  ycoord : dtype = list
    List of y coordinates of detected sources
  mags : dtype = list
    list of magnitudes for corresponding sources
  minsep : dtype = float, optional
    Minimum separation between sources, defaults to 1.5 pixels

  xcoord, ycoord and mags must have same length

  Returns:
  sepcoords : dtype = list of tuples
    list of coordinates for sources sufficiently apart
  sepmags : dtype = list of tuples
    list of magnitudes, corresponding to sources sufficiently apart

  Notes:
  - When two sources are closer than minsep, the first source is retained,
    while the second is rejected
  - Euclidean distance is used
  '''

  sepcoords = []
  sepmags = []
  for idx in range(len(xcoord)):

    y, x = ycoord[idx], xcoord[idx]

    close = False
    for yy, xx in sepcoords:
      dist = np.sqrt((x - xx)**2 + (y - yy)**2)
      if dist < minsep:
        close = True
        break

    if close == False:
      sepcoords.append((y, x))
      sepmags.append(mags[idx])

  return sepcoords, sepmags

#####

def keep_close(xcoord1, ycoord1, mag1, xcoord2, ycoord2, mag2):
  '''
  Function decides whether sources between filters are kept,
  by calculating distance between given coords across filters,
  and checking if said distance is above or below the given threhsold
  (i.e., if the sources are in the same position across filters).

  Input parameters:
  xcoord1 : dtype = list
    List of x coordinates of detected sources in first filter
  ycoord1 : dtype = list
    List of y coordinates of detected sources in first filter
  xcoord2 : dtype = list
    List of x coordinates of detected sources in second filter
  ycoord2 : dtype = list
    List of y coordinates of detected sources in second filter
  mag1 : dtype = list
    list of magnitudes for corresponding sources for first filter
  mag2 : : dtype = list
    list of magnitudes for corresponding sources for first filter

  xcoord1, ycoord1, xcoord2, ycoord2, mag1 and mag2 must have same length

  Returns:
  closecoords : dtype = list of tuples
    list of coordinates for sources close together across filters
  closemags : dtype = list of tuples
    list of magnitudes, corresponding to sources sufficiently close

  Notes:
  - When two sources are closer than minsep, the first source is retained,
    while the second is rejected
  - Euclidean distance is used
  '''
  closecoords = []
  closemags = []
  for idx in range(len(xcoord1)):

    y1, x1 = ycoord1[idx], xcoord1[idx]
    mags1 = mag1[idx]

    close = False
    for iidx in range(len(xcoord2)):
      y2, x2 = ycoord2[iidx], xcoord2[iidx]
      mags2 = mag2[iidx]
      dist = np.sqrt((x1 - x2)**2 + (y1 - y2)**2)

      if dist < 2:
        close = True
        break

    if close == True:
      closecoords.append((y1, x1))
      closemags.append((mags1, mags2))

  return closecoords, closemags

#####

def gauss2d(coordinates, amp, x0, y0, sigma_x, sigma_y):
  '''
  Calculate a 2D Gaussian function over given coordinates

  Input parameters:
  coordinates : dtype = tuple
    Tuple containing (y, x) coordinates where Gaussian is created
  amp : dtype = float
    Amplitude (peak height) of the Gaussian at the center (x0, y0).
  x0 : dtype = float
    X-coordinate of the Gaussian center.
  y0 : dtype = float
    Y-coordinate of the Gaussian center.
  sigma_x : dtype = float
    Standard deviation of the Gaussian in the x-direction. Controls the
    width along the x-axis.
  sigma_y : dtype = float
    Standard deviation of the Gaussian in the y-direction. Controls the
    width along the y-axis.

  Returns:
  Calculated Gaussian : dtype = np.array
    A ravelled 1D array
  '''
  y, x = coordinates

  return amp * np.exp(-(((x - x0)**2) / (2 * sigma_x**2) + ((y - y0)**2) / (2 * sigma_y**2))).ravel()

#####

def fitgauss2d(data, p0 = None):
  '''
  Function fits a 2D Gaussian function to imagery data

  Input parameters:
  data : dtype = array
    Image data to be fitted. Should contain a single source, centered approximately in the frame.
  p0 : dtype = list, optional
    Initial guess for fit parameters in the order:
    [amplitude, x_center, y_center, sigma_x, sigma_y]
    If None, initial values are generated via:
    - amplitude: max(data) - median(data)
    - x_center, y_center: center of the data
    - sigma_x, sigma_y: 0.8 pixels

  Returns:
  popt : dtype = array
    Optimised fit parameters for Gaussian fit:
      [amplitude, x_center, y_center, sigma_x, sigma_y]

  pcov : dtype = array
    Covariance matrix of fit parameters

  Notes:
  Use of bounds:
  - All params must be positive
  - Centroid positions must be reasonable, relative to input data
  - Sigma values are constrained to produce reasonable PSFs

  If the fit fails to converge, it will automatically retry, with max_fev = 2000

  Function used previously defined gauss2d function.
  '''
  if p0 is None:
    a0 = data.max() - np.median(data)
    x0 = data.shape[1] / 2
    y0 = data.shape[0] / 2

    sigx0 = 0.8
    sigy0 = 0.8
    p0 = [a0, x0, y0, sigx0, sigy0]
  y, x = np.indices(data.shape)

  bounds = ([0, 0, 0, 0.2, 0.2], [np.inf, data.shape[1], data.shape[0], 3, 3])
  try:
    popt, pcov = curve_fit(gauss2d, (y, x), data.ravel(), p0 = p0, bounds = bounds)

    return popt, pcov
  except:
    popt, pcov = curve_fit(gauss2d, (y, x), data.ravel(), p0 = p0, bounds = bounds, maxfev = 2000)

    return popt, pcov

#####

def make_info_table(table, popt, f, SNR, mags, coords, stat, reason=None):
  '''
  Function populates an inputted table with given photometric and Gaussian fit parameters

  Input parameters:
  table : dtype = dictionary of lists
    Dictionary containing lists for each column of the table. Must have keys:
    'amplitude', 'xcentroid', 'ycentroid', 'sigx', 'sigy', 'ellipticity',
    'SNR', 'Coord', 'Magnitudes', and optionally 'Reason'.
  popt : dtype = array
    Optimised Gaussian fit parameters. Order:
      [Amplitude, xcenter, ycenter, sigmax, sigmay]
  f : dtype = float
    Measure of ellipticity of source
  SNR : dtype = float
    Calculated signal-to-noise ratio
  mags : dtype = float
    Magnitude of the source
  coords: dtype = tuple
    Pixel coords of the source in question
  stat : dtype = str
    Dictates whether or not the source in question is appended to the
    accepted or rejected candidates table
  Reason : dtype = str, optional
    Reason for rejection if stat = 'Reject'
    Default is None otherwise

  Returns:
  None
  Function modifies the given table, appending the values to the relevant columns
  '''
  amp, xc, yc, sigx, sigy = popt
  table['amplitude'].append(amp)
  table['xcentroid'].append(xc)
  table['ycentroid'].append(yc)
  table['sigx'].append(sigx)
  table['sigy'].append(sigy)
  table['ellipticity'].append(f)
  table['SNR'].append(SNR)
  table['Coord'].append(coords)
  table['Magnitudes'].append(mags)

  if stat == 'Reject':
    table['Reason'].append(reason)

#####

def HRgenerator(mags1, mags2):
  '''
  Function generates a Hertzprung Russel diagram using input magnitude data.
  Specifically, the function plots colour (F336 - F555) vs F555 Magnitude

  Input parameters:
  mags1 : dtype = list of tuples
    List of 1st filter magnitudes in the form of a tuple

  mags2 : dtype = list of tuples
    List of 2nd filter magnitudes in the form of a tuple

  Returns:
  None

  Function produces a matplotlib plot, does not return any values

  Notes:
  - They Y axis is inverted, as is standard convention, with brighter stars at the top
  - Stars are plotted as black points
  - Function uses previously defined apparenttoabs function to convert apparent magnitudes
    to absolute magnitudes
  '''
  mag555 = []
  mag336 = []
  colours = []

  for mag in mags1:
    mag336.append(mag)

  for mag in mags2:
    mag555.append(mag)

  for m336, m555 in zip(mag336, mag555):
    colour = m336 - m555
    colours.append(colour)

  plt.figure(figsize = (8, 8))
  plt.scatter(colours, mag555, s = 10, c = 'black')

  plt.axhline(y = 21.5, linestyle = '--', color = 'green', label = 'Below line: Main Sequence')
  plt.axhline(y = 20.5, linestyle = '--', color = 'purple', label = 'Below line: Subgiant branch')
  plt.axhline(y = 18.7, linestyle = '--', color = 'blue', label = 'Below line: Asymptotic Giant Branch')

  plt.plot([-3, -1.75], [26, 21.5], color='red', linestyle='--', linewidth=1.5, label='Main Sequence', alpha=0.7)

  plt.xlabel('F336 - F555')
  plt.ylabel('F555 Apparent Magnitude')
  plt.title('HR diagram for NGC 1261')
  plt.gca().invert_yaxis()
  plt.grid(True)
  plt.legend(loc = 'best')

  plt.savefig('hrdiagramofNGC1261.png')
  print('HR diagram saved.')

  plt.close()

def normalize_image(x):
    """
    Function normalises data to range (0 - 1), clipping the original data
    between defined percentiles
    Input parameters:
    x: dtype = array
        Data to be normalised
    Returns:
    norm: the normalized data in the range [0, 1]
    """
    vmin, vmax = np.percentile(x, [1, 99])
    norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)

    return norm(x)

def transform_by(x):
    """
    transform/stretch the data by a named function,
    normalizing the data using given percentiles

    Input parameters:
    x: dtype = array
        data to be transformed

    Returns:
    x_norm : dtype = array
        the stretched data
    """
    x_norm = normalize_image(x)

    return np.sqrt(x_norm)

#####

def master(filelistlist, filternames, val, dec, decc = 'No'):
  '''
  Master function, which calls together all previously defined functions, to cosmic ray removal, then
  star candidcacy checks and photometry, and producing a HR diagram (and 3 optional files containing Accepeted
  and rejected sources based on quality, and the final catalogue of star candidates)

  Processes two filters with the following process:
  i) median combine FITS files to remove cosmic rays
  ii) Identify initial candidate sources via percentile thresholding
  iii) Filter sources by edge proximity 
  iv) Aperture background with local background subtraction
  v) Signal-to-Noise ratio filtering
  vi) Fit 2D gaussian Point Spread Functions, filters on ellipticity, size and SNR
  vii) Cross validate source positions across filters
  viii) generate HR diagram of successful candidates

  Input parameters:
  filelistlist : dtype = list of lists of strings
    A list containing two sublists
    - filelist1 = List of file paths for FITS files matching filestring1
    - filelist2 = List of file paths for FITS files matching filestring2

  filternames : dtype = list of strings
    A list containing the names of the HST filters used, for printing to know which filter is being analysed

  val : dtype = int
    Percentile value used for intensity thresholding at the very start. If the pixels are above this threshold
    they are considered potential candidates. Higher values result in fewer, brighter stars from the start.
    This parameter has a massive impact on runtime. Typical value = 92 

  dec : dtype = string
    A decision flag that dictates whether or not a certain file is saved as a csv file
    Any other value will result in no file saving


  Returns:
  Finalcands : dtype = pandas dataframe or None
    Final catalogue containing successful candidate information:
    - 'Coords' : x, y pixel coordinates of candidate
    - 'Mags F555' : Apparent magnitude of source in F555 filter
    - 'Mags F336' : Apparent magnitude of source in F336 filter
    - 'Aperture Radius' : radius of aperture used in circular photometry for generating cutout

  or returns empty dataframe if the code fails somewhere along the way.

  Notes:
  If dec is set to 'savecsv' - Output files generated:
  'Accepeted (Quality) - filter: {filtername}.csv' - Sources which passed quality checks
  'Rejected (Quality) - filter: {filtername}.csv' - Sources which failed quality checks

  'Finalcatalogue.csv' - Cross matched sources detected in both filters

  Function prints status messages containing how many candidates are left after previous filter step

  Function final output is a HR diagram plotting (F336W - F555W) colour vs F555 Magnitude
  '''
  results = []
  combined = []
  for filterind, filelist in enumerate(filelistlist):
    print(f'Analysing filter: {filternames[filterind]}') # add print to see which filter is being analysed
    imgstack, photflam, photzpt = rayremovingfluxconverter(filelist) # perform cosmic ray removal using function
    combined.append(imgstack)
    thresh = np.percentile(imgstack, val) # using percentile thresholding 
    flag = imgstack > thresh # only take values greater than the threshold
    ys, xs = np.where(flag)
    print(f'  After percentile thresholding, there are now {len(xs)} candidates')

    acctable = {
         'Coord' : [],
         'amplitude':[],
         'xcentroid' : [],
         'ycentroid' : [],
         'sigx' : [],
         'sigy': [],
         'ellipticity': [],
         'SNR': [],
         'Magnitudes':[]
      }

    rejtable = {
         'Coord' : [],
         'amplitude':[],
         'xcentroid' : [],
         'ycentroid' : [],
         'sigx' : [],
         'sigy': [],
         'ellipticity': [],
         'SNR': [],
         'Magnitudes': [],
         'Reason': [] # additional column added to rejects table, to which the reason for rejection will be appended
    }

    candsflux = (imgstack[flag] > 0) # filter out any negative pixels that might have somehow survived percentile thresholding 

    # edge filtering - never trust sources at the edge of an image
    if np.any(candsflux):
      ysin = ys[candsflux] < 775
      xsin = xs[candsflux] < 775
      valid = xsin & ysin

      validy = ys[candsflux][valid]
      validx = xs[candsflux][valid]

      print(f'  After edge filtering, there are now {len(validx)} candidates')
    else:
      print('Thats not good, there are no sources')


    # Early crude filtering
    pixels = imgstack[validy, validx] # extract pixel values at valid coordinates
    magsrough = -pixels # create a rough check of magnitude for early proximity filtering only

    sepcoords, sepmags = remove_close(validx, validy, magsrough.tolist()) # remove any sources too close together in the image

    print(f'  After removing sources too close together, there are now {len(sepcoords)} candidates')

    xssep = [x for y, x in sepcoords] # reseperate coords into x and y for use in multicoord function
    yssep = [y for y, x in sepcoords]
    
    # photometry performed using multicoord function
    magssep, fluxsep, cutoutssep, maskssep, backgroundsep = multicoord(xssep, yssep, imgstack, photflam, photzpt)

    # performing SNR filtering

    # Generative AI was used in the below 2 lines. I did not realise I had to make arrays of
    # the flux and background values for further analysis, and so the code continuously failed.
    # GenAI (specifically GPT-4-Turbo) suggested making arrays of these values would solve the issue.
    fluxseparray = np.array(fluxsep)
    backgroundseparray = np.array(backgroundsep)

    SNRs = np.where(backgroundseparray != 0, fluxseparray/ (backgroundseparray), 0) # calculate SNR, avoiding divison by 0

    SNRval = (SNRs >= 3)

    print(f'  After SNR Filtering, there are now {np.sum(SNRval)} candidates')

    indstofit = np.where(SNRval)[0] # finding indices of sources which have been successful so far

    # filtering based on PSF characteristics
    for idx in indstofit:
      cutout = cutoutssep[idx] # aperture cutout
      y, x = sepcoords[idx] # pixel coordinates
      mag = magssep[idx] # source magnitude
      backgroundlvl = backgroundsep[idx] # local background level

      cutoutbckg = cutout - backgroundlvl # performing background subtraction
      cutoutbckg = np.clip(cutoutbckg, 0, None) # remove any negative pixels

      popt, pcov = fitgauss2d(cutoutbckg) # extract Gaussian params and covariant matrix

      if popt is None: # if the fit failed, skip the source
        continue

      sigx, sigy = popt[3], popt[4]

      sourceSNR = SNRs[idx] # SNR value of source

      f = ellipticitycheck(sigx, sigy) # calculate ellipticity of source 

      snrok = (sourceSNR > 3) & (sourceSNR != 0) # check if SNR meets requirements for significance

      ellipseok = f < 0.4 # check if ellipticity is circular (roughly)

      # Threshhold values below obtained from: 
      #    https://hst-docs.stsci.edu/wfc3ihb/chapter-6-uvis-imaging-with-wfc3/6-6-uvis-optical-performance#id-6.6UVISOpticalPerformance-6.6.16.6.1PSFWidthandSharpness

      sizeok =  (0.5 < sigx < 1.6) & (0.5 < sigy < 1.6) # check if source size meets requirements 

      if ellipseok & sizeok & snrok: # use boolean logic to see if all filters passed
         make_info_table(acctable, popt, f, sourceSNR, mag, (int(x), int(y)), 'Success') # create table of successful candidates
      else:
          reason = []
          if not ellipseok: # append specific reason for failing to rejects table
            reason.append('Not elliptical')
          if not sizeok:
            reason.append('Bad sizing')
          if not snrok:
            reason.append('Signal likely to be noise, or the background was bad')

          make_info_table(rejtable, popt, f, sourceSNR, mag, (int(x), int(y)), 'Reject', reason = ','.join(reason))
    # make dataframes of the tables
    dfaccepted = pd.DataFrame(acctable)
    dfrejected = pd.DataFrame(rejtable)

    if dec == 'savecsv': # condition for saving of the files
      dfaccepted.to_csv(f'Accepted (Quality) - filter: {filternames[filterind]}.csv', index=False)
      dfrejected.to_csv(f'Rejected (Quality) - filter: {filternames[filterind]}.csv', index=False)

    print(f'  After Aperture characteristic thresholding, there are now {len(dfaccepted)} candidates')
    print()

    results.append((dfaccepted, dfrejected)) 

  if len(results) >= 2: # the tuple should have both tables
    # GenAI (specifically GPT-4-Turbo) was used here. Initially, I did not realise how to access the
    # specific parts of 'results', to access the accepeted magnitude values, as I was using incorrect indexing.
    # GenAI gave the solution:
    F555 = results[0][0] # accessing the F555 data of the accepeted sources
    F336 = results[1][0] # accessing the F336 data of the accepeted sources

    if (len(F555) > 0) and (len(F336) > 0): # ensuring both lists are populated
      F555coords = F555.Coord
      F336coords = F336.Coord

      F555mags = F555.Magnitudes
      F336mags = F336.Magnitudes

      F555ycs = []
      F555xcs = []

      F336ycs = []
      F336xcs = []

      for coords in F555coords:
        F555ycs.append(coords[0])
        F555xcs.append(coords[1])

      for coords in F336coords:
        F336ycs.append(coords[0])
        F336xcs.append(coords[1])

      # check if sources between 2 filters are close enough together
      finalcandcoords, finalcandmags = keep_close(F555xcs, F555ycs, F555mags, F336xcs, F336ycs, F336mags) 

      finalcands = {
      'Coords':[],
      'Mags F555': [],
      'Mags F336': [],
      'Aperture Radius': [],
      }

      finalcands['Coords'] = finalcandcoords

      for magg in finalcandmags:
        finalcands['Mags F336'].append(magg[1])
        finalcands['Mags F555'].append(magg[0])
      
      finalcands['Aperture Radius'] = 8

      Finalcands =  pd.DataFrame(finalcands)

      if dec == 'savecsv':
        Finalcands.to_csv(f'Finalcatalogue.csv', index=False)

      mags336 = Finalcands['Mags F336'] # access magnitudes for filter from Finalcands dataframe
      mags555 = Finalcands['Mags F555']
      print(f'  After cross checking between F555 and F336 for pixel proximity, there are now {len(mags336)} candidates')

      HRgenerator(mags336, mags555) #use HR generator function to produce HR diagram

      if decc == 'Plot':
        fig, ax = plt.subplots(figsize=(8,8))
        norm = transform_by(combined[0])
        im = ax.imshow(norm, origin='lower')
        
        crds = Finalcands['Coords']

        ycrds = [crd[0] for crd in crds]
        xcrds = [crd[1] for crd in crds]

        ax.scatter(xcrds, ycrds, marker = 'x', color = 'red', label = 'Successful candidates')

        ax.axvline(x=775, color='black', linestyle='--')
        ax.axvline(x=25, color='black', linestyle='--')
        ax.axhline(y=775, color='black', linestyle='--')
        ax.axhline(y=25, color='black', linestyle='--')

        ax.set_title('Final Source Candidates')
        ax.legend(loc='best')
        ax.grid(True)

        plt.savefig('SuccessfulCandPos.png')

        # GenAI was used here. Everytime the script was ran (in Google Colab, so a Jupyter environment),
        # both generated plots would be plotted repeatedly. I could not figure out why, GenAI provided the solution of using
        # plt.close (which I should've been using anyway, as it's good practice) to prevent repeated plotting
        plt.close()

      return Finalcands

    else:
      print('Something is wrong, one or more filter results are empty:')
      print(f'{len(F555)}, {len(F336)}')

  else:
    print('Filter Missing')
    return pd.DataFrame(finalcands) # return empty data frame if the code fails to run

#####

if __name__ == "__main__":
  import sys
  import os

  # defaults for running the function
  percentile = 93
  savecsv = True
  plot = False


  if len(sys.argv) > 1:
    data_folder = sys.argv[1]
  else:
    data_folder = 'data' # defaults to name of data folder if no input name given

  if len(sys.argv) > 2:
    # GenAI was used here. Upon running the script at the simplest level (i.e., %run myscript.py in Jupyter environment),
    # there would be an error thrown attempting to turn sys.argv[2] to a float. This occured when any command line flags were used without
    # a percentile threshold value, the code would try to convert the string to a float. GenAI suggested the use of a try-except-pass block, which
    # solved the issue. 
    try: 
      percentile = float(sys.argv[2]) # percentile value option, as this drastically impacts the runtime
    except:
      pass 

  if 'nocsv' in sys.argv:
    savecsv = False # option to not save csv files

  if 'Plot' in sys.argv:
    plot = True

  if not os.path.exists(data_folder):
    print(f"ERROR: Could not find data folder '{data_folder}'")
    print(f"\nPlease ensure the 'data' folder is in the same directory as this script.")
    print(f"The folder should contain F555W/ and F336W/ subdirectories with FITS files.")
    sys.exit(1)
  # check to ensure files are present as expected

  F555 = glob(f"{data_folder}/F555W/*.fits")
  F336 = glob(f"{data_folder}/F336W/*.fits")

  print('There should be 3 files in each filter folder:')
  print(f"Found: {len(F555)} F555W | {len(F336)} F336W")

  print('\nStarting Analysis:')

  try:
    filelistlist = fileintaker('F555W', 'F336W') # intake FITS files using defined function
    filternames = ['F555W', 'F336W']

    savemode = 'savecsv' if savecsv else 'nocsv'
    plotmode = 'Plot' if plot else 'No' 
    result = master(filelistlist, filternames, percentile, savemode, plotmode) # option for CSV saving

    if savecsv:
      print('Saving CSV files...')
    else:
      print('CSV files will not be saved.')

    if plot:
      print('Saved: SuccessfulCandPos.png')
    else:
      print('Candidate positions will not be displayed.')

    print('DONE')

  except Exception as e:
    print(f'ERROR: {e}')
    sys.exit(1)
