#!/usr/bin/python
'''
Kernel Density Estimation Tools

Functions aid in creating a point density 
distribution using Kernel Density Estimation.
An optimized kernel bandwidth  is calculated in R
using the SAMSE method (Tarn Duong, 2007 
http://www.jstatsoft.org/v21/i07/paper) in 
samse_bandwidth().

Functions reproject data, plot data, and evaluate
point density on three planet surfaces.

Code has not been extensively tested, & comes with no
guarantees. Contact jarichardson@mail.usf.edu if bugs
are found.

AUTHOR: Jacob Richardson (github.com/jarichardson)

REQUIREMENTS: GMT, R with KS library, GDAL, PROJ.4,
  and pylibs: pyproj, numpy, os, scipy, osgeo
'''
import pyproj
import numpy
import os
import scipy.linalg as linalg
from scipy.stats import norm
from osgeo import gdal,osr

#REPROJECT
def reproject(llcoors,planet="earth",utmzone=-999,clon=-999,inverse=False):
	'''
	Reprojects long-lat data into transverse mercator coordinates
	with units of meters. Optionally, set inverse=True to change
	transverse mercator to long-lat.
	
	Input should be numpy 2-col long, lat array.
	Output will be numpy 2-col easting, northing array.	
	
	Planet options: 'earth', 'venus', or 'mars'
	Earth requires a valid UTM zone
	Venus and Mars require a valid center longitude of the dataset
	
	Earth Transverse Mercator fit to WGS84 datum
	Venus Transverse Mercator fit to Spheriod of radius 6051800 m
	Mars Transverse Mercator fit to Spheriod of radius 3396190 m
	'''
	
	if planet == "earth":
		if (utmzone<1) or (utmzone>60):
			print "error in reproject: utm zone not set correctly (1<=utmzone<=60)"
			return 0
		TransMerc = pyproj.Proj('+proj=utm +datum=WGS84 +zone='+str(utmzone))
	elif planet == "venus":
		if (clon<-360) or (clon>360):
			print "error in reproject: center longitude not set correctly (-360<=clon<=360)"
			return 0
		TransMerc = pyproj.Proj('+proj=tmerc +lat_0=0 +lon_0='+str(clon)+' +k=0.9996 +x_0=0 +y_0=0 +a=6051800 +b=6051800 +units=m +no_defs')
	elif planet == "mars":
		if (clon<-360) or (clon>360):
			print "error in reproject: center longitude not set correctly (-360<=clon<=360)"
			return 0
		TransMerc = pyproj.Proj('+proj=tmerc +lat_0=0 +lon_0='+str(clon)+' +k=0.9996 +x_0=0 +y_0=0 +a=3396190 +b=3396190 +units=m +no_defs')
	else:
		print "error in reproject: planet not earth, venus, or mars."
		return 0
	
	if inverse==True:
		reproj = TransMerc(llcoors)
		mcoors = numpy.transpose(reproj)
	else:
		reproj = TransMerc(llcoors)
		mcoors = numpy.transpose(reproj)
	
	return mcoors

#FIND A BUFFERED RANGE OF N-DIMENSIONAL DATA
def rangeBuffer(coords,B=0):
	'''
	Creates a buffer of B% [default 0%, no buffer] around 
	N-dimensional data. Input should be a numpy array.
	Output will be 2xN array, with min, max of each dimension in
	columns 1 and 2, respectively.
	
	ex: 
	data         range output
	[[1,5],      [[0,2],
	 [2,5],  =>   [4,9]]
	 [0,4],
	 [1,9]]
	'''
	extents = numpy.ones([numpy.shape(coords)[1],2])
	for dim in range(numpy.shape(coords)[1]):
		data = coords[:,dim]
		dataRange = data.max() - data.min()
		bufsize = dataRange*(B/100.0)
		extents[dim,0] = data.min() - bufsize
		extents[dim,1] = data.max() + bufsize
		
	return extents
	
def samse_bandwidth(coords):
	'''
	Evaluates the SAMSE Kernel in R using a coordinate list (coords).
	Returns 2x2 bandwidth covariance matrix.
	Requires: R, KS library in R.
	'''
	
	bandwidthfile='tmpbdR.out'
	datafile     ='tmpcrs.out'
	#Writes the data to a file for R
	numpy.savetxt(datafile,coords)
	
	#Writes the batch file that R will use
	with open('samse_batch.r','w') as f:
		f.write('library(ks)\n')
		f.write('data<-read.table("'+datafile+'")\n')
		f.write('bd <- Hpi(x=data,nstage=2,pilot="samse",pre="sphere")\n')
		f.write('sink("'+bandwidthfile+'")\n')
		f.write('show(bd)\n')
		f.write('sink()')
	
	#command to run the batch file
	os.system('R CMD BATCH samse_batch.r')
	
	#Extract the bandwidth matrix from the bandwidth txt file
	bandwidth = numpy.loadtxt(bandwidthfile,skiprows=1,usecols=(1,2))
	
	#remove all these temporary files
	os.remove('samse_batch.r')
	os.remove('samse_batch.r.Rout')
	os.remove(bandwidthfile)
	os.remove(datafile)
	
	return bandwidth


def KDE(bd,coors,ranges,spacings):
	'''
	Estimates point density using:
	bd       - a kernel bandwidth (2x2 covariance	matrix)
	coors    - 2xN list of coordinates for N points.
	ranges   - a 2x2 [[W,E],[S,N]] array
	spacings - a 1x2 [X-resolution,Y-resolution] array
	
	Outputs X,Y,D: Eastings, Northings, and Densities in a Meshgrid
	format (i.e. X will be tiled, Y will be repeated)
	'''
	
	ventct = len(coors)
	detH = linalg.det(linalg.sqrtm(bd)) #determinate sqrt bandwidth
	invH = linalg.inv(linalg.sqrtm(bd)) #inverse sqrt bandwidth
	
	constant = 2.0*numpy.pi*detH*ventct

	#define map grid
	x = numpy.arange(ranges[0][0],(ranges[0][1]+spacings[0]),spacings[0])
	y = numpy.arange(ranges[1][0],(ranges[1][1]+spacings[1]),spacings[1])
	X,Y = numpy.meshgrid(x,y)	#X and Y are now tiled to grid
	D = numpy.zeros(numpy.shape(X)) #Density Grid
	dist = numpy.zeros(numpy.shape(X)) #distance matrix grid
	
	for v in coors:
		for i,e in enumerate(x):
			for j,n in enumerate(y):
				dx = e-v[0]
				dy = n-v[1]
				dxdy = numpy.dot(invH,numpy.array[[dx],[dy]])
				dist[i][j] = numpy.dot(numpy.transpose(dxdy),dxdy)[0][0]
		D += numpy.exp(-0.5 * dist)
	
	D /= constant #normalize
	
	return X,Y,D
	

def ellipseGen(bd,eps=False,epsfilename='bandwidth_ellipse.eps'):
	'''
	Identifies the major and minor axes directions and standard
	deviations of a Gaussian ellipse defined by a 2x2 covariance
	matrix. Precision is to the nearest degree.
	
	Prints out solution, and optionally uses GMT to draw the ellipse
	to epsfilename, if eps=True.
	
	Outputs major-axis direction, major-axis standard deviation, and
	        minor-axis standard-deviation.
	'''
	detH = linalg.det(linalg.sqrtm(bd)) #determinate sqrt bandwidth
	invH = linalg.inv(linalg.sqrtm(bd)) #inverse sqrt bandwidth
	constant = 2.0*numpy.pi*detH*ventct
	
	radius = 20
	angles = numpy.arange(0,numpy.pi,(numpy.pi/180.0))
	D = numpy.zeros(len(angles))
	
	for i,phi in enumerate(angles):
		dx = radius*cos(phi)
		dy = radius*sin(phi)
		dxdy = numpy.dot(invH,numpy.array[[dx],[dy]])
		dist = numpy.dot(numpy.transpose(dxdy),dxdy)[0][0]
		D[i] = numpy.exp(-0.5*dist)/(detH*constant)
	
	maxaz = az[numpy.where(D=max(D))]
	minaz = az[numpy.where(D=min(D))]
	
	#Calculate Density at vent location
	dxdy = numpy.dot(invH,numpy.array[[0],[0]])
	dist = numpy.dot(numpy.transpose(dxdy),dxdy)[0][0]
	ventD = numpy.exp(-0.5*dist)/(detH*constant)
	
	#Calculate standard deviations
	#For the major axis
	majsd = (10*(2**0.5))/(-1*log(max(D)/ventD))**0.5 #(radius=20 units)
	majdir = 90-numpy.degrees(maxaz) #Gives direction from North. East is +
	#For the minor axis
	minsd = (10*(2**0.5))/(-1*log(min(D)/ventD))**0.5
	mindir = 90-numpy.degrees(minaz)
		
	print 'Quick check:'
	print ('Density at vent             = %0.3e' % ventD)
	print ('Theoretical density at vent = %0.3e' % (1/(2*numpy.pi*majsd*minsd)))

	#Print out the results
	print '\nBandwidth Ellipse Information'
	print 'major axis:'
	print ('  degrees from north - %0.1f' % majdir)
	print ('  standard deviation - %0.4f' % majsd)
	print 'minor axis:'
	print ('  degrees from north - %0.1f' % mindir)
	print ('  standard deviation - %0.4f' % minsd)
	
	if eps==True:
		majaxis = 2*majsd
		minaxis = 2*minsd
		
		with open('ellipseGMT.xy','w') as f:
			f.write('0\t0\t%0.0f\t%0.4f\t%0.4f' % majdir, majaxis, minaxis)
		os.system('psxy ellipseGMT.xy -SE -Wblack -JX6i -R-%0.4f/%0.4f/-%0.4f/%0.4f -Ba%0.4f -K >'+epsfilename % majaxis,majaxis,majaxis,majaxis,majsd)
		os.remove('ellipseGMT.xy')			
		print ('\nPlotted ellipse at '+epsfilename)
		
		return majdir,majsd,minsd

def contourBySigma(Z,sigmas=[0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0],gridspacings=[1,1]):
	'''
	Identifies Density contour levels given Sigma thresholds.
	Contours define areas within which lay X% of the total field density
	e.g.: Within the Sigma-2.0 contour lies 95.4% of total field density.
	      The density value of each contour decreases with increasing
	      sigma.
	Requires a numpy array of density values (Z), with any shape.
	Optionally, provide a list of requested sigma thresholds, and the
	grid size as a 2 item list, to normalize the density.
	
	Outputs a dictionary of {sigma-level: density value}. If sigma-levels
	are not found (off the grid if the grid is too small), they will not
	be included in the dictionary.
	'''
	#find cumulative density that is used to pass the given
	#sigma thresholds in "contours"
	cum_thresholds = 2*(norm.cdf(sigmas)-norm.cdf(0))

	integrate = 0.0
	curcontour = 0
	densitycontours = {}
	
	#sort and reverse Z
	Z = numpy.sort(Z,kind='quicksort',axis=None)[::-1]
	#augment Z by grid spacing, assuming density units are m^-2, but
	#spacing is not 1 cell m^-2
	Z *= gridspacings[0]*gridspacings[1]
	
	for d in Z:
		integrate+=d
		#if the elapsed density surpasses the next contour
		if (integrate >= cum_thresholds[curcontour]):
			densitycontours[sigmas[curcontour]] = d
			
			curcontour += 1
			if (curcontour>=len(sigmas)):
				break
	
	return densitycontours

def densityToRaster(griddata,ranges,spacings,outrastername,clon=-999,utmzone=-999,planet='earth',driver='GTiff',outproj="tm"):
	'''
	Outputs a 2-D array to a gdal-readable raster. Input expected to be
	in a transverse mercator projection.
	griddata: 2D data array
	outrastername: file name of raster output. If meter output is desired,
	   it would be good practice to define clon or utm zone
	planet: 'earth','venus', or 'mars'. This is only needed to translate to
	   latlong projections
	clon: center longitude of non-earth transverse mercator data
	utmzone: utm zone of earth data
	driver: gdal-readable raster short name [GTiff]
	outproj: 'tm' or 'll' for transverse mercator (no tranformation occurs)
	   or latlong (gdalwarp will be implemented). [tm]
	'''
	
	gdaldriver = gdal.GetDriverByName(driver)
	driver.Register()
	cols = numpy.shape(griddata)[0]
	rows = numpy.shape(griddata)[1]
	bands = 1
	
	
	dest_raster = gdaldriver.Create(outrastername, cols, rows, bands, gdal.GDT_Float64 )

	#adfGeoTransform[0] /* top left x */
  #adfGeoTransform[1] /* w-e pixel resolution */
  #adfGeoTransform[2] /* rotation, 0 if image is "north up" */
  #adfGeoTransform[3] /* top left y */
  #adfGeoTransform[4] /* rotation, 0 if image is "north up" */
  #adfGeoTransform[5] /* n-s pixel resolution */	
	geotrans = [ranges[0][0],spacings[0],0,ranges[1][1],0,(-1*spacings[1])]
	dest_raster.SetGeoTransform(geotrans)
	
	#set transverse mercator projection
	if (utmzone>=1 and utmzone<=60):
		srs = osr.SpatialReference()
		srs.SetUTM( utmzone, 1 ) #1 means north, this could be problematic
		srs.SetWellKnownGeogCS( 'WGS84' );
		dest_raster.SetProjection( srs.ExportToWkt() )

	elif (clon>=-360 and clon<=360):
		srs = osr.SpatialReference()
		
		if (planet == 'venus'):
			srs.ImportFromProj4( '+proj=tmerc +lat_0=0 +lon_0='+str(clon)+' +k=0.9996 +x_0=0 +y_0=0 +a=6051800 +b=6051800 +units=m +no_defs' )
		
		elif (planet == 'mars'):
			srs.ImportFromProj4( '+proj=tmerc +lat_0=0 +lon_0='+str(clon)+' +k=0.9996 +x_0=0 +y_0=0 +a=3396190 +b=3396190 +units=m +no_defs' )
		else:
			print 'error: clon set but planet is not venus or mars.'
			return 0
		dest_raster.SetProjection( srs.ExportToWkt() )
	
	dest_raster.GetRasterBand(1).WriteArray( griddata )
	
	#warp to ll if necessary
	if outproj=='ll':
		if planet=='earth':
			#catch invalid utmzone
			if ((utmzone<1) or (utmzone>60)):
				print 'utm zone not valid (1-60). Cannot create latlong raster.'
				print 'utm raster saved at '+outrastername
				return 0
			
			#reproject the transverse mercator grid
			os.system('gdalwarp -t_srs "+proj=longlat +datum=WGS84" '+outrastername+' tmpLL.tif')
			
		#if not earth, catch invalid center_lon
		elif ((clon<-360) or (clon>360)):
			print 'center longitude not valid (-360 to 360). Cannot create latlong raster.'
			print 'transverse mercator raster saved at '+outrastername
			return 0
		else:
			if planet=='mars':
				radius='3396190'
			elif planet=='venus':
				radius='6051800'
			else:
				print 'planet not either earth, venus, or mars. cannot create latlong raster.'
				print 'transverse mercator raster saved at '+outrastername
				return 0
			
			#reproject the transverse mercator grid
			os.system('gdalwarp -t_srs "+proj=longlat +k=0.9996 +x_0=0 +y_0=0 +a='+radius+' +b='+radius+' +no_defs" '+outrastername+' tmpLL.tif')
			
		#overwrite the transverse meter raster with the longlat raster file
		if driver=='GTiff':
			os.system('mv tmpLL.tif '+outrastername)
		else:
			os.system('gdal_translate -of '+driver+' tmpLL.tif '+outrastername)
		os.remove('tmpLL.tif')

