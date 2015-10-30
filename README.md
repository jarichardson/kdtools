# KDtools
## Kernel Density Estimation Tools

## DESCRIPTION
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

## AUTHOR
Jacob Richardson (github.com/jarichardson)

## REQUIREMENTS
* Python Libraries: 
	* pyproj
	* numpy
	* os
	* time
	* osgeo
	* matplotlib - For main (test) function
* PROJ.4	- For reproject, densityToRaster functions
* GDAL	- For densityToRaster
* R	- For samse_bandwidth function
	* Also requires the KS library in R. See [Installing KS Library Section](#installing-ks-library-in-r) below.
* GMT	-For ellipseGen function

## FUNCTIONS
**contourBySigma**(Z,sigmas=[0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0],gridspacings=[1,1])
	
	Identifies Density contour levels given Sigma thresholds.
	Contours define areas within which lay X% of the total field density
	e.g.: Within the Sigma-2.0 contour lies 95.4% of total field density.
	      The density value of each contour decreases with increasing
	      sigma.
	
	Requires a numpy array of density values (*Z*), with any shape.
	Optionally, provide a list of requested sigma thresholds, and the
	grid size as a 2 item list, to normalize the density.
	
	Outputs a dictionary of {sigma-level: density value}. If sigma-levels
	are not found (off the grid if the grid is too small), they will not
	be included in the dictionary.
		
**densityToRaster**(griddata,ranges,spacings,outrastername,clon=-999,utmzone=-999,planet='earth',driver='GTiff',outproj="tm")
	
	Outputs a 2-D array to a gdal-readable raster. Input expected to be
	in a transverse mercator projection.
	
	griddata: 2D data array
	outrastername - file name of raster output. If meter output is desired,
	                it would be good practice to define clon or utm zone
	planet        - 'earth','venus', or 'mars'. This is only needed to 
	                translate to latlong projections
	clon          - center longitude of non-earth transverse mercator data
	utmzone       - utm zone of earth data
	driver        - gdal-readable raster short name [GTiff]
	outproj       - 'tm' or 'll' for transverse mercator (no tranformation
	                occurs) or latlong (gdalwarp will be implemented). [tm]
	   
	ISSUES: If values are very low (normal for density grids), gdalwarp 
	     doesn't work, so it is suggested that output remain in meters.
	     A workaround might be to supply log10 values of griddata.
		
**ellipseGen**(bd,eps=False,epsfilename='bandwidth_ellipse.eps')
	
	Identifies the major and minor axes directions and standard
	deviations of a Gaussian ellipse defined by a 2x2 covariance
	matrix. Precision is to the nearest degree.
	
	Prints out solution, and optionally uses GMT to draw the ellipse
	to epsfilename, if *eps*=True.
	
	Outputs major-axis direction, major-axis standard deviation, and
	        minor-axis standard-deviation.
	
**KD**(bd,coors,ranges,spacings):
	
	Estimates point density using:
	bd       - a kernel bandwidth (2x2 covariance	matrix)
	coors    - 2xN list of coordinates for N points.
	ranges   - a 2x2 [[W,E],[S,N]] array
	spacings - a 1x2 [X-resolution,Y-resolution] array
	weights  - a 2xN list of weights for N points (default: empty [])
	
	Outputs X,Y,D: Eastings, Northings, and Densities in a Meshgrid
	format (i.e. X will be tiled, Y will be repeated)

**main**()

	Runs tests for kdtools functions using a synthetic dataset.
	Some tests are visual and require matplotlib.

**rangeBuffer**(coords,B=0)
	
	Creates a buffer of *B*% [default 0%, no buffer] around 
	N-dimensional data. Input should be a numpy array.
	Output will be 2xN array, with min, max of each dimension in
	columns 1 and 2, respectively.
	
	ex: 
	data         range output
	[[1,5],      [[0,2],
	 [2,5],  =>   [4,9]]
	 [0,4],
	 [1,9]]
	
**reproject**(llcoors,planet="earth",utmzone=-999,clon=-999,inverse=False)
	
	Reprojects long-lat data into transverse mercator coordinates
	with units of meters. Optionally, set *inverse*=True to change
	transverse mercator to long-lat.
	
	Input should be numpy 2-col long, lat array.
	Output will be numpy 2-col easting, northing array.	
	
	*planet* options: 'earth', 'venus', or 'mars'
	Earth requires a valid UTM zone
	Venus and Mars require a valid center longitude of the dataset
	
	Earth Transverse Mercator fit to WGS84 datum
	Venus Transverse Mercator fit to Spheriod of radius 6051800 m
	Mars Transverse Mercator fit to Spheriod of radius 3396190 m
	
**samse_bandwidth**(coords)
	
	Evaluates the SAMSE Kernel in R using a coordinate list (coords).
	Returns 2x2 bandwidth covariance matrix.
	Requires: R, KS library in R.
	
## INSTALLING KS LIBRARY IN R

To install ks in R on a unix machine, you will need the following programs installed:

* R-base
* R-base-devel
* Mesa
* Mesa-devel
* gcc
* gcc-fortran
* libpng14-devel
* libpng12-devel
* libpng12-compat-devel
* possibly others depending on the errors given while installing packages in R.

Then in R (in superuser), run:

	install.packages("ks")

and choose a mirror. The KS package requires other packages:
* rgl
* mvtnorm

If the above utilities are installed, these two packages should install easily.
If the two packages are installed, KS should install easily. 
To check KS installation, run (in R, out of superuser)

	library(ks)
