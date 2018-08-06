import hsc.pipe.tasks.processCcd
assert type(root)==hsc.pipe.tasks.processCcd.SubaruProcessCcdConfig, 'config is of type %s.%s instead of hsc.pipe.tasks.processCcd.SubaruProcessCcdConfig' % (type(root).__module__, type(root).__name__)
import lsst.meas.extensions.photometryKron.version
import hsc.meas.astrom.astrom
import lsst.meas.multifit.measureMulti
import lsst.meas.extensions.multiShapelet.multiShapeletLib
import lsst.obs.subaru.astrometry
import lsst.meas.mosaic
import hsc.meas.astrom
import lsst.meas.extensions.photometryKron.kronLib
import lsst.pipe.tasks.coaddInputRecorder
import lsst.meas.extensions.shapeHSM.version
import lsst.meas.multifit
import lsst.pipe.tasks.selectImages
import lsst.meas.extensions.shapeHSM.hsmLib
import lsst.shapelet.shapeletLib
import lsst.pipe.tasks.coaddBase
import hsc.meas.astrom.astromLib
import lsst.meas.multifit.measureCcd
import hsc.meas
import lsst.obs.subaru.crosstalkYagi
import lsst.shapelet
import lsst.meas.multifit.priors
import lsst.meas.extensions.psfex
import lsst.meas.extensions.psfex.psfexLib
import lsst.meas.mosaic.updateExposure
import lsst.obs.hsc.vignette
import lsst.shapelet.tractor
import lsst.meas.multifit.baseMeasure
import lsst.obs.subaru.subaruLib
import lsst.meas.multifit.measureCoadd
import lsst.obs.subaru.isr
import lsst.meas.multifit.models
import lsst.meas.multifit.version
import lsst.meas.extensions.photometryKron
import lsst.meas.extensions.psfex.psfexPsfDeterminer
import lsst.meas.extensions
import lsst.meas.multifit.fitRegion
import lsst.afw.display.rgb
import lsst.meas.extensions.multiShapelet.version
import lsst.obs.subaru
import lsst.obs.subaru.crosstalk
import lsst.meas.extensions.shapeHSM
import lsst.meas.multifit.multifitLib
import lsst.meas.multifit.measureImage
import lsst.meas.multifit.samplers
import lsst.meas.multifit.optimizer
import lsst.meas.mosaic.mosaicLib
import lsst.meas.extensions.multiShapelet
'''Perform ISR? '''
root.doIsr=True

'''Write the denormalized match table as well as the normalized match table (in the final write)? '''
root.doWriteUnpackedMatches=True

'''Include HeavyFootprint data in source table? '''
root.doWriteHeavyFootprintsInSources=False

'''Write all outputs at end?  If so, you also likely want doWriteCalibrate=False. '''
root.doFinalWrite=True

'''Write sources? '''
root.doWriteSources=True

import lsst.obs.subaru.isr
root.isr.retarget(target=lsst.obs.subaru.isr.SubaruIsrTask, ConfigClass=lsst.obs.subaru.isr.SubaruIsrConfig)
'''Number of stdev below the background to set thumbnail minimum '''
root.isr.thumbnailStdev=3.0

'''Assemble detrend/calibration frames? '''
root.isr.doAssembleDetrends=False

'''How to estimate the average value for BAD regions.
Allowed values:
	None	Field is optional
	MEDIAN	Correct using the median of the good data
	MEANCLIP	Correct using the (clipped) mean of good data
 '''
root.isr.badStatistic='MEANCLIP'

'''Widen bleed trails based on their width? '''
root.isr.doWidenSaturationTrails=True

import lsst.obs.subaru.crosstalk
root.isr.crosstalk.retarget(target=lsst.obs.subaru.crosstalk.CrosstalkTask, ConfigClass=lsst.obs.subaru.crosstalk.CrosstalkConfig)
'''Shape of coeffs array '''
root.isr.crosstalk.coeffs.shape=[4, 4]

'''Crosstalk coefficients '''
root.isr.crosstalk.coeffs.values=[0.0, -0.000125, -0.000149, -0.000156, -0.000124, 0.0, -0.000132, -0.000157, -0.000171, -0.000134, 0.0, -0.000153, -0.000157, -0.000151, -0.000137, 0.0]

'''Name for crosstalk mask plane '''
root.isr.crosstalk.crosstalkMaskPlane='CROSSTALK'

'''Set crosstalk mask plane for pixels over this value '''
root.isr.crosstalk.minPixelToMask=45000.0

'''Order of polynomial or to fit if overscan fit type is a polynomial, or number of spline knots if overscan fit type is a spline. '''
root.isr.overscanPolyOrder=30

'''Apply dark frame correction? '''
root.isr.doDark=True

'''update exposure metadata in the assembled ccd to reflect the effective gain of the assembled chip '''
root.isr.setGainAssembledCcd=True

'''Fallback default filter name for calibrations '''
root.isr.fallbackFilterName=None

'''Correct for crosstalk '''
root.isr.doCrosstalk=True

'''FWHM of PSF used when interpolating over bad columns (arcsec) '''
root.isr.fwhmForBadColumnInterpolation=1.0

'''Number of points to define the Vignette polygon '''
root.isr.numPolygonPoints=100

'''FWHM of PSF (arcsec) '''
root.isr.fwhm=1.0

'''Maximum number of iterations for the brighter fatter correction '''
root.isr.brighterFatterMaxIter=10

'''If flatScalingType is 'USER' then scale flat by this amount; ignored otherwise '''
root.isr.flatUserScale=1.0

'''Rejection threshold (sigma) for collapsing overscan before fit '''
root.isr.overscanRej=3.0

'''Name of mask plane to use in saturation detection and interpolation '''
root.isr.saturatedMaskName='SAT'

'''Should the gain be applied when applying the brighter fatter correction? '''
root.isr.brighterFatterApplyGain=True

'''Should we set the level of all BAD patches of the chip to the chip's average value? '''
root.isr.doSetBadRegions=True

'''Correct for nonlinearity of the detector's response (ignored if coefficients are 0.0) '''
root.isr.doLinearize=True

'''Kernel file used for the brighter fatter correction '''
root.isr.brighterFatterKernelFile='/projects/p025_swin/pipes/hscpipe/4.0.5/Linux64/obs_subaru/HSC-4.0.3/hsc/brighter_fatter_kernel.pkl'

'''Apply the brighter fatter correction '''
root.isr.doBrighterFatter=True

'''Apply bias frame correction? '''
root.isr.doBias=True

'''Apply flat field correction? '''
root.isr.doFlat=True

'''Remove any PC cards in the header '''
root.isr.removePcCards=True

'''Apply fringe correction? '''
root.isr.doFringe=True

'''trim out non-data regions? '''
root.isr.assembleCcd.doTrim=True

'''FITS headers to remove (in addition to DATASEC, BIASSEC, TRIMSEC and perhaps GAIN) '''
root.isr.assembleCcd.keysToRemove=[]

'''renormalize to a gain of 1? (ignored if setGain false) '''
root.isr.assembleCcd.doRenorm=False

'''set gain? '''
root.isr.assembleCcd.setGain=True

'''Calculate variance? '''
root.isr.doVariance=True

'''Default value for fluxMag0T1 (for an unrecognised filter) '''
root.isr.defaultFluxMag0T1=28.0

'''Softening parameter for thumbnail mapping '''
root.isr.thumbnailQ=20.0

'''fields to remove from the metadata of the assembled ccd. '''
root.isr.keysToRemoveFromAssembledCcd=[]

'''Center of vignetting pattern, in x (focal plane coords) '''
root.isr.vignette.xCenter=-100.0

'''Radius of vignetting pattern, in focal plane coords '''
root.isr.vignette.radius=17500.0

'''Center of vignetting pattern, in y (focal plane coords) '''
root.isr.vignette.yCenter=100.0

'''The approximate flux of a zero-magnitude object in a one-second exposure, per filter '''
root.isr.fluxMag0T1={'g': 398107170553.49854, 'N816': 15848931924.611174, 'i': 275422870333.81744, 'r': 398107170553.49854, 'N921': 19054607179.632523, 'N515': 20892961308.54041, 'y': 91201083935.59116, 'z': 120226443461.74132}

'''Do overscan subtraction? '''
root.isr.doOverscan=True

'''Binning factor for thumbnail '''
root.isr.thumbnailBinning=4

'''Do fringe subtraction after flat-fielding? '''
root.isr.fringeAfterFlat=True

'''Border around saturated pixels for thumbnail '''
root.isr.thumbnailSatBorder=2

'''Mask saturated pixels? '''
root.isr.doSaturation=True

'''Trim guider shadow '''
root.isr.doGuider=False

'''The method for scaling the flat on the fly.
Allowed values:
	None	Field is optional
	MEDIAN	Scale by the inverse of the median
	USER	Scale by flatUserScale
	MEAN	Scale by the inverse of the mean
 '''
root.isr.flatScalingType='USER'

'''Number of pixels by which to grow the saturation footprints '''
root.isr.growSaturationFootprintSize=1

'''Correct the amplifiers for their gains

N.b. this is intended to be used *instead* of doFlat; it's useful if you're measuring system throughput
 '''
root.isr.doApplyGains=False

'''Persist Polygon used to define vignetted region? '''
root.isr.doWriteVignettePolygon=True

'''Offset to the random number generator seed (full seed includes exposure ID) '''
root.isr.fringe.stats.rngSeedOffset=0

'''Ignore pixels with these masks '''
root.isr.fringe.stats.badMaskPlanes=['SAT', 'NO_DATA']

'''Statistic to use '''
root.isr.fringe.stats.stat=32

'''Number of fitting iterations '''
root.isr.fringe.stats.iterations=3

'''Sigma clip threshold '''
root.isr.fringe.stats.clip=3.0

'''Only fringe-subtract these filters '''
root.isr.fringe.filters=['y', 'N921']

'''Sigma clip threshold '''
root.isr.fringe.clip=3.0

'''Half-size of large (background) measurements (pixels) '''
root.isr.fringe.large=30

'''Number of fringe measurements '''
root.isr.fringe.num=30000

'''Number of fitting iterations '''
root.isr.fringe.iterations=20

'''Half-size of small (fringe) measurements (pixels) '''
root.isr.fringe.small=3

'''Remove fringe pedestal? '''
root.isr.fringe.pedestal=False

'''Threshold used to stop iterating the brighter fatter correction.  It is the  absolute value of the difference between the current corrected image and the one from the previous iteration summed over all the pixels. '''
root.isr.brighterFatterThreshold=1000.0

'''Write OverScan-Subtracted thumbnail? '''
root.isr.qa.doThumbnailOss=False

'''Mesh size in X (pix) to calculate count statistics '''
root.isr.qa.flatness.meshX=256

'''How many times do we iterate clipping outliers in calculate count statistics? '''
root.isr.qa.flatness.nIter=3

'''Do we clip outliers in calculate count statistics? '''
root.isr.qa.flatness.doClip=True

'''How many sigma is used to clip outliers in calculate count statistics? '''
root.isr.qa.flatness.clipSigma=3.0

'''Mesh size in Y (pix) to calculate count statistics '''
root.isr.qa.flatness.meshY=256

'''Write flattened thumbnail? '''
root.isr.qa.doThumbnailFlattened=False

'''Write OverScan-Subtracted image? '''
root.isr.qa.doWriteOss=False

'''Write flattened image? '''
root.isr.qa.doWriteFlattened=False

'''Range for thumbnail mapping '''
root.isr.thumbnailRange=5.0

'''Persist postISRCCD? '''
root.isr.doWrite=False

'''Mask defect pixels? '''
root.isr.doDefect=True

'''Normalize all the amplifiers in each CCD to have the same gain

This does not measure the gains, it simply forces the median of each amplifier to be equal
after applying the nominal gain
 '''
root.isr.normalizeGains=False

'''Maximum deviation from the median for overscan '''
root.isr.overscanMaxDev=1000.0

'''The method for fitting the overscan bias level.
Allowed values:
	None	Field is optional
	LEG	Fit Legendre polynomial to the longest axis of the overscan region
	CUBIC_SPLINE	Fit cubic spline to the longest axis of the overscan region
	MEDIAN	Correct using the median of the overscan region
	POLY	Fit ordinary polynomial to the longest axis of the overscan region
	CHEB	Fit Chebyshev polynomial to the longest axis of the overscan region
	AKIMA_SPLINE	Fit Akima spline to the longest axis of the overscan region
	NATURAL_SPLINE	Fit natural spline to the longest axis of the overscan region
	MEAN	Correct using the mean of the overscan region
 '''
root.isr.overscanFitType='AKIMA_SPLINE'

'''Tweak flats to match observed amplifier ratios? '''
root.isr.doTweakFlat=False

'''Detect sources? '''
root.doDetection=False

'''Make plots? '''
root.qa.seeing.doPlots=False

'''Size of grid (pixels) '''
root.qa.seeing.gridSize=1024.0

'''Type of psf for contour plot (psfcand:psf candidates, model:psf model or both:both)
Allowed values:
	both	will plot both of the psf candidates and psf model separately
	model	Psf profiles directly taken from determined Psf model
	None	Field is optional
	psfcand	Psf profiles generated by stacking psf candidates in each grid
 '''
root.qa.seeing.psfContourType='both'

'''Undistort when evaluating the 2nd moments of sources? '''
root.qa.seeing.starsel.doUndistort=False

'''How many sigmas around the peak fwhm are used for calculating statistics of PSF sequence '''
root.qa.seeing.starsel.fwhmMarginNsigma=3.0

'''Size of grid (pixels) '''
root.qa.seeing.starsel.gridSize=1024.0

'''Minimum fwhm allowed in estimation of seeing (pix) '''
root.qa.seeing.starsel.fwhmMin=1.5

'''Number of smallest objects which are used to determine rough-interim seeing '''
root.qa.seeing.starsel.nSmallSampleRoughFwhm=50

'''How many sigmas around the peak fwhm are used for calculating statistics of PSF sequence '''
root.qa.seeing.starsel.psfSeqStatNsigma=3.0

'''Number of bins for number counting as a fn of instrumnetal mag '''
root.qa.seeing.starsel.nbinMagHist=80

'''What fraction of sources from the brightest is to be included for initial guess of seeing to avoid cosmic rays which dominate faint magnitudes '''
root.qa.seeing.starsel.fracSrcIni=0.15

'''How many times do we iterate calculating statistics of PSF sequence '''
root.qa.seeing.starsel.psfSeqStatNiter=3

'''size of the kernel to create '''
root.qa.seeing.starsel.kernelSize=21

'''Bin size of FWHM histogram '''
root.qa.seeing.starsel.fwhmBinSize=0.2

'''Faintest mag for number counting as a fn of instrumnetal mag '''
root.qa.seeing.starsel.magMaxHist=0.0

'''number of pixels to ignore around the edge of PSF candidate postage stamps '''
root.qa.seeing.starsel.borderWidth=0

'''Number of brightest (non-saturated) objects which are used to determine rough-interim seeing '''
root.qa.seeing.starsel.nBrightSampleRoughFwhm=30

'''Brightest mag for number counting as a fn of instrumnetal mag '''
root.qa.seeing.starsel.magMinHist=-20.0

'''How many pixels around the peak are used for calculating scatter of psf candidates '''
root.qa.seeing.starsel.fwhmMarginFinal=0.75

'''Make plots? '''
root.qa.seeing.starsel.doPlots=False

'''Statistical algorithm to derive rough Fwhm in the 1st step seeing estimation
Allowed values:
	None	Field is optional
	MEDIAN	median of sample
	MEANCLIP	clipped mean of sample with 3-sigma clip + 3-times iteration
 '''
root.qa.seeing.starsel.statAlgRoughFwhm='MEDIAN'

'''How many magnitudes to extend the faint-end limit for extracting PSF sources, from the base magnitude determined by fracSrcIni. '''
root.qa.seeing.starsel.magLimitFaintExtension=0.0

'''Maxmum fwhm allowed in estimation of seeing (pix) '''
root.qa.seeing.starsel.fwhmMax=12.0

'''Use icsources(calib.sources) rather than final sources '''
root.qa.useIcsources=False

'''Compute and write src to reference matches? '''
root.doWriteSourceMatches=True

'''Estimate the background again after final source detection? '''
root.detection.reEstimateBackground=True

'''detected sources with fewer than the specified number of pixels will be ignored
	Valid Range = [0,inf) '''
root.detection.minPixels=1

'''Grow detections by nSigmaToGrow * sigma; if 0 then do not grow '''
root.detection.nSigmaToGrow=2.4

'''Names of mask planes to ignore while estimating the background '''
root.detection.footprintBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

'''behaviour if there are too few points in grid for requested interpolation style
Allowed values:
	None	Field is optional
	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
	THROW_EXCEPTION	throw an exception if there are too few points
	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
 '''
root.detection.footprintBackground.undersampleStyle='REDUCE_INTERP_ORDER'

'''how to interpolate the background values. This maps to an enum; see afw::math::Background
Allowed values:
	None	Field is optional
	CONSTANT	Use a single constant value
	LINEAR	Use linear interpolation
	NONE	No background estimation is to be attempted
	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
 '''
root.detection.footprintBackground.algorithm='AKIMA_SPLINE'

'''how large a region of the sky should be used for each background point
	Valid Range = [10,inf) '''
root.detection.footprintBackground.binSize=64

'''Ignore NaNs when estimating the background '''
root.detection.footprintBackground.isNanSafe=False

'''type of statistic to use for grid points
Allowed values:
	None	Field is optional
	MEDIAN	median
	MEANCLIP	clipped mean
	MEAN	unclipped mean
 '''
root.detection.footprintBackground.statisticsProperty='MEANCLIP'

'''Apprimation order for background Chebyshev (valid only with useApprox=True) '''
root.detection.footprintBackground.approxOrder=6

'''Use Approximate (Chebyshev) to model background. '''
root.detection.footprintBackground.useApprox=False

'''Do background subtraction before footprint detection? '''
root.detection.doFootprintBackground=False

'''Include threshold relative to thresholdValue
	Valid Range = [0.0,inf) '''
root.detection.includeThresholdMultiplier=1.0

'''Pixels should be grown as isotropically as possible (slower) '''
root.detection.isotropicGrow=True

'''Fiddle factor to add to the background; debugging only '''
root.detection.adjustBackground=0.0

'''specifies the desired flavor of Threshold
Allowed values:
	pixel_stdev	threshold applied to per-pixel std deviation
	variance	threshold applied to image variance
	value	threshold applied to image value
	stdev	threshold applied to image std deviation
 '''
root.detection.thresholdType='stdev'

'''Names of mask planes to ignore while estimating the background '''
root.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

'''behaviour if there are too few points in grid for requested interpolation style
Allowed values:
	None	Field is optional
	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
	THROW_EXCEPTION	throw an exception if there are too few points
	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
 '''
root.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

'''how to interpolate the background values. This maps to an enum; see afw::math::Background
Allowed values:
	None	Field is optional
	CONSTANT	Use a single constant value
	LINEAR	Use linear interpolation
	NONE	No background estimation is to be attempted
	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
 '''
root.detection.background.algorithm='NATURAL_SPLINE'

'''how large a region of the sky should be used for each background point
	Valid Range = [10,inf) '''
root.detection.background.binSize=128

'''Ignore NaNs when estimating the background '''
root.detection.background.isNanSafe=False

'''type of statistic to use for grid points
Allowed values:
	None	Field is optional
	MEDIAN	median
	MEANCLIP	clipped mean
	MEAN	unclipped mean
 '''
root.detection.background.statisticsProperty='MEANCLIP'

'''Apprimation order for background Chebyshev (valid only with useApprox=True) '''
root.detection.background.approxOrder=6

'''Use Approximate (Chebyshev) to model background. '''
root.detection.background.useApprox=True

'''Grow detections to set the image mask bits, but return the original (not-grown) footprints '''
root.detection.returnOriginalFootprints=False

'''specifies whether to detect positive, or negative sources, or both
Allowed values:
	positive	detect only positive sources
	negative	detect only negative sources
	both	detect both positive and negative sources
 '''
root.detection.thresholdPolarity='positive'

'''Threshold for footprints
	Valid Range = [0.0,inf) '''
root.detection.thresholdValue=5.0

'''Write icSrc to reference matches? '''
root.doWriteCalibrateMatches=True

'''Deblend sources? '''
root.doDeblend=False

'''PSF model type
Allowed values:
	None	Field is optional
	DoubleGaussian	Double Gaussian model
	SingleGaussian	Single Gaussian model
 '''
root.calibrate.initialPsf.model='SingleGaussian'

'''FWHM of PSF model (arcsec) '''
root.calibrate.initialPsf.fwhm=1.0

'''Size of PSF model (pixels) '''
root.calibrate.initialPsf.size=15

'''Names of mask planes to ignore while estimating the background '''
root.calibrate.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

'''behaviour if there are too few points in grid for requested interpolation style
Allowed values:
	None	Field is optional
	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
	THROW_EXCEPTION	throw an exception if there are too few points
	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
 '''
root.calibrate.background.undersampleStyle='REDUCE_INTERP_ORDER'

'''how to interpolate the background values. This maps to an enum; see afw::math::Background
Allowed values:
	None	Field is optional
	CONSTANT	Use a single constant value
	LINEAR	Use linear interpolation
	NONE	No background estimation is to be attempted
	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
 '''
root.calibrate.background.algorithm='NATURAL_SPLINE'

'''how large a region of the sky should be used for each background point
	Valid Range = [10,inf) '''
root.calibrate.background.binSize=128

'''Ignore NaNs when estimating the background '''
root.calibrate.background.isNanSafe=False

'''type of statistic to use for grid points
Allowed values:
	None	Field is optional
	MEDIAN	median
	MEANCLIP	clipped mean
	MEAN	unclipped mean
 '''
root.calibrate.background.statisticsProperty='MEANCLIP'

'''Apprimation order for background Chebyshev (valid only with useApprox=True) '''
root.calibrate.background.approxOrder=6

'''Use Approximate (Chebyshev) to model background. '''
root.calibrate.background.useApprox=True

'''Interpolate over defects? (ignored unless you provide a list of defects) '''
root.calibrate.repair.doInterpolate=True

'''Smoothly taper (on the PSF scale) to the fallback value at the edge of the image? '''
root.calibrate.repair.interp.useFallbackValueAtEdge=True

'''Interpolation kernel size = interpFwhm converted to pixels * interpKernelSizeFactor. '''
root.calibrate.repair.interp.interpKernelSizeFactor=3.0

'''Find and mask out cosmic rays? '''
root.calibrate.repair.doCosmicRay=True

'''Don't interpolate over CR pixels '''
root.calibrate.repair.cosmicray.keepCRs=False

'''used in condition 3 for CR; see CR.cc code '''
root.calibrate.repair.cosmicray.cond3_fac=2.5

'''used in condition 3 for CR; see CR.cc code '''
root.calibrate.repair.cosmicray.cond3_fac2=0.4

'''Names of mask planes to ignore while estimating the background '''
root.calibrate.repair.cosmicray.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

'''behaviour if there are too few points in grid for requested interpolation style
Allowed values:
	None	Field is optional
	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
	THROW_EXCEPTION	throw an exception if there are too few points
	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
 '''
root.calibrate.repair.cosmicray.background.undersampleStyle='REDUCE_INTERP_ORDER'

'''how to interpolate the background values. This maps to an enum; see afw::math::Background
Allowed values:
	None	Field is optional
	CONSTANT	Use a single constant value
	LINEAR	Use linear interpolation
	NONE	No background estimation is to be attempted
	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
 '''
root.calibrate.repair.cosmicray.background.algorithm='AKIMA_SPLINE'

'''how large a region of the sky should be used for each background point
	Valid Range = [10,inf) '''
root.calibrate.repair.cosmicray.background.binSize=100000

'''Ignore NaNs when estimating the background '''
root.calibrate.repair.cosmicray.background.isNanSafe=False

'''type of statistic to use for grid points
Allowed values:
	None	Field is optional
	MEDIAN	median
	MEANCLIP	clipped mean
	MEAN	unclipped mean
 '''
root.calibrate.repair.cosmicray.background.statisticsProperty='MEDIAN'

'''Apprimation order for background Chebyshev (valid only with useApprox=True) '''
root.calibrate.repair.cosmicray.background.approxOrder=6

'''Use Approximate (Chebyshev) to model background. '''
root.calibrate.repair.cosmicray.background.useApprox=False

'''number of times to look for contaminated pixels near known CR pixels '''
root.calibrate.repair.cosmicray.niteration=3

'''maximum number of contaminated pixels '''
root.calibrate.repair.cosmicray.nCrPixelMax=1000000

'''CRs must be > this many sky-sig above sky '''
root.calibrate.repair.cosmicray.minSigma=6.0

'''CRs must have > this many DN (== electrons/gain) in initial detection '''
root.calibrate.repair.cosmicray.min_DN=150.0

'''Minimum number of degrees of freedom (# of valid data points - # of parameters); if this is exceeded, the order of the fit is decreased (in both dimensions), and if we can't decrease it enough, we'll raise ValueError. '''
root.calibrate.measureApCorr.minDegreesOfFreedom=1

'''if true, only include terms where the sum of the x and y order is less than or equal to max(orderX, orderY) '''
root.calibrate.measureApCorr.fit.triangular=True

'''maximum Chebyshev function order in x '''
root.calibrate.measureApCorr.fit.orderX=2

'''maximum Chebyshev function order in y '''
root.calibrate.measureApCorr.fit.orderY=2

'''Name of a flag field that indicates that a source should be used to constrain the aperture corrections '''
root.calibrate.measureApCorr.inputFilterFlag='calib.psf.used'

'''Number of standard devisations to clip at '''
root.calibrate.measureApCorr.numSigmaClip=3.0

'''Number of iterations for sigma clipping '''
root.calibrate.measureApCorr.numIter=4

'''Name of the flux field other measurements should be corrected to match '''
root.calibrate.measureApCorr.reference='flux.naive'

'''Perform PSF fitting? '''
root.calibrate.doPsf=False

import lsst.obs.subaru.astrometry
root.calibrate.astrometry.retarget(target=lsst.obs.subaru.astrometry.SubaruAstrometryTask, ConfigClass=lsst.obs.subaru.astrometry.SubaruAstrometryConfig)
'''Rejection threshold for Wcs fitting
	Valid Range = (0.0,inf) '''
root.calibrate.astrometry.rejectThresh=3.0

'''number of points to define a shape for matching '''
root.calibrate.astrometry.solver.numPointsForShape=6

'''Use the parity (flip / handedness) of the image from the input exposure's WCS headers? '''
root.calibrate.astrometry.solver.useWcsParity=True

'''Minimum number of matched pairs
	Valid Range = [2,inf) '''
root.calibrate.astrometry.solver.minMatchedPairNumber=30

'''Maximum CPU time to spend solving, in seconds
	Valid Range = [0.0,inf) '''
root.calibrate.astrometry.solver.maxCpuTime=0.0

'''Difference of angle between x and y from 90 degree allowed (degree)
	Valid Range = [-inf,45.0) '''
root.calibrate.astrometry.solver.angleDiffFrom90=0.2

'''Offset between sources and catalog allowed (pixel)
	Valid Range = [-inf,4000) '''
root.calibrate.astrometry.solver.offsetAllowedInPixel=300

'''Polynomial order of SIP distortion terms
	Valid Range = [1,inf) '''
root.calibrate.astrometry.solver.sipOrder=3

'''Matching threshold for Astrometry.net solver (log-odds)
	Valid Range = [13.815510558,inf) '''
root.calibrate.astrometry.solver.matchThreshold=27.631021115928547

'''Use the pixel scale from the input exposure's WCS headers? '''
root.calibrate.astrometry.solver.useWcsPixelScale=True

'''Padding to add to image size (pixels)
	Valid Range = [0,inf) '''
root.calibrate.astrometry.solver.pixelMargin=50

'''Roation angle allowed between sources and catalog (radian)
	Valid Range = [-inf,0.1) '''
root.calibrate.astrometry.solver.rotationAllowedInRad=0.02

'''Sigma-clipping parameter in sip/cleanBadPoints.py
	Valid Range = [0.0,inf) '''
root.calibrate.astrometry.solver.cleaningParameter=3.0

'''Number of bright stars to use
	Valid Range = [2,inf) '''
root.calibrate.astrometry.solver.numBrightStars=50

'''Use the RA,Dec center information from the input exposure's WCS headers? '''
root.calibrate.astrometry.solver.useWcsRaDecCenter=True

'''Maximum number of stars to use in Astrometry.net solving
	Valid Range = [10,inf) '''
root.calibrate.astrometry.solver.maxStars=50

'''When matching image to reference catalog stars, how big should
        the matching radius be?
	Valid Range = [0.0,inf) '''
root.calibrate.astrometry.solver.catalogMatchDist=2.0

'''limit on determinant of linear transforming matrix '''
root.calibrate.astrometry.solver.limitOnDeterminant=0.02

'''Retrieve all available fluxes (and errors) from catalog? '''
root.calibrate.astrometry.solver.allFluxes=True

'''Minimum number of matched pairs, expressed as a fraction of the reference catalogue size
	Valid Range = [0,1) '''
root.calibrate.astrometry.solver.minMatchedPairFrac=0.3

'''Range of pixel scales, around the value in the WCS header, to search.  If the value of this field is X and the nominal scale is S, the range searched will be  S/X to S*X
	Valid Range = [1.001,inf) '''
root.calibrate.astrometry.solver.pixelScaleUncertainty=1.1

'''When useWcsRaDecCenter=True, this is the radius, in degrees, around the RA,Dec center specified in the input exposure's WCS to search for a solution.
	Valid Range = [0.0,inf) '''
root.calibrate.astrometry.solver.raDecSearchRadius=1.0

'''Compute polynomial SIP distortion terms? '''
root.calibrate.astrometry.solver.calculateSip=True

'''Mapping from input filter to catalogue filter '''
root.calibrate.astrometry.solver.filterMap={'B': 'g', 'r2': 'r', 'N1010': 'z', 'N816': 'i', 'I': 'i', 'N387': 'g', 'i2': 'i', 'R': 'r', 'N921': 'z', 'N515': 'g', 'V': 'r'}

'''Fail over from hscAstrom to meas_astrom? '''
root.calibrate.astrometry.failover=False

'''Assume that the input image's WCS is correct, without comparing it to any external reality. (In contrast to using Astrometry.net).  NOTE, if you set this, you probably also want to un-set 'solver.calculateSip'; otherwise we'll still try to find a TAN-SIP WCS starting  from the existing WCS '''
root.calibrate.astrometry.forceKnownWcs=False

'''Rejection iterations for Wcs fitting
	Valid Range = [0,inf) '''
root.calibrate.astrometry.rejectIter=3

'''Proceed even if astrometry fails? '''
root.calibrate.astrometry.allowFailedAstrometry=False

'''Subtract background (after computing it, if not supplied)? '''
root.calibrate.doBackground=True

'''Estimate the background again after final source detection? '''
root.calibrate.detection.reEstimateBackground=True

'''detected sources with fewer than the specified number of pixels will be ignored
	Valid Range = [0,inf) '''
root.calibrate.detection.minPixels=1

'''Grow detections by nSigmaToGrow * sigma; if 0 then do not grow '''
root.calibrate.detection.nSigmaToGrow=2.4

'''Names of mask planes to ignore while estimating the background '''
root.calibrate.detection.footprintBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

'''behaviour if there are too few points in grid for requested interpolation style
Allowed values:
	None	Field is optional
	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
	THROW_EXCEPTION	throw an exception if there are too few points
	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
 '''
root.calibrate.detection.footprintBackground.undersampleStyle='REDUCE_INTERP_ORDER'

'''how to interpolate the background values. This maps to an enum; see afw::math::Background
Allowed values:
	None	Field is optional
	CONSTANT	Use a single constant value
	LINEAR	Use linear interpolation
	NONE	No background estimation is to be attempted
	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
 '''
root.calibrate.detection.footprintBackground.algorithm='AKIMA_SPLINE'

'''how large a region of the sky should be used for each background point
	Valid Range = [10,inf) '''
root.calibrate.detection.footprintBackground.binSize=64

'''Ignore NaNs when estimating the background '''
root.calibrate.detection.footprintBackground.isNanSafe=False

'''type of statistic to use for grid points
Allowed values:
	None	Field is optional
	MEDIAN	median
	MEANCLIP	clipped mean
	MEAN	unclipped mean
 '''
root.calibrate.detection.footprintBackground.statisticsProperty='MEANCLIP'

'''Apprimation order for background Chebyshev (valid only with useApprox=True) '''
root.calibrate.detection.footprintBackground.approxOrder=6

'''Use Approximate (Chebyshev) to model background. '''
root.calibrate.detection.footprintBackground.useApprox=False

'''Do background subtraction before footprint detection? '''
root.calibrate.detection.doFootprintBackground=False

'''Include threshold relative to thresholdValue
	Valid Range = [0.0,inf) '''
root.calibrate.detection.includeThresholdMultiplier=10.0

'''Pixels should be grown as isotropically as possible (slower) '''
root.calibrate.detection.isotropicGrow=False

'''Fiddle factor to add to the background; debugging only '''
root.calibrate.detection.adjustBackground=0.0

'''specifies the desired flavor of Threshold
Allowed values:
	pixel_stdev	threshold applied to per-pixel std deviation
	variance	threshold applied to image variance
	value	threshold applied to image value
	stdev	threshold applied to image std deviation
 '''
root.calibrate.detection.thresholdType='stdev'

'''Names of mask planes to ignore while estimating the background '''
root.calibrate.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

'''behaviour if there are too few points in grid for requested interpolation style
Allowed values:
	None	Field is optional
	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
	THROW_EXCEPTION	throw an exception if there are too few points
	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
 '''
root.calibrate.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

'''how to interpolate the background values. This maps to an enum; see afw::math::Background
Allowed values:
	None	Field is optional
	CONSTANT	Use a single constant value
	LINEAR	Use linear interpolation
	NONE	No background estimation is to be attempted
	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
 '''
root.calibrate.detection.background.algorithm='NATURAL_SPLINE'

'''how large a region of the sky should be used for each background point
	Valid Range = [10,inf) '''
root.calibrate.detection.background.binSize=128

'''Ignore NaNs when estimating the background '''
root.calibrate.detection.background.isNanSafe=False

'''type of statistic to use for grid points
Allowed values:
	None	Field is optional
	MEDIAN	median
	MEANCLIP	clipped mean
	MEAN	unclipped mean
 '''
root.calibrate.detection.background.statisticsProperty='MEANCLIP'

'''Apprimation order for background Chebyshev (valid only with useApprox=True) '''
root.calibrate.detection.background.approxOrder=6

'''Use Approximate (Chebyshev) to model background. '''
root.calibrate.detection.background.useApprox=True

'''Grow detections to set the image mask bits, but return the original (not-grown) footprints '''
root.calibrate.detection.returnOriginalFootprints=True

'''specifies whether to detect positive, or negative sources, or both
Allowed values:
	positive	detect only positive sources
	negative	detect only negative sources
	both	detect both positive and negative sources
 '''
root.calibrate.detection.thresholdPolarity='positive'

'''Threshold for footprints
	Valid Range = [0.0,inf) '''
root.calibrate.detection.thresholdValue=5.0

'''Compute photometric zeropoint? '''
root.calibrate.doPhotoCal=True

'''Whether to compute quantities related to the Gaussian-weighted shape '''
root.calibrate.initialMeasurement.blendedness.doShape=True

'''Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter) '''
root.calibrate.initialMeasurement.blendedness.doOld=True

'''Whether to compute quantities related to the Gaussian-weighted flux '''
root.calibrate.initialMeasurement.blendedness.doFlux=True

'''Radius factor that sets the maximum extent of the weight function (and hence the flux measurements) '''
root.calibrate.initialMeasurement.blendedness.nSigmaWeightMax=3.0

'''Whether to compute blendedness metrics '''
root.calibrate.initialMeasurement.doBlendedness=False

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.centroider['centroid.sdss'].priority=0.0

'''if the peak's less thatn this insist on binning at least once '''
root.calibrate.initialMeasurement.centroider['centroid.sdss'].peakMin=-1.0

'''fiddle factor for adjusting the binning '''
root.calibrate.initialMeasurement.centroider['centroid.sdss'].wfac=1.5

'''maximum allowed binning '''
root.calibrate.initialMeasurement.centroider['centroid.sdss'].binmax=16

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.centroider['centroid.naive'].priority=0.0

'''FIXME! NEVER DOCUMENTED! '''
root.calibrate.initialMeasurement.centroider['centroid.naive'].background=0.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.centroider['centroid.gaussian'].priority=0.0

root.calibrate.initialMeasurement.centroider.name='centroid.sdss'
'''prefix for all measurement fields '''
root.calibrate.initialMeasurement.prefix='initial.'

'''Largest aperture for which to use the slow, accurate, sinc aperture code '''
root.calibrate.initialMeasurement.algorithms['flux.kron'].maxSincRadius=10.0

'''Number of times to iterate when setting the Kron radius '''
root.calibrate.initialMeasurement.algorithms['flux.kron'].nIterForRadius=1

'''Use the Footprint size as part of initial estimate of Kron radius '''
root.calibrate.initialMeasurement.algorithms['flux.kron'].useFootprintRadius=False

'''Minimum Kron radius (if == 0.0 use PSF's Kron radius) if enforceMinimumRadius. Also functions as fallback aperture radius if set. '''
root.calibrate.initialMeasurement.algorithms['flux.kron'].minimumRadius=0.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['flux.kron'].priority=2.0

'''Multiplier of rms size for aperture used to initially estimate the Kron radius '''
root.calibrate.initialMeasurement.algorithms['flux.kron'].nSigmaForRadius=6.0

'''If true check that the Kron radius exceeds some minimum '''
root.calibrate.initialMeasurement.algorithms['flux.kron'].enforceMinimumRadius=True

'''if true, use existing shape and centroid measurements instead of fitting '''
root.calibrate.initialMeasurement.algorithms['flux.kron'].fixed=False

'''Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K '''
root.calibrate.initialMeasurement.algorithms['flux.kron'].smoothingSigma=-1.0

'''Number of Kron radii for Kron flux '''
root.calibrate.initialMeasurement.algorithms['flux.kron'].nRadiusForFlux=2.5

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['flux.naive'].priority=2.0

'''FIXME! NEVER DOCUMENTED! '''
root.calibrate.initialMeasurement.algorithms['flux.naive'].radius=7.0

'''Shapelet order of inner expansion (0 == Gaussian) '''
root.calibrate.initialMeasurement.algorithms['multishapelet.psf'].innerOrder=2

'''Initial radius of inner component in pixels '''
root.calibrate.initialMeasurement.algorithms['multishapelet.psf'].initialRadius=1.5

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.psf'].priority=2.0

'''Minimum inner radius in pixels. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.psf'].minRadius=0.1

'''outer radius divided by inner radius (fixed) '''
root.calibrate.initialMeasurement.algorithms['multishapelet.psf'].radiusRatio=2.0

'''Minimum axis ratio for ellipse (b/a). '''
root.calibrate.initialMeasurement.algorithms['multishapelet.psf'].minAxisRatio=0.1

'''outer Gaussian peak height divided by inner Gaussian peak height; held fixed in double-Gaussian ellipse fit, then allowed to vary when shapelets coefficients are fit and ellipses are held fixed. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.psf'].peakRatio=0.1

'''Shapelet order of outer expansion (0 == Gaussian) '''
root.calibrate.initialMeasurement.algorithms['multishapelet.psf'].outerOrder=1

'''Use fast approximate exponential (good to ~1E-4) '''
root.calibrate.initialMeasurement.algorithms['multishapelet.psf'].useApproximateExp=False

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['flux.peakLikelihood'].priority=2.0

'''Name of warping kernel (e.g. "lanczos4") used to compute the peak '''
root.calibrate.initialMeasurement.algorithms['flux.peakLikelihood'].warpingKernelName='lanczos4'

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['classification.extendedness'].priority=5.0

'''correction factor for psfFlux error '''
root.calibrate.initialMeasurement.algorithms['classification.extendedness'].psfErrFactor=0.0

'''correction factor for modelFlux error '''
root.calibrate.initialMeasurement.algorithms['classification.extendedness'].modelErrFactor=0.0

'''critical ratio of model to psf flux '''
root.calibrate.initialMeasurement.algorithms['classification.extendedness'].fluxRatio=0.95

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['flags.pixel'].priority=0.0

'''List of mask planes for which to search entire footprint '''
root.calibrate.initialMeasurement.algorithms['flags.pixel'].any=[]

'''List of mask planes for which to search center of footprint '''
root.calibrate.initialMeasurement.algorithms['flags.pixel'].center=[]

'''Root name of the FitProfileAlgorithm dev comoonent fields. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.combo'].devName='multishapelet.dev'

'''Number of pixels to grow the footprint by. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.combo'].growFootprint=5

'''Number of half-light radii used to determine the pixels to fit '''
root.calibrate.initialMeasurement.algorithms['multishapelet.combo'].radiusInputFactor=4.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.combo'].priority=2.6

'''Mask planes that indicate pixels that should be ignored in the fit. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.combo'].badMaskPlanes=['EDGE', 'SAT']

'''Root name of the FitProfileAlgorithm exp component fields. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.combo'].expName='multishapelet.exp'

'''If true, individually weigh pixels using the variance image. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.combo'].usePixelWeights=False

'''Root name of the FitPsfAlgorithm fields. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.combo'].psfName='multishapelet.psf'

'''Use fast approximate exponential (good to ~1E-4) '''
root.calibrate.initialMeasurement.algorithms['multishapelet.combo'].useApproximateExp=False

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['shape.hsm.moments'].priority=1.0

'''Mask planes used to reject bad pixels. '''
root.calibrate.initialMeasurement.algorithms['shape.hsm.moments'].badMaskPlanes=[]

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['flux.sinc'].priority=2.0

'''major axis of inner boundary (pixels) '''
root.calibrate.initialMeasurement.algorithms['flux.sinc'].radius1=0.0

'''major axis of outer boundary (pixels) '''
root.calibrate.initialMeasurement.algorithms['flux.sinc'].radius2=7.0

'''measured from x anti-clockwise; radians '''
root.calibrate.initialMeasurement.algorithms['flux.sinc'].angle=0.0

'''1 - b/a '''
root.calibrate.initialMeasurement.algorithms['flux.sinc'].ellipticity=0.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['jacobian'].priority=3.0

'''Nominal pixel size (arcsec) '''
root.calibrate.initialMeasurement.algorithms['jacobian'].pixelScale=0.5

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['shape.hsm.regauss'].priority=1.0

'''Mask planes used to reject bad pixels. '''
root.calibrate.initialMeasurement.algorithms['shape.hsm.regauss'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

'''Field name for number of deblend children '''
root.calibrate.initialMeasurement.algorithms['shape.hsm.regauss'].deblendNChild=''

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['skycoord'].priority=5.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['flux.psf'].priority=2.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['countInputs'].priority=2.0

'''Name of a registered multi-Gaussian profile. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.dev'].profile='tractor-devaucouleur'

'''Minimum half-light radius in units of PSF inner radius for initial parameters. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.dev'].minInitialRadius=None

'''Number of pixels to grow the footprint by. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.dev'].growFootprint=None

'''Number of half-light radii used to determine the pixels to fit '''
root.calibrate.initialMeasurement.algorithms['multishapelet.dev'].radiusInputFactor=None

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.dev'].priority=None

'''Mask planes that indicate pixels that should be ignored in the fit. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.dev'].badMaskPlanes=None

'''Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.dev'].maxBadPixelFraction=None

'''Attempt to approximately deconvolve the canonical shape before using it to set the initial parameters. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.dev'].deconvolveShape=None

'''Minimum axis ratio for ellipse (b/a). '''
root.calibrate.initialMeasurement.algorithms['multishapelet.dev'].minAxisRatio=None

'''If true, individually weigh pixels using the variance image. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.dev'].usePixelWeights=None

'''Use fast approximate exponential (good to ~1E-4) '''
root.calibrate.initialMeasurement.algorithms['multishapelet.dev'].useApproximateExp=None

'''Root name of the FitPsfAlgorithm fields. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.dev'].psfName=None

'''Minimum half-light radius in units of PSF inner radius. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.dev'].minRadius=None

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['flux.aperture'].priority=2.0

'''Maximum number of radial annuli to measure '''
root.calibrate.initialMeasurement.algorithms['flux.aperture'].nApertureMax=10

'''vector of radii for apertures (in pixels) '''
root.calibrate.initialMeasurement.algorithms['flux.aperture'].radii=[1.0, 1.5625, 2.44140625, 3.814697265625, 5.9604644775390625, 9.313225746154785, 14.551915228366852, 22.737367544323206, 35.52713678800501, 55.51115123125783]

'''Largest aperture for which to use the slow, accurate, sinc aperture code '''
root.calibrate.initialMeasurement.algorithms['flux.aperture'].maxSincRadius=10.0

'''Name of a registered multi-Gaussian profile. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.exp'].profile='tractor-exponential'

'''Minimum half-light radius in units of PSF inner radius for initial parameters. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.exp'].minInitialRadius=None

'''Number of pixels to grow the footprint by. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.exp'].growFootprint=None

'''Number of half-light radii used to determine the pixels to fit '''
root.calibrate.initialMeasurement.algorithms['multishapelet.exp'].radiusInputFactor=None

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.exp'].priority=None

'''Mask planes that indicate pixels that should be ignored in the fit. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.exp'].badMaskPlanes=None

'''Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.exp'].maxBadPixelFraction=None

'''Attempt to approximately deconvolve the canonical shape before using it to set the initial parameters. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.exp'].deconvolveShape=None

'''Minimum axis ratio for ellipse (b/a). '''
root.calibrate.initialMeasurement.algorithms['multishapelet.exp'].minAxisRatio=None

'''If true, individually weigh pixels using the variance image. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.exp'].usePixelWeights=None

'''Use fast approximate exponential (good to ~1E-4) '''
root.calibrate.initialMeasurement.algorithms['multishapelet.exp'].useApproximateExp=None

'''Root name of the FitPsfAlgorithm fields. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.exp'].psfName=None

'''Minimum half-light radius in units of PSF inner radius. '''
root.calibrate.initialMeasurement.algorithms['multishapelet.exp'].minRadius=None

'''suffix of shape field flag to check if fixed is true '''
root.calibrate.initialMeasurement.algorithms['flux.gaussian'].shapeFlag='.flags'

'''Convergence tolerance for FWHM '''
root.calibrate.initialMeasurement.algorithms['flux.gaussian'].tol2=9.999999747378752e-05

'''Convergence tolerance for e1,e2 '''
root.calibrate.initialMeasurement.algorithms['flux.gaussian'].tol1=9.999999747378752e-06

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['flux.gaussian'].priority=2.0

'''name of shape field to use if fixed is true '''
root.calibrate.initialMeasurement.algorithms['flux.gaussian'].shape='initial.shape.sdss'

'''name of centroid field to use if fixed is true '''
root.calibrate.initialMeasurement.algorithms['flux.gaussian'].centroid='initial.shape.sdss.centroid'

'''FIXME! NEVER DOCUMENTED! '''
root.calibrate.initialMeasurement.algorithms['flux.gaussian'].background=0.0

'''Maximum number of iterations '''
root.calibrate.initialMeasurement.algorithms['flux.gaussian'].maxIter=100

'''if true, use existing shape and centroid measurements instead of fitting '''
root.calibrate.initialMeasurement.algorithms['flux.gaussian'].fixed=True

'''FIXME! NEVER DOCUMENTED! '''
root.calibrate.initialMeasurement.algorithms['flux.gaussian'].shiftmax=10.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['shape.hsm.psfMoments'].priority=1.0

'''Convergence tolerance for FWHM '''
root.calibrate.initialMeasurement.algorithms['shape.sdss'].tol2=9.999999747378752e-05

'''Convergence tolerance for e1,e2 '''
root.calibrate.initialMeasurement.algorithms['shape.sdss'].tol1=9.999999747378752e-06

'''Whether to also compute the shape of the PSF model '''
root.calibrate.initialMeasurement.algorithms['shape.sdss'].doMeasurePsf=True

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['shape.sdss'].priority=1.0

'''Additional value to add to background '''
root.calibrate.initialMeasurement.algorithms['shape.sdss'].background=0.0

'''Maximum number of iterations '''
root.calibrate.initialMeasurement.algorithms['shape.sdss'].maxIter=100

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['flux.scaled'].priority=2.0

'''scaling factor of PSF FWHM for aperture radius '''
root.calibrate.initialMeasurement.algorithms['flux.scaled'].scale=3.14

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['centroid.sdss'].priority=0.0

'''if the peak's less thatn this insist on binning at least once '''
root.calibrate.initialMeasurement.algorithms['centroid.sdss'].peakMin=-1.0

'''fiddle factor for adjusting the binning '''
root.calibrate.initialMeasurement.algorithms['centroid.sdss'].wfac=1.5

'''maximum allowed binning '''
root.calibrate.initialMeasurement.algorithms['centroid.sdss'].binmax=16

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['centroid.record'].priority=5.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['shape.hsm.linear'].priority=1.0

'''Mask planes used to reject bad pixels. '''
root.calibrate.initialMeasurement.algorithms['shape.hsm.linear'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

'''Field name for number of deblend children '''
root.calibrate.initialMeasurement.algorithms['shape.hsm.linear'].deblendNChild=''

'''Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].minInitialRadius=0.1

'''Use this multiple of the initial fit ellipse then grow by the PSF width to determine the maximum final fit region size. '''
root.calibrate.initialMeasurement.algorithms['cmodel'].region.nFitRadiiMax=3.0

'''Use this multiple of the Kron ellipse to set the fit region (for the final fit region, subject to the nFitRadiiMin and nFitRadiiMax constraints). '''
root.calibrate.initialMeasurement.algorithms['cmodel'].region.nKronRadii=1.5

'''Mask planes that indicate pixels that should be ignored in the fit. '''
root.calibrate.initialMeasurement.algorithms['cmodel'].region.badMaskPlanes=['EDGE', 'SAT', 'BAD', 'NO_DATA']

'''If the Kron radius is less than this multiple of the PSF width, ignore it and fall back to a PSF-oriented ellipse scaled to match the area of the footprint or this radius (whichever is larger). '''
root.calibrate.initialMeasurement.algorithms['cmodel'].region.nPsfSigmaMin=4.0

'''Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try. '''
root.calibrate.initialMeasurement.algorithms['cmodel'].region.maxBadPixelFraction=0.1

'''Use this multiple of the initial fit ellipse then grow by the PSF width to determine the minimum final fit region size. '''
root.calibrate.initialMeasurement.algorithms['cmodel'].region.nFitRadiiMin=1.0

'''Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region '''
root.calibrate.initialMeasurement.algorithms['cmodel'].region.nPsfSigmaGrow=2.0

'''Abort if the fit region grows beyond this many pixels. '''
root.calibrate.initialMeasurement.algorithms['cmodel'].region.maxArea=100000

'''If the maximum of the gradient falls below this threshold, consider the algorithm converged '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.gradientThreshold=0.01

'''steps with reduction radio less than this will decrease the trust radius '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.trustRegionShrinkReductionRatio=0.25

'''when increase the trust region size, multiply the radius by this factor '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.trustRegionGrowFactor=2.0

'''steps with length this fraction of the trust radius may increase the trust radius '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.trustRegionGrowStepFraction=0.8

'''whether to save all iterations for debugging purposes '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.doSaveIterations=False

'''steps with reduction radio greater than this may increase the trust radius '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.trustRegionGrowReductionRatio=0.75

'''when reducing the trust region size, multiply the radius by this factor '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.trustRegionShrinkFactor=0.3333333333333333

'''steps with reduction ratio greater than this are accepted '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.stepAcceptThreshold=0.0

'''the initial trust region will be set to this value '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.trustRegionInitialSize=1.0

'''If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.noSR1Term=False

'''relative step size used for numerical derivatives (added to other steps) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.numDiffRelStep=0.0

'''step size (in units of trust radius) used for numerical derivatives (added to relative step) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.numDiffTrustRadiusStep=0.1

'''maximum number of steps '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.maxOuterIterations=500

'''Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.skipSR1UpdateThreshold=1e-08

'''maximum number of iterations (i.e. function evaluations and trust region subproblems) per step '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.maxInnerIterations=20

'''If the trust radius falls below this threshold, consider the algorithm converged '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.minTrustRadiusThreshold=0.01

'''value passed as the tolerance to solveTrustRegion '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.trustRegionSolverTolerance=1e-08

'''absolute step size used for numerical derivatives (added to other steps) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.optimizer.numDiffAbsStep=0.0

'''Number of degrees of freedom for the Student's T distribution on ln(radius) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.empiricalPriorConfig.logRadiusNu=50.0

'''Width of the Student's T distribution in ln(radius). '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.empiricalPriorConfig.logRadiusSigma=0.45

'''Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T. '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.empiricalPriorConfig.logRadiusMu=-1.0

'''Width of exponential ellipticity distribution (conformal shear units) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.empiricalPriorConfig.ellipticitySigma=0.3

'''Softened core width for ellipticity distribution (conformal shear units '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.empiricalPriorConfig.ellipticityCore=0.001

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.empiricalPriorConfig.logRadiusMinInner=-6.0

'''Minimum ln(radius) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.empiricalPriorConfig.logRadiusMinOuter=-6.001

'''One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.priorSource='EMPIRICAL'

'''Number of Gaussian used to approximate the profile '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.nComponents=3

'''ln(radius) at which the softened cutoff begins towards the maximum '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.linearPriorConfig.logRadiusMaxInner=3.0

'''The ratio P(logRadiusMinInner)/P(logRadiusMaxInner) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.linearPriorConfig.logRadiusMinMaxRatio=1.0

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.linearPriorConfig.logRadiusMinInner=-6.0

'''Maximum ln(radius) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.linearPriorConfig.logRadiusMaxOuter=3.001

'''Ellipticity magnitude (conformal shear units) at which the softened cutoff begins '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.linearPriorConfig.ellipticityMaxInner=2.0

'''Minimum ln(radius) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.linearPriorConfig.logRadiusMinOuter=-6.001

'''Maximum ellipticity magnitude (conformal shear units) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.linearPriorConfig.ellipticityMaxOuter=2.001

'''Name of the Prior that defines the model to fit (a filename in $MEAS_MULTIFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting. '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.priorName=''

'''Whether to record the time spent in this stage '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.doRecordTime=True

'''Name of the shapelet.RadialProfile that defines the model to fit '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.profileName='lux'

'''Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.usePixelWeights=True

'''Maximum radius used in approximating profile with Gaussians (0=default for this profile) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.maxRadius=0

'''Whether to record the steps the optimizer takes (or just the number, if running as a plugin) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].initial.doRecordHistory=True

'''If the maximum of the gradient falls below this threshold, consider the algorithm converged '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.gradientThreshold=1e-05

'''steps with reduction radio less than this will decrease the trust radius '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.trustRegionShrinkReductionRatio=0.25

'''when increase the trust region size, multiply the radius by this factor '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.trustRegionGrowFactor=2.0

'''steps with length this fraction of the trust radius may increase the trust radius '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.trustRegionGrowStepFraction=0.8

'''whether to save all iterations for debugging purposes '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.doSaveIterations=False

'''steps with reduction radio greater than this may increase the trust radius '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.trustRegionGrowReductionRatio=0.75

'''when reducing the trust region size, multiply the radius by this factor '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.trustRegionShrinkFactor=0.3333333333333333

'''steps with reduction ratio greater than this are accepted '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.stepAcceptThreshold=0.0

'''the initial trust region will be set to this value '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.trustRegionInitialSize=1.0

'''If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.noSR1Term=False

'''relative step size used for numerical derivatives (added to other steps) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.numDiffRelStep=0.0

'''step size (in units of trust radius) used for numerical derivatives (added to relative step) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.numDiffTrustRadiusStep=0.1

'''maximum number of steps '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.maxOuterIterations=500

'''Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.skipSR1UpdateThreshold=1e-08

'''maximum number of iterations (i.e. function evaluations and trust region subproblems) per step '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.maxInnerIterations=20

'''If the trust radius falls below this threshold, consider the algorithm converged '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.minTrustRadiusThreshold=1e-05

'''value passed as the tolerance to solveTrustRegion '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.trustRegionSolverTolerance=1e-08

'''absolute step size used for numerical derivatives (added to other steps) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.optimizer.numDiffAbsStep=0.0

'''Number of degrees of freedom for the Student's T distribution on ln(radius) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.empiricalPriorConfig.logRadiusNu=50.0

'''Width of the Student's T distribution in ln(radius). '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.empiricalPriorConfig.logRadiusSigma=0.45

'''Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T. '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.empiricalPriorConfig.logRadiusMu=-1.0

'''Width of exponential ellipticity distribution (conformal shear units) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.empiricalPriorConfig.ellipticitySigma=0.3

'''Softened core width for ellipticity distribution (conformal shear units '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.empiricalPriorConfig.ellipticityCore=0.001

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.empiricalPriorConfig.logRadiusMinInner=-6.0

'''Minimum ln(radius) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.empiricalPriorConfig.logRadiusMinOuter=-6.001

'''One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.priorSource='EMPIRICAL'

'''Number of Gaussian used to approximate the profile '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.nComponents=8

'''ln(radius) at which the softened cutoff begins towards the maximum '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.linearPriorConfig.logRadiusMaxInner=3.0

'''The ratio P(logRadiusMinInner)/P(logRadiusMaxInner) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.linearPriorConfig.logRadiusMinMaxRatio=1.0

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.linearPriorConfig.logRadiusMinInner=-6.0

'''Maximum ln(radius) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.linearPriorConfig.logRadiusMaxOuter=3.001

'''Ellipticity magnitude (conformal shear units) at which the softened cutoff begins '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.linearPriorConfig.ellipticityMaxInner=2.0

'''Minimum ln(radius) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.linearPriorConfig.logRadiusMinOuter=-6.001

'''Maximum ellipticity magnitude (conformal shear units) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.linearPriorConfig.ellipticityMaxOuter=2.001

'''Name of the Prior that defines the model to fit (a filename in $MEAS_MULTIFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting. '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.priorName=''

'''Whether to record the time spent in this stage '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.doRecordTime=True

'''Name of the shapelet.RadialProfile that defines the model to fit '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.profileName='luv'

'''Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.usePixelWeights=False

'''Maximum radius used in approximating profile with Gaussians (0=default for this profile) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.maxRadius=0

'''Whether to record the steps the optimizer takes (or just the number, if running as a plugin) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].dev.doRecordHistory=True

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['cmodel'].priority=2.5

'''If the maximum of the gradient falls below this threshold, consider the algorithm converged '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.gradientThreshold=1e-05

'''steps with reduction radio less than this will decrease the trust radius '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.trustRegionShrinkReductionRatio=0.25

'''when increase the trust region size, multiply the radius by this factor '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.trustRegionGrowFactor=2.0

'''steps with length this fraction of the trust radius may increase the trust radius '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.trustRegionGrowStepFraction=0.8

'''whether to save all iterations for debugging purposes '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.doSaveIterations=False

'''steps with reduction radio greater than this may increase the trust radius '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.trustRegionGrowReductionRatio=0.75

'''when reducing the trust region size, multiply the radius by this factor '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.trustRegionShrinkFactor=0.3333333333333333

'''steps with reduction ratio greater than this are accepted '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.stepAcceptThreshold=0.0

'''the initial trust region will be set to this value '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.trustRegionInitialSize=1.0

'''If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.noSR1Term=False

'''relative step size used for numerical derivatives (added to other steps) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.numDiffRelStep=0.0

'''step size (in units of trust radius) used for numerical derivatives (added to relative step) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.numDiffTrustRadiusStep=0.1

'''maximum number of steps '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.maxOuterIterations=250

'''Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.skipSR1UpdateThreshold=1e-08

'''maximum number of iterations (i.e. function evaluations and trust region subproblems) per step '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.maxInnerIterations=20

'''If the trust radius falls below this threshold, consider the algorithm converged '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.minTrustRadiusThreshold=1e-05

'''value passed as the tolerance to solveTrustRegion '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.trustRegionSolverTolerance=1e-08

'''absolute step size used for numerical derivatives (added to other steps) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.optimizer.numDiffAbsStep=0.0

'''Number of degrees of freedom for the Student's T distribution on ln(radius) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.empiricalPriorConfig.logRadiusNu=50.0

'''Width of the Student's T distribution in ln(radius). '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.empiricalPriorConfig.logRadiusSigma=0.45

'''Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T. '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.empiricalPriorConfig.logRadiusMu=-1.0

'''Width of exponential ellipticity distribution (conformal shear units) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.empiricalPriorConfig.ellipticitySigma=0.3

'''Softened core width for ellipticity distribution (conformal shear units '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.empiricalPriorConfig.ellipticityCore=0.001

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.empiricalPriorConfig.logRadiusMinInner=-6.0

'''Minimum ln(radius) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.empiricalPriorConfig.logRadiusMinOuter=-6.001

'''One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.priorSource='EMPIRICAL'

'''Number of Gaussian used to approximate the profile '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.nComponents=6

'''ln(radius) at which the softened cutoff begins towards the maximum '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.linearPriorConfig.logRadiusMaxInner=3.0

'''The ratio P(logRadiusMinInner)/P(logRadiusMaxInner) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.linearPriorConfig.logRadiusMinMaxRatio=1.0

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.linearPriorConfig.logRadiusMinInner=-6.0

'''Maximum ln(radius) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.linearPriorConfig.logRadiusMaxOuter=3.001

'''Ellipticity magnitude (conformal shear units) at which the softened cutoff begins '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.linearPriorConfig.ellipticityMaxInner=2.0

'''Minimum ln(radius) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.linearPriorConfig.logRadiusMinOuter=-6.001

'''Maximum ellipticity magnitude (conformal shear units) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.linearPriorConfig.ellipticityMaxOuter=2.001

'''Name of the Prior that defines the model to fit (a filename in $MEAS_MULTIFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting. '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.priorName=''

'''Whether to record the time spent in this stage '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.doRecordTime=True

'''Name of the shapelet.RadialProfile that defines the model to fit '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.profileName='lux'

'''Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.usePixelWeights=False

'''Maximum radius used in approximating profile with Gaussians (0=default for this profile) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.maxRadius=0

'''Whether to record the steps the optimizer takes (or just the number, if running as a plugin) '''
root.calibrate.initialMeasurement.algorithms['cmodel'].exp.doRecordHistory=True

'''If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this.  If <= 0.0, abort the fit early instead. '''
root.calibrate.initialMeasurement.algorithms['cmodel'].fallbackInitialMomentsPsfFactor=1.5

'''Root name of the FitPsfAlgorithm fields. '''
root.calibrate.initialMeasurement.algorithms['cmodel'].psfName='multishapelet.psf'

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['focalplane'].priority=3.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['centroid.naive'].priority=0.0

'''FIXME! NEVER DOCUMENTED! '''
root.calibrate.initialMeasurement.algorithms['centroid.naive'].background=0.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['shape.hsm.bj'].priority=1.0

'''Mask planes used to reject bad pixels. '''
root.calibrate.initialMeasurement.algorithms['shape.hsm.bj'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

'''Field name for number of deblend children '''
root.calibrate.initialMeasurement.algorithms['shape.hsm.bj'].deblendNChild=''

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['centroid.gaussian'].priority=0.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['flux.aperture.elliptical'].priority=1.899999976158142

'''Maximum number of radial annuli to measure '''
root.calibrate.initialMeasurement.algorithms['flux.aperture.elliptical'].nApertureMax=10

'''vector of radii for apertures (in pixels) '''
root.calibrate.initialMeasurement.algorithms['flux.aperture.elliptical'].radii=[1.0, 1.5625, 2.44140625, 3.814697265625, 5.9604644775390625, 9.313225746154785, 14.551915228366852, 22.737367544323206, 35.52713678800501, 55.51115123125783]

'''Largest aperture for which to use the slow, accurate, sinc aperture code '''
root.calibrate.initialMeasurement.algorithms['flux.aperture.elliptical'].maxSincRadius=10.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['variance'].priority=2.0

'''scale factor to apply to shape for aperture '''
root.calibrate.initialMeasurement.algorithms['variance'].scale=5.0

'''mask planes to ignore '''
root.calibrate.initialMeasurement.algorithms['variance'].mask=['DETECTED', 'DETECTED_NEGATIVE']

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['correctfluxes'].priority=3.0

'''List of flux fields that should not be corrected (otherwise all fields in getApCorrRegistry() will be) '''
root.calibrate.initialMeasurement.algorithms['correctfluxes'].ignored=[]

'''Whether to propagate aperture correction uncertainties into flux uncertainties '''
root.calibrate.initialMeasurement.algorithms['correctfluxes'].doPropagateErrors=False

'''Whether to set the general failure flag for a flux when it cannot be aperture-corrected '''
root.calibrate.initialMeasurement.algorithms['correctfluxes'].doFlagApCorrFailures=True

'''Whether to save the per-source per-flux aperture corrections and their errors '''
root.calibrate.initialMeasurement.algorithms['correctfluxes'].doRecordApCorr=True

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.initialMeasurement.algorithms['shape.hsm.ksb'].priority=1.0

'''Mask planes used to reject bad pixels. '''
root.calibrate.initialMeasurement.algorithms['shape.hsm.ksb'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

'''Field name for number of deblend children '''
root.calibrate.initialMeasurement.algorithms['shape.hsm.ksb'].deblendNChild=''

root.calibrate.initialMeasurement.algorithms.names=['flux.psf', 'flags.pixel', 'flux.naive', 'flux.gaussian', 'centroid.naive', 'flux.sinc', 'shape.sdss', 'skycoord', 'classification.extendedness']
'''The seed value to use for random number generation. '''
root.calibrate.initialMeasurement.replaceWithNoise.noiseSeed=0

'''Add ann offset to the generated noise. '''
root.calibrate.initialMeasurement.replaceWithNoise.noiseOffset=0.0

'''How do we choose the mean and variance of the Gaussian noise we generate?
Allowed values:
	variance	Mean = 0, variance = the image's variance
	meta	Mean = 0, variance = the "BGMEAN" metadata entry
	measure	Measure clipped mean and variance from the whole image
 '''
root.calibrate.initialMeasurement.replaceWithNoise.noiseSource='measure'

'''the name of the aperture flux algorithm used for calibration '''
root.calibrate.initialMeasurement.slots.calibFlux='flux.naive'

'''the name of the algorithm used to set the source aperture flux slot '''
root.calibrate.initialMeasurement.slots.apFlux='flux.sinc'

'''the name of the algorithm used to set the source inst flux slot '''
root.calibrate.initialMeasurement.slots.instFlux='flux.gaussian'

'''the name of the algorithm used to set source moments parameters '''
root.calibrate.initialMeasurement.slots.shape='shape.sdss'

'''the name of the centroiding algorithm used to set source x,y '''
root.calibrate.initialMeasurement.slots.centroid='centroid.sdss'

'''the name of the algorithm used to set the source model flux slot '''
root.calibrate.initialMeasurement.slots.modelFlux='flux.gaussian'

'''the name of the algorithm used to set the source psf flux slot '''
root.calibrate.initialMeasurement.slots.psfFlux='flux.psf'

'''When measuring, replace other detected footprints with noise? '''
root.calibrate.initialMeasurement.doReplaceWithNoise=True

'''List of values of reduced chi^2 that should be applied in order to clip sources

        N.b. these values are so large because of contamination in the annuli by to e.g. the faint wings of
        neighbouring objects.  The real solution here is to write cleverer aperture flux code, 'a la SDSS
         '''
root.calibrate.measureCurveOfGrowth.maxRChi2=[10000.0, 1000.0]

'''After perfoming an ML estimation of the curve of growth to estimate the relative
        weights we may undertake a per-radial-point re-estimation of the curve of growth;
        finalEstimationAlgorithm dictates how to estimate the value of each point.
        
Allowed values:
	None	no clipping
	median	median
	mean	weighted mean
 '''
root.calibrate.measureCurveOfGrowth.finalEstimationAlgorithm='mean'

'''Maximum fraction of the pixels in an annulus which may be labelled INTRP

        Annuli with more than fracInterpolatedMax interpolated pixels are rejected
         '''
root.calibrate.measureCurveOfGrowth.fracInterpolatedMax=0.1

'''Maximum number of aperture fluxes to use (all, if None) '''
root.calibrate.measureCurveOfGrowth.nAperture=8

'''How many of the brightest saturated candidates to use (0 => all) '''
root.calibrate.measureCurveOfGrowth.nSaturated=500

'''Minimum value of psf.flux for an object to be included in curve of growth '''
root.calibrate.measureCurveOfGrowth.psfFluxMin=50000.0

'''Maximum value of classification.extendedness for an object to be included in curve of growth '''
root.calibrate.measureCurveOfGrowth.classificationMax=float('nan')

'''List of flags which cause a source to be rejected as bad

        N.b. may contain globs such as *.flags.pixel.edge '''
root.calibrate.measureCurveOfGrowth.badFlags=['*.flags.pixel.edge']

'''Minimum flux measured in an annulus to be included in the curve of growth

        Annuli with fluxes of less than minAnnularFlux are rejected
         '''
root.calibrate.measureCurveOfGrowth.minAnnularFlux=0.0

'''How many of the brightest non-saturated candidates to use (0 => all) '''
root.calibrate.measureCurveOfGrowth.nNonSaturated=200

'''The per-pixel std. dev. of our knowledge of the sky, to be added in quadrature '''
root.calibrate.measureCurveOfGrowth.skyNoiseFloor=0.0

'''Compute astrometric solution? '''
root.calibrate.doAstrometry=True

'''floor for variance is lam*data '''
root.calibrate.measurePsf.psfDeterminer['pca'].lam=0.05

'''for psf candidate evaluation '''
root.calibrate.measurePsf.psfDeterminer['pca'].reducedChi2ForPsfCandidates=2.0

'''Use non-linear fitter for spatial variation of Kernel '''
root.calibrate.measurePsf.psfDeterminer['pca'].nonLinearSpatialFit=False

'''Mask blends in image? '''
root.calibrate.measurePsf.psfDeterminer['pca'].doMaskBlends=True

'''size of cell used to determine PSF (pixels, column direction) '''
root.calibrate.measurePsf.psfDeterminer['pca'].sizeCellX=256

'''size of cell used to determine PSF (pixels, row direction) '''
root.calibrate.measurePsf.psfDeterminer['pca'].sizeCellY=256

'''Rejection threshold (stdev) for candidates based on spatial fit '''
root.calibrate.measurePsf.psfDeterminer['pca'].spatialReject=3.0

'''radius of the kernel to create, relative to the square root of the stellar quadrupole moments '''
root.calibrate.measurePsf.psfDeterminer['pca'].kernelSize=10.0

'''Reject candidates that are blended? '''
root.calibrate.measurePsf.psfDeterminer['pca'].doRejectBlends=False

'''number of iterations of PSF candidate star list '''
root.calibrate.measurePsf.psfDeterminer['pca'].nIterForPsf=3

'''Should each PSF candidate be given the same weight, independent of magnitude? '''
root.calibrate.measurePsf.psfDeterminer['pca'].constantWeight=True

'''Maximum radius of the kernel '''
root.calibrate.measurePsf.psfDeterminer['pca'].kernelSizeMax=45

'''number of stars per psf Cell for spatial fitting '''
root.calibrate.measurePsf.psfDeterminer['pca'].nStarPerCellSpatialFit=5

'''Number of pixels to ignore around the edge of PSF candidate postage stamps '''
root.calibrate.measurePsf.psfDeterminer['pca'].borderWidth=0

'''number of eigen components for PSF kernel creation '''
root.calibrate.measurePsf.psfDeterminer['pca'].nEigenComponents=4

'''Threshold (stdev) for rejecting extraneous pixels around candidate; applied if positive '''
root.calibrate.measurePsf.psfDeterminer['pca'].pixelThreshold=0.0

'''specify spatial order for PSF kernel creation '''
root.calibrate.measurePsf.psfDeterminer['pca'].spatialOrder=2

'''tolerance of spatial fitting '''
root.calibrate.measurePsf.psfDeterminer['pca'].tolerance=0.01

'''Minimum radius of the kernel '''
root.calibrate.measurePsf.psfDeterminer['pca'].kernelSizeMin=25

'''number of stars per psf cell for PSF kernel creation '''
root.calibrate.measurePsf.psfDeterminer['pca'].nStarPerCell=3

'''floor for variance is lam*data '''
root.calibrate.measurePsf.psfDeterminer['psfex'].lam=0.05

'''for psf candidate evaluation '''
root.calibrate.measurePsf.psfDeterminer['psfex'].reducedChi2ForPsfCandidates=2.0

'''specify spatial order for PSF kernel creation '''
root.calibrate.measurePsf.psfDeterminer['psfex'].spatialOrder=2

'''number of eigen components for PSF kernel creation '''
root.calibrate.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nEigenComponents=4

'''Should PSFEX be permitted to recentroid PSF candidates? '''
root.calibrate.measurePsf.psfDeterminer['psfex'].recentroid=False

'''Resolution of the internal PSF model relative to the pixel size; e.g. 0.5 is equal to 2x oversampling '''
root.calibrate.measurePsf.psfDeterminer['psfex'].samplingSize=1.0

'''size of cell used to determine PSF (pixels, row direction) '''
root.calibrate.measurePsf.psfDeterminer['psfex'].sizeCellY=256

'''List of mask bits which cause a source to be rejected as bad
N.b. INTRP is used specially in PsfCandidateSet; it means "Contaminated by neighbour"
 '''
root.calibrate.measurePsf.psfDeterminer['psfex'].badMaskBits=['INTRP', 'SAT']

'''Rejection threshold (stdev) for candidates based on spatial fit '''
root.calibrate.measurePsf.psfDeterminer['psfex'].spatialReject=3.0

'''Number of pixels to ignore around the edge of PSF candidate postage stamps '''
root.calibrate.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__borderWidth=0

'''radius of the kernel to create, relative to the square root of the stellar quadrupole moments '''
root.calibrate.measurePsf.psfDeterminer['psfex'].kernelSize=41.0

'''size of cell used to determine PSF (pixels, column direction) '''
root.calibrate.measurePsf.psfDeterminer['psfex'].sizeCellX=256

'''Should each PSF candidate be given the same weight, independent of magnitude? '''
root.calibrate.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__constantWeight=True

'''Maximum radius of the kernel '''
root.calibrate.measurePsf.psfDeterminer['psfex'].kernelSizeMax=45

'''number of stars per psf Cell for spatial fitting '''
root.calibrate.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nStarPerCellSpatialFit=5

'''number of stars per psf cell for PSF kernel creation '''
root.calibrate.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nStarPerCell=3

'''number of iterations of PSF candidate star list '''
root.calibrate.measurePsf.psfDeterminer['psfex']._PsfexPsfDeterminerConfig__nIterForPsf=3

'''tolerance of spatial fitting '''
root.calibrate.measurePsf.psfDeterminer['psfex'].tolerance=0.01

'''Minimum radius of the kernel '''
root.calibrate.measurePsf.psfDeterminer['psfex'].kernelSizeMin=25

root.calibrate.measurePsf.psfDeterminer.name='psfex'
'''maximum width to include in histogram '''
root.calibrate.measurePsf.starSelector['objectSize'].widthMax=10.0

'''Standard deviation of width allowed to be interpreted as good stars '''
root.calibrate.measurePsf.starSelector['objectSize'].widthStdAllowed=0.15

'''minimum width to include in histogram '''
root.calibrate.measurePsf.starSelector['objectSize'].widthMin=0.9

'''size of the Psf kernel to create '''
root.calibrate.measurePsf.starSelector['objectSize'].kernelSize=21

'''specify the minimum psfFlux for good Psf Candidates '''
root.calibrate.measurePsf.starSelector['objectSize'].fluxMin=4000.0

'''number of pixels to ignore around the edge of PSF candidate postage stamps '''
root.calibrate.measurePsf.starSelector['objectSize'].borderWidth=0

'''List of flags which cause a source to be rejected as bad '''
root.calibrate.measurePsf.starSelector['objectSize'].badFlags=['initial.flags.pixel.edge', 'initial.flags.pixel.interpolated.center', 'initial.flags.pixel.saturated.center', 'initial.flags.pixel.cr.center', 'initial.flags.pixel.bad', 'initial.flags.pixel.interpolated.any']

'''specify the maximum psfFlux for good Psf Candidates (ignored if == 0) '''
root.calibrate.measurePsf.starSelector['objectSize'].fluxMax=0.0

'''Name of field in Source to use for flux measurement '''
root.calibrate.measurePsf.starSelector['objectSize'].sourceFluxField='initial.flux.psf'

'''Keep objects within this many sigma of cluster 0's median '''
root.calibrate.measurePsf.starSelector['objectSize'].nSigmaClip=2.0

'''Undistort when evaluating the 2nd moments of sources? '''
root.calibrate.measurePsf.starSelector['mitaka'].doUndistort=False

'''How many sigmas around the peak fwhm are used for calculating statistics of PSF sequence '''
root.calibrate.measurePsf.starSelector['mitaka'].fwhmMarginNsigma=3.0

'''Size of grid (pixels) '''
root.calibrate.measurePsf.starSelector['mitaka'].gridSize=1024.0

'''Minimum fwhm allowed in estimation of seeing (pix) '''
root.calibrate.measurePsf.starSelector['mitaka'].fwhmMin=1.5

'''Number of smallest objects which are used to determine rough-interim seeing '''
root.calibrate.measurePsf.starSelector['mitaka'].nSmallSampleRoughFwhm=50

'''How many sigmas around the peak fwhm are used for calculating statistics of PSF sequence '''
root.calibrate.measurePsf.starSelector['mitaka'].psfSeqStatNsigma=3.0

'''Number of bins for number counting as a fn of instrumnetal mag '''
root.calibrate.measurePsf.starSelector['mitaka'].nbinMagHist=80

'''What fraction of sources from the brightest is to be included for initial guess of seeing to avoid cosmic rays which dominate faint magnitudes '''
root.calibrate.measurePsf.starSelector['mitaka'].fracSrcIni=0.15

'''How many times do we iterate calculating statistics of PSF sequence '''
root.calibrate.measurePsf.starSelector['mitaka'].psfSeqStatNiter=3

'''size of the kernel to create '''
root.calibrate.measurePsf.starSelector['mitaka'].kernelSize=21

'''Bin size of FWHM histogram '''
root.calibrate.measurePsf.starSelector['mitaka'].fwhmBinSize=0.2

'''Faintest mag for number counting as a fn of instrumnetal mag '''
root.calibrate.measurePsf.starSelector['mitaka'].magMaxHist=0.0

'''number of pixels to ignore around the edge of PSF candidate postage stamps '''
root.calibrate.measurePsf.starSelector['mitaka'].borderWidth=0

'''Number of brightest (non-saturated) objects which are used to determine rough-interim seeing '''
root.calibrate.measurePsf.starSelector['mitaka'].nBrightSampleRoughFwhm=30

'''Brightest mag for number counting as a fn of instrumnetal mag '''
root.calibrate.measurePsf.starSelector['mitaka'].magMinHist=-20.0

'''How many pixels around the peak are used for calculating scatter of psf candidates '''
root.calibrate.measurePsf.starSelector['mitaka'].fwhmMarginFinal=0.75

'''Make plots? '''
root.calibrate.measurePsf.starSelector['mitaka'].doPlots=True

'''Statistical algorithm to derive rough Fwhm in the 1st step seeing estimation
Allowed values:
	None	Field is optional
	MEDIAN	median of sample
	MEANCLIP	clipped mean of sample with 3-sigma clip + 3-times iteration
 '''
root.calibrate.measurePsf.starSelector['mitaka'].statAlgRoughFwhm='MEDIAN'

'''How many magnitudes to extend the faint-end limit for extracting PSF sources, from the base magnitude determined by fracSrcIni. '''
root.calibrate.measurePsf.starSelector['mitaka'].magLimitFaintExtension=0.0

'''Maxmum fwhm allowed in estimation of seeing (pix) '''
root.calibrate.measurePsf.starSelector['mitaka'].fwhmMax=12.0

'''specify the minimum psfFlux for good Psf Candidates '''
root.calibrate.measurePsf.starSelector['catalog'].fluxLim=0.0

'''specify the maximum psfFlux for good Psf Candidates (ignored if == 0) '''
root.calibrate.measurePsf.starSelector['catalog'].fluxMax=0.0

'''PSF candidate objects may not have any of these bits set '''
root.calibrate.measurePsf.starSelector['catalog'].badStarPixelFlags=['flags.pixel.edge', 'flags.pixel.interpolated.center', 'flags.pixel.saturated.center', 'initial.flags.pixel.edge', 'initial.flags.pixel.interpolated.center', 'initial.flags.pixel.saturated.center']

'''size of the kernel to create '''
root.calibrate.measurePsf.starSelector['catalog'].kernelSize=21

'''number of pixels to ignore around the edge of PSF candidate postage stamps '''
root.calibrate.measurePsf.starSelector['catalog'].borderWidth=0

'''size of the kernel to create '''
root.calibrate.measurePsf.starSelector['secondMoment'].kernelSize=21

'''Multiplier of mean for minimum moments histogram range '''
root.calibrate.measurePsf.starSelector['secondMoment'].histMomentMinMultiplier=2.0

'''Number of bins in moment histogram '''
root.calibrate.measurePsf.starSelector['secondMoment'].histSize=64

'''Clipping threshold for moments histogram range '''
root.calibrate.measurePsf.starSelector['secondMoment'].histMomentClip=5.0

'''number of pixels to ignore around the edge of PSF candidate postage stamps '''
root.calibrate.measurePsf.starSelector['secondMoment'].borderWidth=0

'''List of flags which cause a source to be rejected as bad '''
root.calibrate.measurePsf.starSelector['secondMoment'].badFlags=['initial.flags.pixel.edge', 'initial.flags.pixel.interpolated.center', 'initial.flags.pixel.saturated.center', 'initial.flags.pixel.cr.center']

'''Multiplier of mean for maximum moments histogram range '''
root.calibrate.measurePsf.starSelector['secondMoment'].histMomentMaxMultiplier=5.0

'''specify the maximum psfFlux for good Psf Candidates (ignored if == 0) '''
root.calibrate.measurePsf.starSelector['secondMoment'].fluxMax=0.0

'''candidate PSF's shapes must lie within this many sigma of the average shape '''
root.calibrate.measurePsf.starSelector['secondMoment'].clumpNSigma=2.0

'''specify the minimum psfFlux for good Psf Candidates '''
root.calibrate.measurePsf.starSelector['secondMoment'].fluxLim=12500.0

'''Maximum moment to consider '''
root.calibrate.measurePsf.starSelector['secondMoment'].histMomentMax=100.0

'''Fraction of objects to use in first pass '''
root.calibrate.measurePsf.starSelector['sizeMagnitude'].startn1=0.1

'''Minimum size to use '''
root.calibrate.measurePsf.starSelector['sizeMagnitude'].minsize=0.0

'''Order of polynomial of fit of size(x,y) '''
root.calibrate.measurePsf.starSelector['sizeMagnitude'].fitorder=1

'''Minimum magnitude to use '''
root.calibrate.measurePsf.starSelector['sizeMagnitude'].minmag=0.0

'''What fraction of objects are likely stars? '''
root.calibrate.measurePsf.starSelector['sizeMagnitude'].starfrac=0.5

'''Maximum magnitude to use '''
root.calibrate.measurePsf.starSelector['sizeMagnitude'].maxmag=1e+100

'''Are sizes already log(size)? '''
root.calibrate.measurePsf.starSelector['sizeMagnitude'].logsize=False

'''Perform size(x,y) fit with fitStars brightest stars '''
root.calibrate.measurePsf.starSelector['sizeMagnitude'].starsperbin=30

'''nSigma to reject a star as an outlier '''
root.calibrate.measurePsf.starSelector['sizeMagnitude'].fitsigclip=4.0

'''Maximum size to use '''
root.calibrate.measurePsf.starSelector['sizeMagnitude'].maxsize=1e+100

'''nSigma to reject a star as an outlier '''
root.calibrate.measurePsf.starSelector['sizeMagnitude'].aperture=5.0

'''Smaller = purer smaple of stars, larger = more stars '''
root.calibrate.measurePsf.starSelector['sizeMagnitude'].purityratio=0.05

root.calibrate.measurePsf.starSelector.name='objectSize'
'''This number will be multplied by the exposure ID to set the random seed for reserving candidates '''
root.calibrate.measurePsf.reserveSeed=1

'''Fraction PSF candidates to reserve from fitting '''
root.calibrate.measurePsf.reserveFraction=0.2

'''Require photometric calibration to succeed? '''
root.calibrate.requirePhotoCal=True

'''Compute aperture corrections? '''
root.calibrate.doMeasureApCorr=False

'''number of iterations '''
root.calibrate.photocal.nIter=20

'''Name of the source flux field to use.  The associated flag field
('<name>.flags') will be implicitly included in badFlags.
 '''
root.calibrate.photocal.fluxField='flux.naive'

root.calibrate.photocal.colorterms.library={}
root.calibrate.photocal.colorterms.library['hsc*']=lsst.meas.photocal.colorterms.ColortermGroupConfig()
root.calibrate.photocal.colorterms.library['hsc*'].group={}
root.calibrate.photocal.colorterms.library['hsc*'].group['i']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['hsc*'].group['i'].c2=0.0

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['hsc*'].group['i'].c1=0.0

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['hsc*'].group['i'].c0=0.0

'''Primary band '''
root.calibrate.photocal.colorterms.library['hsc*'].group['i'].primary='i'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['hsc*'].group['i'].secondary='i'

root.calibrate.photocal.colorterms.library['hsc*'].group['y']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['hsc*'].group['y'].c2=0.0

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['hsc*'].group['y'].c1=0.0

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['hsc*'].group['y'].c0=0.0

'''Primary band '''
root.calibrate.photocal.colorterms.library['hsc*'].group['y'].primary='y'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['hsc*'].group['y'].secondary='y'

root.calibrate.photocal.colorterms.library['hsc*'].group['r']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['hsc*'].group['r'].c2=0.0

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['hsc*'].group['r'].c1=0.0

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['hsc*'].group['r'].c0=0.0

'''Primary band '''
root.calibrate.photocal.colorterms.library['hsc*'].group['r'].primary='r'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['hsc*'].group['r'].secondary='r'

root.calibrate.photocal.colorterms.library['hsc*'].group['z']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['hsc*'].group['z'].c2=0.0

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['hsc*'].group['z'].c1=0.0

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['hsc*'].group['z'].c0=0.0

'''Primary band '''
root.calibrate.photocal.colorterms.library['hsc*'].group['z'].primary='z'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['hsc*'].group['z'].secondary='z'

root.calibrate.photocal.colorterms.library['hsc*'].group['g']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['hsc*'].group['g'].c2=0.0

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['hsc*'].group['g'].c1=0.0

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['hsc*'].group['g'].c0=0.0

'''Primary band '''
root.calibrate.photocal.colorterms.library['hsc*'].group['g'].primary='g'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['hsc*'].group['g'].secondary='g'

root.calibrate.photocal.colorterms.library['sdss*']=lsst.meas.photocal.colorterms.ColortermGroupConfig()
root.calibrate.photocal.colorterms.library['sdss*'].group={}
root.calibrate.photocal.colorterms.library['sdss*'].group['r2']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['r2'].c2=-0.0284842

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['r2'].c1=-0.00830543

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['r2'].c0=0.00074087

'''Primary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['r2'].primary='r'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['r2'].secondary='i'

root.calibrate.photocal.colorterms.library['sdss*'].group['g']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['g'].c2=-0.00726883

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['g'].c1=-0.08366937

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['g'].c0=-0.00816446

'''Primary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['g'].primary='g'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['g'].secondary='r'

root.calibrate.photocal.colorterms.library['sdss*'].group['N816']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['N816'].c2=-0.05474862

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['N816'].c1=-0.63558358

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['N816'].c0=0.00927133

'''Primary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['N816'].primary='i'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['N816'].secondary='z'

root.calibrate.photocal.colorterms.library['sdss*'].group['i']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['i'].c2=-0.01374245

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['i'].c1=-0.16922042

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['i'].c0=0.00130204

'''Primary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['i'].primary='i'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['i'].secondary='z'

root.calibrate.photocal.colorterms.library['sdss*'].group['i2']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['i2'].c2=-0.01067212

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['i2'].c1=-0.20739606

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['i2'].c0=0.00124676

'''Primary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['i2'].primary='i'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['i2'].secondary='z'

root.calibrate.photocal.colorterms.library['sdss*'].group['r']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['r'].c2=-0.03068248

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['r'].c1=0.01284177

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['r'].c0=0.0023181

'''Primary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['r'].primary='r'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['r'].secondary='i'

root.calibrate.photocal.colorterms.library['sdss*'].group['N921']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['N921'].c2=-0.05451118

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['N921'].c1=0.0986353

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['N921'].c0=0.00752972

'''Primary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['N921'].primary='z'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['N921'].secondary='i'

root.calibrate.photocal.colorterms.library['sdss*'].group['y']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['y'].c2=0.00574408

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['y'].c1=0.35652971

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['y'].c0=0.01739708

'''Primary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['y'].primary='z'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['y'].secondary='i'

root.calibrate.photocal.colorterms.library['sdss*'].group['z']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['z'].c2=0.01479369

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['z'].c1=0.01353969

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['sdss*'].group['z'].c0=-0.0068062

'''Primary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['z'].primary='z'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['sdss*'].group['z'].secondary='i'

root.calibrate.photocal.colorterms.library['ps1*']=lsst.meas.photocal.colorterms.ColortermGroupConfig()
root.calibrate.photocal.colorterms.library['ps1*'].group={}
root.calibrate.photocal.colorterms.library['ps1*'].group['r2']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['r2'].c2=-0.01667794

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['r2'].c1=3.996e-05

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['r2'].c0=0.0011769

'''Primary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['r2'].primary='r'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['r2'].secondary='i'

root.calibrate.photocal.colorterms.library['ps1*'].group['g']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['g'].c2=-0.0151057

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['g'].c1=0.06508481

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['g'].c0=0.00730066

'''Primary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['g'].primary='g'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['g'].secondary='r'

root.calibrate.photocal.colorterms.library['ps1*'].group['N816']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['N816'].c2=-0.10781564

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['N816'].c1=-0.68757034

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['N816'].c0=0.01191062

'''Primary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['N816'].primary='i'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['N816'].secondary='z'

root.calibrate.photocal.colorterms.library['ps1*'].group['i']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['i'].c2=-0.03034094

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['i'].c1=-0.13944659

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['i'].c0=0.00166891

'''Primary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['i'].primary='i'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['i'].secondary='z'

root.calibrate.photocal.colorterms.library['ps1*'].group['i2']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['i2'].c2=-0.02675511

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['i2'].c1=-0.18483562

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['i2'].c0=0.00180361

'''Primary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['i2'].primary='i'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['i2'].secondary='z'

root.calibrate.photocal.colorterms.library['ps1*'].group['r']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['r'].c2=-0.01877566

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['r'].c1=0.02093734

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['r'].c0=0.00279757

'''Primary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['r'].primary='r'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['r'].secondary='i'

root.calibrate.photocal.colorterms.library['ps1*'].group['N921']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['N921'].c2=-0.25059679

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['N921'].c1=-0.59278367

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['N921'].c0=0.00142051

'''Primary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['N921'].primary='z'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['N921'].secondary='y'

root.calibrate.photocal.colorterms.library['ps1*'].group['y']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['y'].c2=0.02880125

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['y'].c1=0.14747401

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['y'].c0=-0.00156858

'''Primary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['y'].primary='y'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['y'].secondary='z'

root.calibrate.photocal.colorterms.library['ps1*'].group['z']=lsst.meas.photocal.colorterms.ColortermConfig()
'''Second-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['z'].c2=-0.00316369

'''First-order parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['z'].c1=-0.28840221

'''Constant parameter '''
root.calibrate.photocal.colorterms.library['ps1*'].group['z'].c0=-0.00907517

'''Primary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['z'].primary='z'

'''Secondary band '''
root.calibrate.photocal.colorterms.library['ps1*'].group['z'].secondary='y'

'''Don't use objects fainter than this magnitude '''
root.calibrate.photocal.magLimit=22.0

'''maximum sigma to use when clipping '''
root.calibrate.photocal.sigmaMax=0.25

'''List of source flag fields that must be set for a source to be used. '''
root.calibrate.photocal.goodFlags=[]

'''clip at nSigma '''
root.calibrate.photocal.nSigma=3.0

'''Name of the flag field that is set for sources that matched the photometric catalog '''
root.calibrate.photocal.candidateSourceField='calib.photocal.candidate'

'''Name of the flag field that is set for sources used in photometric calibration '''
root.calibrate.photocal.usedSourceField='calib.photocal.used'

'''Apply photometric colour terms (if available) to reference stars '''
root.calibrate.photocal.applyColorTerms=True

'''List of source flag fields that will cause a source to be rejected when they are set. '''
root.calibrate.photocal.badFlags=['flags.pixel.edge', 'flags.pixel.interpolated.any', 'flags.pixel.saturated.any']

'''use median instead of mean to compute zeropoint '''
root.calibrate.photocal.useMedian=True

'''Measure and apply curve of growth? '''
root.calibrate.doCurveOfGrowth=False

import lsst.meas.extensions.photometryKron.version
import lsst.meas.multifit.measureMulti
import lsst.meas.extensions.multiShapelet.multiShapeletLib
import lsst.meas.mosaic.updateExposure
import lsst.meas.extensions.photometryKron.kronLib
import lsst.pipe.tasks.coaddInputRecorder
import lsst.meas.extensions.shapeHSM.version
import lsst.pipe.tasks.selectImages
import lsst.meas.extensions.photometryKron
import lsst.meas.multifit.measureImage
import lsst.meas.extensions.shapeHSM.hsmLib
import lsst.shapelet.shapeletLib
import lsst.pipe.tasks.coaddBase
import lsst.meas.multifit
import lsst.meas.multifit.measureCcd
import lsst.meas.mosaic
import lsst.shapelet
import lsst.meas.multifit.priors
import lsst.meas.mosaic.mosaicLib
import lsst.shapelet.tractor
import lsst.meas.multifit.baseMeasure
import lsst.meas.multifit.measureCoadd
import lsst.meas.multifit.models
import lsst.meas.multifit.version
import lsst.meas.multifit.fitRegion
import lsst.meas.extensions.multiShapelet.version
import lsst.meas.extensions.shapeHSM
import lsst.meas.multifit.multifitLib
import lsst.meas.multifit.samplers
import lsst.meas.multifit.optimizer
import lsst.meas.extensions.multiShapelet
'''Whether to compute quantities related to the Gaussian-weighted shape '''
root.calibrate.measurement.blendedness.doShape=True

'''Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter) '''
root.calibrate.measurement.blendedness.doOld=True

'''Whether to compute quantities related to the Gaussian-weighted flux '''
root.calibrate.measurement.blendedness.doFlux=True

'''Radius factor that sets the maximum extent of the weight function (and hence the flux measurements) '''
root.calibrate.measurement.blendedness.nSigmaWeightMax=3.0

'''Whether to compute blendedness metrics '''
root.calibrate.measurement.doBlendedness=False

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.centroider['centroid.sdss'].priority=0.0

'''if the peak's less thatn this insist on binning at least once '''
root.calibrate.measurement.centroider['centroid.sdss'].peakMin=-1.0

'''fiddle factor for adjusting the binning '''
root.calibrate.measurement.centroider['centroid.sdss'].wfac=1.5

'''maximum allowed binning '''
root.calibrate.measurement.centroider['centroid.sdss'].binmax=16

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.centroider['centroid.naive'].priority=0.0

'''FIXME! NEVER DOCUMENTED! '''
root.calibrate.measurement.centroider['centroid.naive'].background=0.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.centroider['centroid.gaussian'].priority=0.0

root.calibrate.measurement.centroider.name='centroid.sdss'
'''prefix for all measurement fields '''
root.calibrate.measurement.prefix=None

'''Largest aperture for which to use the slow, accurate, sinc aperture code '''
root.calibrate.measurement.algorithms['flux.kron'].maxSincRadius=10.0

'''Number of times to iterate when setting the Kron radius '''
root.calibrate.measurement.algorithms['flux.kron'].nIterForRadius=1

'''Use the Footprint size as part of initial estimate of Kron radius '''
root.calibrate.measurement.algorithms['flux.kron'].useFootprintRadius=False

'''Minimum Kron radius (if == 0.0 use PSF's Kron radius) if enforceMinimumRadius. Also functions as fallback aperture radius if set. '''
root.calibrate.measurement.algorithms['flux.kron'].minimumRadius=0.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['flux.kron'].priority=2.0

'''Multiplier of rms size for aperture used to initially estimate the Kron radius '''
root.calibrate.measurement.algorithms['flux.kron'].nSigmaForRadius=6.0

'''If true check that the Kron radius exceeds some minimum '''
root.calibrate.measurement.algorithms['flux.kron'].enforceMinimumRadius=True

'''if true, use existing shape and centroid measurements instead of fitting '''
root.calibrate.measurement.algorithms['flux.kron'].fixed=False

'''Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K '''
root.calibrate.measurement.algorithms['flux.kron'].smoothingSigma=-1.0

'''Number of Kron radii for Kron flux '''
root.calibrate.measurement.algorithms['flux.kron'].nRadiusForFlux=2.5

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['flux.naive'].priority=2.0

'''FIXME! NEVER DOCUMENTED! '''
root.calibrate.measurement.algorithms['flux.naive'].radius=12.0

'''Shapelet order of inner expansion (0 == Gaussian) '''
root.calibrate.measurement.algorithms['multishapelet.psf'].innerOrder=2

'''Initial radius of inner component in pixels '''
root.calibrate.measurement.algorithms['multishapelet.psf'].initialRadius=1.5

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['multishapelet.psf'].priority=2.0

'''Minimum inner radius in pixels. '''
root.calibrate.measurement.algorithms['multishapelet.psf'].minRadius=0.1

'''outer radius divided by inner radius (fixed) '''
root.calibrate.measurement.algorithms['multishapelet.psf'].radiusRatio=2.0

'''Minimum axis ratio for ellipse (b/a). '''
root.calibrate.measurement.algorithms['multishapelet.psf'].minAxisRatio=0.1

'''outer Gaussian peak height divided by inner Gaussian peak height; held fixed in double-Gaussian ellipse fit, then allowed to vary when shapelets coefficients are fit and ellipses are held fixed. '''
root.calibrate.measurement.algorithms['multishapelet.psf'].peakRatio=0.1

'''Shapelet order of outer expansion (0 == Gaussian) '''
root.calibrate.measurement.algorithms['multishapelet.psf'].outerOrder=1

'''Use fast approximate exponential (good to ~1E-4) '''
root.calibrate.measurement.algorithms['multishapelet.psf'].useApproximateExp=False

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['flux.peakLikelihood'].priority=2.0

'''Name of warping kernel (e.g. "lanczos4") used to compute the peak '''
root.calibrate.measurement.algorithms['flux.peakLikelihood'].warpingKernelName='lanczos4'

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['classification.extendedness'].priority=5.0

'''correction factor for psfFlux error '''
root.calibrate.measurement.algorithms['classification.extendedness'].psfErrFactor=0.0

'''correction factor for modelFlux error '''
root.calibrate.measurement.algorithms['classification.extendedness'].modelErrFactor=0.0

'''critical ratio of model to psf flux '''
root.calibrate.measurement.algorithms['classification.extendedness'].fluxRatio=0.985

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['flags.pixel'].priority=0.0

'''List of mask planes for which to search entire footprint '''
root.calibrate.measurement.algorithms['flags.pixel'].any=[]

'''List of mask planes for which to search center of footprint '''
root.calibrate.measurement.algorithms['flags.pixel'].center=[]

'''Root name of the FitProfileAlgorithm dev comoonent fields. '''
root.calibrate.measurement.algorithms['multishapelet.combo'].devName='multishapelet.dev'

'''Number of pixels to grow the footprint by. '''
root.calibrate.measurement.algorithms['multishapelet.combo'].growFootprint=5

'''Number of half-light radii used to determine the pixels to fit '''
root.calibrate.measurement.algorithms['multishapelet.combo'].radiusInputFactor=4.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['multishapelet.combo'].priority=2.6

'''Mask planes that indicate pixels that should be ignored in the fit. '''
root.calibrate.measurement.algorithms['multishapelet.combo'].badMaskPlanes=['EDGE', 'SAT']

'''Root name of the FitProfileAlgorithm exp component fields. '''
root.calibrate.measurement.algorithms['multishapelet.combo'].expName='multishapelet.exp'

'''If true, individually weigh pixels using the variance image. '''
root.calibrate.measurement.algorithms['multishapelet.combo'].usePixelWeights=False

'''Root name of the FitPsfAlgorithm fields. '''
root.calibrate.measurement.algorithms['multishapelet.combo'].psfName='multishapelet.psf'

'''Use fast approximate exponential (good to ~1E-4) '''
root.calibrate.measurement.algorithms['multishapelet.combo'].useApproximateExp=False

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['shape.hsm.moments'].priority=1.0

'''Mask planes used to reject bad pixels. '''
root.calibrate.measurement.algorithms['shape.hsm.moments'].badMaskPlanes=[]

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['flux.sinc'].priority=2.0

'''major axis of inner boundary (pixels) '''
root.calibrate.measurement.algorithms['flux.sinc'].radius1=0.0

'''major axis of outer boundary (pixels) '''
root.calibrate.measurement.algorithms['flux.sinc'].radius2=12.0

'''measured from x anti-clockwise; radians '''
root.calibrate.measurement.algorithms['flux.sinc'].angle=0.0

'''1 - b/a '''
root.calibrate.measurement.algorithms['flux.sinc'].ellipticity=0.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['jacobian'].priority=3.0

'''Nominal pixel size (arcsec) '''
root.calibrate.measurement.algorithms['jacobian'].pixelScale=0.5

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['shape.hsm.regauss'].priority=1.0

'''Mask planes used to reject bad pixels. '''
root.calibrate.measurement.algorithms['shape.hsm.regauss'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

'''Field name for number of deblend children '''
root.calibrate.measurement.algorithms['shape.hsm.regauss'].deblendNChild=''

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['skycoord'].priority=5.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['flux.psf'].priority=2.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['countInputs'].priority=2.0

'''Name of a registered multi-Gaussian profile. '''
root.calibrate.measurement.algorithms['multishapelet.dev'].profile='tractor-devaucouleur'

'''Minimum half-light radius in units of PSF inner radius for initial parameters. '''
root.calibrate.measurement.algorithms['multishapelet.dev'].minInitialRadius=None

'''Number of pixels to grow the footprint by. '''
root.calibrate.measurement.algorithms['multishapelet.dev'].growFootprint=None

'''Number of half-light radii used to determine the pixels to fit '''
root.calibrate.measurement.algorithms['multishapelet.dev'].radiusInputFactor=None

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['multishapelet.dev'].priority=None

'''Mask planes that indicate pixels that should be ignored in the fit. '''
root.calibrate.measurement.algorithms['multishapelet.dev'].badMaskPlanes=None

'''Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try. '''
root.calibrate.measurement.algorithms['multishapelet.dev'].maxBadPixelFraction=None

'''Attempt to approximately deconvolve the canonical shape before using it to set the initial parameters. '''
root.calibrate.measurement.algorithms['multishapelet.dev'].deconvolveShape=None

'''Minimum axis ratio for ellipse (b/a). '''
root.calibrate.measurement.algorithms['multishapelet.dev'].minAxisRatio=None

'''If true, individually weigh pixels using the variance image. '''
root.calibrate.measurement.algorithms['multishapelet.dev'].usePixelWeights=None

'''Use fast approximate exponential (good to ~1E-4) '''
root.calibrate.measurement.algorithms['multishapelet.dev'].useApproximateExp=None

'''Root name of the FitPsfAlgorithm fields. '''
root.calibrate.measurement.algorithms['multishapelet.dev'].psfName=None

'''Minimum half-light radius in units of PSF inner radius. '''
root.calibrate.measurement.algorithms['multishapelet.dev'].minRadius=None

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['flux.aperture'].priority=2.0

'''Maximum number of radial annuli to measure '''
root.calibrate.measurement.algorithms['flux.aperture'].nApertureMax=10

'''vector of radii for apertures (in pixels) '''
root.calibrate.measurement.algorithms['flux.aperture'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

'''Largest aperture for which to use the slow, accurate, sinc aperture code '''
root.calibrate.measurement.algorithms['flux.aperture'].maxSincRadius=10.0

'''Name of a registered multi-Gaussian profile. '''
root.calibrate.measurement.algorithms['multishapelet.exp'].profile='tractor-exponential'

'''Minimum half-light radius in units of PSF inner radius for initial parameters. '''
root.calibrate.measurement.algorithms['multishapelet.exp'].minInitialRadius=None

'''Number of pixels to grow the footprint by. '''
root.calibrate.measurement.algorithms['multishapelet.exp'].growFootprint=None

'''Number of half-light radii used to determine the pixels to fit '''
root.calibrate.measurement.algorithms['multishapelet.exp'].radiusInputFactor=None

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['multishapelet.exp'].priority=None

'''Mask planes that indicate pixels that should be ignored in the fit. '''
root.calibrate.measurement.algorithms['multishapelet.exp'].badMaskPlanes=None

'''Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try. '''
root.calibrate.measurement.algorithms['multishapelet.exp'].maxBadPixelFraction=None

'''Attempt to approximately deconvolve the canonical shape before using it to set the initial parameters. '''
root.calibrate.measurement.algorithms['multishapelet.exp'].deconvolveShape=None

'''Minimum axis ratio for ellipse (b/a). '''
root.calibrate.measurement.algorithms['multishapelet.exp'].minAxisRatio=None

'''If true, individually weigh pixels using the variance image. '''
root.calibrate.measurement.algorithms['multishapelet.exp'].usePixelWeights=None

'''Use fast approximate exponential (good to ~1E-4) '''
root.calibrate.measurement.algorithms['multishapelet.exp'].useApproximateExp=None

'''Root name of the FitPsfAlgorithm fields. '''
root.calibrate.measurement.algorithms['multishapelet.exp'].psfName=None

'''Minimum half-light radius in units of PSF inner radius. '''
root.calibrate.measurement.algorithms['multishapelet.exp'].minRadius=None

'''suffix of shape field flag to check if fixed is true '''
root.calibrate.measurement.algorithms['flux.gaussian'].shapeFlag='.flags'

'''Convergence tolerance for FWHM '''
root.calibrate.measurement.algorithms['flux.gaussian'].tol2=9.999999747378752e-05

'''Convergence tolerance for e1,e2 '''
root.calibrate.measurement.algorithms['flux.gaussian'].tol1=9.999999747378752e-06

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['flux.gaussian'].priority=2.0

'''name of shape field to use if fixed is true '''
root.calibrate.measurement.algorithms['flux.gaussian'].shape='shape.sdss'

'''name of centroid field to use if fixed is true '''
root.calibrate.measurement.algorithms['flux.gaussian'].centroid='shape.sdss.centroid'

'''FIXME! NEVER DOCUMENTED! '''
root.calibrate.measurement.algorithms['flux.gaussian'].background=0.0

'''Maximum number of iterations '''
root.calibrate.measurement.algorithms['flux.gaussian'].maxIter=100

'''if true, use existing shape and centroid measurements instead of fitting '''
root.calibrate.measurement.algorithms['flux.gaussian'].fixed=True

'''FIXME! NEVER DOCUMENTED! '''
root.calibrate.measurement.algorithms['flux.gaussian'].shiftmax=10.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['shape.hsm.psfMoments'].priority=1.0

'''Convergence tolerance for FWHM '''
root.calibrate.measurement.algorithms['shape.sdss'].tol2=9.999999747378752e-05

'''Convergence tolerance for e1,e2 '''
root.calibrate.measurement.algorithms['shape.sdss'].tol1=9.999999747378752e-06

'''Whether to also compute the shape of the PSF model '''
root.calibrate.measurement.algorithms['shape.sdss'].doMeasurePsf=True

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['shape.sdss'].priority=1.0

'''Additional value to add to background '''
root.calibrate.measurement.algorithms['shape.sdss'].background=0.0

'''Maximum number of iterations '''
root.calibrate.measurement.algorithms['shape.sdss'].maxIter=100

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['flux.scaled'].priority=2.0

'''scaling factor of PSF FWHM for aperture radius '''
root.calibrate.measurement.algorithms['flux.scaled'].scale=3.14

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['centroid.sdss'].priority=0.0

'''if the peak's less thatn this insist on binning at least once '''
root.calibrate.measurement.algorithms['centroid.sdss'].peakMin=-1.0

'''fiddle factor for adjusting the binning '''
root.calibrate.measurement.algorithms['centroid.sdss'].wfac=1.5

'''maximum allowed binning '''
root.calibrate.measurement.algorithms['centroid.sdss'].binmax=16

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['centroid.record'].priority=5.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['shape.hsm.linear'].priority=1.0

'''Mask planes used to reject bad pixels. '''
root.calibrate.measurement.algorithms['shape.hsm.linear'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

'''Field name for number of deblend children '''
root.calibrate.measurement.algorithms['shape.hsm.linear'].deblendNChild=''

'''Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution) '''
root.calibrate.measurement.algorithms['cmodel'].minInitialRadius=0.1

'''Use this multiple of the initial fit ellipse then grow by the PSF width to determine the maximum final fit region size. '''
root.calibrate.measurement.algorithms['cmodel'].region.nFitRadiiMax=3.0

'''Use this multiple of the Kron ellipse to set the fit region (for the final fit region, subject to the nFitRadiiMin and nFitRadiiMax constraints). '''
root.calibrate.measurement.algorithms['cmodel'].region.nKronRadii=1.5

'''Mask planes that indicate pixels that should be ignored in the fit. '''
root.calibrate.measurement.algorithms['cmodel'].region.badMaskPlanes=['EDGE', 'SAT', 'BAD', 'NO_DATA']

'''If the Kron radius is less than this multiple of the PSF width, ignore it and fall back to a PSF-oriented ellipse scaled to match the area of the footprint or this radius (whichever is larger). '''
root.calibrate.measurement.algorithms['cmodel'].region.nPsfSigmaMin=4.0

'''Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try. '''
root.calibrate.measurement.algorithms['cmodel'].region.maxBadPixelFraction=0.1

'''Use this multiple of the initial fit ellipse then grow by the PSF width to determine the minimum final fit region size. '''
root.calibrate.measurement.algorithms['cmodel'].region.nFitRadiiMin=1.0

'''Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region '''
root.calibrate.measurement.algorithms['cmodel'].region.nPsfSigmaGrow=2.0

'''Abort if the fit region grows beyond this many pixels. '''
root.calibrate.measurement.algorithms['cmodel'].region.maxArea=100000

'''If the maximum of the gradient falls below this threshold, consider the algorithm converged '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.gradientThreshold=0.01

'''steps with reduction radio less than this will decrease the trust radius '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.trustRegionShrinkReductionRatio=0.25

'''when increase the trust region size, multiply the radius by this factor '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.trustRegionGrowFactor=2.0

'''steps with length this fraction of the trust radius may increase the trust radius '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.trustRegionGrowStepFraction=0.8

'''whether to save all iterations for debugging purposes '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.doSaveIterations=False

'''steps with reduction radio greater than this may increase the trust radius '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.trustRegionGrowReductionRatio=0.75

'''when reducing the trust region size, multiply the radius by this factor '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.trustRegionShrinkFactor=0.3333333333333333

'''steps with reduction ratio greater than this are accepted '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.stepAcceptThreshold=0.0

'''the initial trust region will be set to this value '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.trustRegionInitialSize=1.0

'''If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.noSR1Term=False

'''relative step size used for numerical derivatives (added to other steps) '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.numDiffRelStep=0.0

'''step size (in units of trust radius) used for numerical derivatives (added to relative step) '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.numDiffTrustRadiusStep=0.1

'''maximum number of steps '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.maxOuterIterations=500

'''Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.skipSR1UpdateThreshold=1e-08

'''maximum number of iterations (i.e. function evaluations and trust region subproblems) per step '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.maxInnerIterations=20

'''If the trust radius falls below this threshold, consider the algorithm converged '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.minTrustRadiusThreshold=0.01

'''value passed as the tolerance to solveTrustRegion '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.trustRegionSolverTolerance=1e-08

'''absolute step size used for numerical derivatives (added to other steps) '''
root.calibrate.measurement.algorithms['cmodel'].initial.optimizer.numDiffAbsStep=0.0

'''Number of degrees of freedom for the Student's T distribution on ln(radius) '''
root.calibrate.measurement.algorithms['cmodel'].initial.empiricalPriorConfig.logRadiusNu=50.0

'''Width of the Student's T distribution in ln(radius). '''
root.calibrate.measurement.algorithms['cmodel'].initial.empiricalPriorConfig.logRadiusSigma=0.45

'''Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T. '''
root.calibrate.measurement.algorithms['cmodel'].initial.empiricalPriorConfig.logRadiusMu=-1.0

'''Width of exponential ellipticity distribution (conformal shear units) '''
root.calibrate.measurement.algorithms['cmodel'].initial.empiricalPriorConfig.ellipticitySigma=0.3

'''Softened core width for ellipticity distribution (conformal shear units '''
root.calibrate.measurement.algorithms['cmodel'].initial.empiricalPriorConfig.ellipticityCore=0.001

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.calibrate.measurement.algorithms['cmodel'].initial.empiricalPriorConfig.logRadiusMinInner=-6.0

'''Minimum ln(radius) '''
root.calibrate.measurement.algorithms['cmodel'].initial.empiricalPriorConfig.logRadiusMinOuter=-6.001

'''One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None '''
root.calibrate.measurement.algorithms['cmodel'].initial.priorSource='EMPIRICAL'

'''Number of Gaussian used to approximate the profile '''
root.calibrate.measurement.algorithms['cmodel'].initial.nComponents=3

'''ln(radius) at which the softened cutoff begins towards the maximum '''
root.calibrate.measurement.algorithms['cmodel'].initial.linearPriorConfig.logRadiusMaxInner=3.0

'''The ratio P(logRadiusMinInner)/P(logRadiusMaxInner) '''
root.calibrate.measurement.algorithms['cmodel'].initial.linearPriorConfig.logRadiusMinMaxRatio=1.0

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.calibrate.measurement.algorithms['cmodel'].initial.linearPriorConfig.logRadiusMinInner=-6.0

'''Maximum ln(radius) '''
root.calibrate.measurement.algorithms['cmodel'].initial.linearPriorConfig.logRadiusMaxOuter=3.001

'''Ellipticity magnitude (conformal shear units) at which the softened cutoff begins '''
root.calibrate.measurement.algorithms['cmodel'].initial.linearPriorConfig.ellipticityMaxInner=2.0

'''Minimum ln(radius) '''
root.calibrate.measurement.algorithms['cmodel'].initial.linearPriorConfig.logRadiusMinOuter=-6.001

'''Maximum ellipticity magnitude (conformal shear units) '''
root.calibrate.measurement.algorithms['cmodel'].initial.linearPriorConfig.ellipticityMaxOuter=2.001

'''Name of the Prior that defines the model to fit (a filename in $MEAS_MULTIFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting. '''
root.calibrate.measurement.algorithms['cmodel'].initial.priorName=''

'''Whether to record the time spent in this stage '''
root.calibrate.measurement.algorithms['cmodel'].initial.doRecordTime=True

'''Name of the shapelet.RadialProfile that defines the model to fit '''
root.calibrate.measurement.algorithms['cmodel'].initial.profileName='lux'

'''Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances) '''
root.calibrate.measurement.algorithms['cmodel'].initial.usePixelWeights=True

'''Maximum radius used in approximating profile with Gaussians (0=default for this profile) '''
root.calibrate.measurement.algorithms['cmodel'].initial.maxRadius=0

'''Whether to record the steps the optimizer takes (or just the number, if running as a plugin) '''
root.calibrate.measurement.algorithms['cmodel'].initial.doRecordHistory=True

'''If the maximum of the gradient falls below this threshold, consider the algorithm converged '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.gradientThreshold=1e-05

'''steps with reduction radio less than this will decrease the trust radius '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.trustRegionShrinkReductionRatio=0.25

'''when increase the trust region size, multiply the radius by this factor '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.trustRegionGrowFactor=2.0

'''steps with length this fraction of the trust radius may increase the trust radius '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.trustRegionGrowStepFraction=0.8

'''whether to save all iterations for debugging purposes '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.doSaveIterations=False

'''steps with reduction radio greater than this may increase the trust radius '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.trustRegionGrowReductionRatio=0.75

'''when reducing the trust region size, multiply the radius by this factor '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.trustRegionShrinkFactor=0.3333333333333333

'''steps with reduction ratio greater than this are accepted '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.stepAcceptThreshold=0.0

'''the initial trust region will be set to this value '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.trustRegionInitialSize=1.0

'''If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.noSR1Term=False

'''relative step size used for numerical derivatives (added to other steps) '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.numDiffRelStep=0.0

'''step size (in units of trust radius) used for numerical derivatives (added to relative step) '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.numDiffTrustRadiusStep=0.1

'''maximum number of steps '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.maxOuterIterations=500

'''Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.skipSR1UpdateThreshold=1e-08

'''maximum number of iterations (i.e. function evaluations and trust region subproblems) per step '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.maxInnerIterations=20

'''If the trust radius falls below this threshold, consider the algorithm converged '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.minTrustRadiusThreshold=1e-05

'''value passed as the tolerance to solveTrustRegion '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.trustRegionSolverTolerance=1e-08

'''absolute step size used for numerical derivatives (added to other steps) '''
root.calibrate.measurement.algorithms['cmodel'].dev.optimizer.numDiffAbsStep=0.0

'''Number of degrees of freedom for the Student's T distribution on ln(radius) '''
root.calibrate.measurement.algorithms['cmodel'].dev.empiricalPriorConfig.logRadiusNu=50.0

'''Width of the Student's T distribution in ln(radius). '''
root.calibrate.measurement.algorithms['cmodel'].dev.empiricalPriorConfig.logRadiusSigma=0.45

'''Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T. '''
root.calibrate.measurement.algorithms['cmodel'].dev.empiricalPriorConfig.logRadiusMu=-1.0

'''Width of exponential ellipticity distribution (conformal shear units) '''
root.calibrate.measurement.algorithms['cmodel'].dev.empiricalPriorConfig.ellipticitySigma=0.3

'''Softened core width for ellipticity distribution (conformal shear units '''
root.calibrate.measurement.algorithms['cmodel'].dev.empiricalPriorConfig.ellipticityCore=0.001

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.calibrate.measurement.algorithms['cmodel'].dev.empiricalPriorConfig.logRadiusMinInner=-6.0

'''Minimum ln(radius) '''
root.calibrate.measurement.algorithms['cmodel'].dev.empiricalPriorConfig.logRadiusMinOuter=-6.001

'''One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None '''
root.calibrate.measurement.algorithms['cmodel'].dev.priorSource='EMPIRICAL'

'''Number of Gaussian used to approximate the profile '''
root.calibrate.measurement.algorithms['cmodel'].dev.nComponents=8

'''ln(radius) at which the softened cutoff begins towards the maximum '''
root.calibrate.measurement.algorithms['cmodel'].dev.linearPriorConfig.logRadiusMaxInner=3.0

'''The ratio P(logRadiusMinInner)/P(logRadiusMaxInner) '''
root.calibrate.measurement.algorithms['cmodel'].dev.linearPriorConfig.logRadiusMinMaxRatio=1.0

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.calibrate.measurement.algorithms['cmodel'].dev.linearPriorConfig.logRadiusMinInner=-6.0

'''Maximum ln(radius) '''
root.calibrate.measurement.algorithms['cmodel'].dev.linearPriorConfig.logRadiusMaxOuter=3.001

'''Ellipticity magnitude (conformal shear units) at which the softened cutoff begins '''
root.calibrate.measurement.algorithms['cmodel'].dev.linearPriorConfig.ellipticityMaxInner=2.0

'''Minimum ln(radius) '''
root.calibrate.measurement.algorithms['cmodel'].dev.linearPriorConfig.logRadiusMinOuter=-6.001

'''Maximum ellipticity magnitude (conformal shear units) '''
root.calibrate.measurement.algorithms['cmodel'].dev.linearPriorConfig.ellipticityMaxOuter=2.001

'''Name of the Prior that defines the model to fit (a filename in $MEAS_MULTIFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting. '''
root.calibrate.measurement.algorithms['cmodel'].dev.priorName=''

'''Whether to record the time spent in this stage '''
root.calibrate.measurement.algorithms['cmodel'].dev.doRecordTime=True

'''Name of the shapelet.RadialProfile that defines the model to fit '''
root.calibrate.measurement.algorithms['cmodel'].dev.profileName='luv'

'''Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances) '''
root.calibrate.measurement.algorithms['cmodel'].dev.usePixelWeights=False

'''Maximum radius used in approximating profile with Gaussians (0=default for this profile) '''
root.calibrate.measurement.algorithms['cmodel'].dev.maxRadius=0

'''Whether to record the steps the optimizer takes (or just the number, if running as a plugin) '''
root.calibrate.measurement.algorithms['cmodel'].dev.doRecordHistory=True

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['cmodel'].priority=2.5

'''If the maximum of the gradient falls below this threshold, consider the algorithm converged '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.gradientThreshold=1e-05

'''steps with reduction radio less than this will decrease the trust radius '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.trustRegionShrinkReductionRatio=0.25

'''when increase the trust region size, multiply the radius by this factor '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.trustRegionGrowFactor=2.0

'''steps with length this fraction of the trust radius may increase the trust radius '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.trustRegionGrowStepFraction=0.8

'''whether to save all iterations for debugging purposes '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.doSaveIterations=False

'''steps with reduction radio greater than this may increase the trust radius '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.trustRegionGrowReductionRatio=0.75

'''when reducing the trust region size, multiply the radius by this factor '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.trustRegionShrinkFactor=0.3333333333333333

'''steps with reduction ratio greater than this are accepted '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.stepAcceptThreshold=0.0

'''the initial trust region will be set to this value '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.trustRegionInitialSize=1.0

'''If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.noSR1Term=False

'''relative step size used for numerical derivatives (added to other steps) '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.numDiffRelStep=0.0

'''step size (in units of trust radius) used for numerical derivatives (added to relative step) '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.numDiffTrustRadiusStep=0.1

'''maximum number of steps '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.maxOuterIterations=250

'''Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.skipSR1UpdateThreshold=1e-08

'''maximum number of iterations (i.e. function evaluations and trust region subproblems) per step '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.maxInnerIterations=20

'''If the trust radius falls below this threshold, consider the algorithm converged '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.minTrustRadiusThreshold=1e-05

'''value passed as the tolerance to solveTrustRegion '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.trustRegionSolverTolerance=1e-08

'''absolute step size used for numerical derivatives (added to other steps) '''
root.calibrate.measurement.algorithms['cmodel'].exp.optimizer.numDiffAbsStep=0.0

'''Number of degrees of freedom for the Student's T distribution on ln(radius) '''
root.calibrate.measurement.algorithms['cmodel'].exp.empiricalPriorConfig.logRadiusNu=50.0

'''Width of the Student's T distribution in ln(radius). '''
root.calibrate.measurement.algorithms['cmodel'].exp.empiricalPriorConfig.logRadiusSigma=0.45

'''Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T. '''
root.calibrate.measurement.algorithms['cmodel'].exp.empiricalPriorConfig.logRadiusMu=-1.0

'''Width of exponential ellipticity distribution (conformal shear units) '''
root.calibrate.measurement.algorithms['cmodel'].exp.empiricalPriorConfig.ellipticitySigma=0.3

'''Softened core width for ellipticity distribution (conformal shear units '''
root.calibrate.measurement.algorithms['cmodel'].exp.empiricalPriorConfig.ellipticityCore=0.001

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.calibrate.measurement.algorithms['cmodel'].exp.empiricalPriorConfig.logRadiusMinInner=-6.0

'''Minimum ln(radius) '''
root.calibrate.measurement.algorithms['cmodel'].exp.empiricalPriorConfig.logRadiusMinOuter=-6.001

'''One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None '''
root.calibrate.measurement.algorithms['cmodel'].exp.priorSource='EMPIRICAL'

'''Number of Gaussian used to approximate the profile '''
root.calibrate.measurement.algorithms['cmodel'].exp.nComponents=6

'''ln(radius) at which the softened cutoff begins towards the maximum '''
root.calibrate.measurement.algorithms['cmodel'].exp.linearPriorConfig.logRadiusMaxInner=3.0

'''The ratio P(logRadiusMinInner)/P(logRadiusMaxInner) '''
root.calibrate.measurement.algorithms['cmodel'].exp.linearPriorConfig.logRadiusMinMaxRatio=1.0

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.calibrate.measurement.algorithms['cmodel'].exp.linearPriorConfig.logRadiusMinInner=-6.0

'''Maximum ln(radius) '''
root.calibrate.measurement.algorithms['cmodel'].exp.linearPriorConfig.logRadiusMaxOuter=3.001

'''Ellipticity magnitude (conformal shear units) at which the softened cutoff begins '''
root.calibrate.measurement.algorithms['cmodel'].exp.linearPriorConfig.ellipticityMaxInner=2.0

'''Minimum ln(radius) '''
root.calibrate.measurement.algorithms['cmodel'].exp.linearPriorConfig.logRadiusMinOuter=-6.001

'''Maximum ellipticity magnitude (conformal shear units) '''
root.calibrate.measurement.algorithms['cmodel'].exp.linearPriorConfig.ellipticityMaxOuter=2.001

'''Name of the Prior that defines the model to fit (a filename in $MEAS_MULTIFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting. '''
root.calibrate.measurement.algorithms['cmodel'].exp.priorName=''

'''Whether to record the time spent in this stage '''
root.calibrate.measurement.algorithms['cmodel'].exp.doRecordTime=True

'''Name of the shapelet.RadialProfile that defines the model to fit '''
root.calibrate.measurement.algorithms['cmodel'].exp.profileName='lux'

'''Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances) '''
root.calibrate.measurement.algorithms['cmodel'].exp.usePixelWeights=False

'''Maximum radius used in approximating profile with Gaussians (0=default for this profile) '''
root.calibrate.measurement.algorithms['cmodel'].exp.maxRadius=0

'''Whether to record the steps the optimizer takes (or just the number, if running as a plugin) '''
root.calibrate.measurement.algorithms['cmodel'].exp.doRecordHistory=True

'''If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this.  If <= 0.0, abort the fit early instead. '''
root.calibrate.measurement.algorithms['cmodel'].fallbackInitialMomentsPsfFactor=1.5

'''Root name of the FitPsfAlgorithm fields. '''
root.calibrate.measurement.algorithms['cmodel'].psfName='multishapelet.psf'

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['focalplane'].priority=3.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['centroid.naive'].priority=0.0

'''FIXME! NEVER DOCUMENTED! '''
root.calibrate.measurement.algorithms['centroid.naive'].background=0.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['shape.hsm.bj'].priority=1.0

'''Mask planes used to reject bad pixels. '''
root.calibrate.measurement.algorithms['shape.hsm.bj'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

'''Field name for number of deblend children '''
root.calibrate.measurement.algorithms['shape.hsm.bj'].deblendNChild=''

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['centroid.gaussian'].priority=0.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['flux.aperture.elliptical'].priority=1.899999976158142

'''Maximum number of radial annuli to measure '''
root.calibrate.measurement.algorithms['flux.aperture.elliptical'].nApertureMax=10

'''vector of radii for apertures (in pixels) '''
root.calibrate.measurement.algorithms['flux.aperture.elliptical'].radii=[1.0, 1.5625, 2.44140625, 3.814697265625, 5.9604644775390625, 9.313225746154785, 14.551915228366852, 22.737367544323206, 35.52713678800501, 55.51115123125783]

'''Largest aperture for which to use the slow, accurate, sinc aperture code '''
root.calibrate.measurement.algorithms['flux.aperture.elliptical'].maxSincRadius=10.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['variance'].priority=2.0

'''scale factor to apply to shape for aperture '''
root.calibrate.measurement.algorithms['variance'].scale=5.0

'''mask planes to ignore '''
root.calibrate.measurement.algorithms['variance'].mask=['DETECTED', 'DETECTED_NEGATIVE']

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['correctfluxes'].priority=3.0

'''List of flux fields that should not be corrected (otherwise all fields in getApCorrRegistry() will be) '''
root.calibrate.measurement.algorithms['correctfluxes'].ignored=[]

'''Whether to propagate aperture correction uncertainties into flux uncertainties '''
root.calibrate.measurement.algorithms['correctfluxes'].doPropagateErrors=False

'''Whether to set the general failure flag for a flux when it cannot be aperture-corrected '''
root.calibrate.measurement.algorithms['correctfluxes'].doFlagApCorrFailures=True

'''Whether to save the per-source per-flux aperture corrections and their errors '''
root.calibrate.measurement.algorithms['correctfluxes'].doRecordApCorr=True

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.calibrate.measurement.algorithms['shape.hsm.ksb'].priority=1.0

'''Mask planes used to reject bad pixels. '''
root.calibrate.measurement.algorithms['shape.hsm.ksb'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

'''Field name for number of deblend children '''
root.calibrate.measurement.algorithms['shape.hsm.ksb'].deblendNChild=''

root.calibrate.measurement.algorithms.names=['flux.psf', 'flags.pixel', 'flux.gaussian', 'flux.aperture', 'flux.naive', 'cmodel', 'flux.kron', 'shape.hsm.moments', 'centroid.naive', 'flux.sinc', 'shape.sdss', 'shape.hsm.regauss', 'multishapelet.psf', 'correctfluxes', 'classification.extendedness', 'skycoord', 'shape.hsm.psfMoments']
'''The seed value to use for random number generation. '''
root.calibrate.measurement.replaceWithNoise.noiseSeed=0

'''Add ann offset to the generated noise. '''
root.calibrate.measurement.replaceWithNoise.noiseOffset=0.0

'''How do we choose the mean and variance of the Gaussian noise we generate?
Allowed values:
	variance	Mean = 0, variance = the image's variance
	meta	Mean = 0, variance = the "BGMEAN" metadata entry
	measure	Measure clipped mean and variance from the whole image
 '''
root.calibrate.measurement.replaceWithNoise.noiseSource='measure'

'''the name of the aperture flux algorithm used for calibration '''
root.calibrate.measurement.slots.calibFlux='flux.naive'

'''the name of the algorithm used to set the source aperture flux slot '''
root.calibrate.measurement.slots.apFlux='flux.sinc'

'''the name of the algorithm used to set the source inst flux slot '''
root.calibrate.measurement.slots.instFlux='flux.gaussian'

'''the name of the algorithm used to set source moments parameters '''
root.calibrate.measurement.slots.shape='shape.hsm.moments'

'''the name of the centroiding algorithm used to set source x,y '''
root.calibrate.measurement.slots.centroid='centroid.sdss'

'''the name of the algorithm used to set the source model flux slot '''
root.calibrate.measurement.slots.modelFlux='flux.gaussian'
root.calibrate.measurement.algorithms.names.remove('cmodel')

'''the name of the algorithm used to set the source psf flux slot '''
root.calibrate.measurement.slots.psfFlux='flux.psf'

'''When measuring, replace other detected footprints with noise? '''
root.calibrate.measurement.doReplaceWithNoise=True

'''Require astrometry to succeed, if activated? '''
root.calibrate.requireAstrometry=True

'''Whether to compute quantities related to the Gaussian-weighted shape '''
root.measurement.blendedness.doShape=True

'''Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter) '''
root.measurement.blendedness.doOld=True

'''Whether to compute quantities related to the Gaussian-weighted flux '''
root.measurement.blendedness.doFlux=True

'''Radius factor that sets the maximum extent of the weight function (and hence the flux measurements) '''
root.measurement.blendedness.nSigmaWeightMax=3.0

'''Whether to compute blendedness metrics '''
root.measurement.doBlendedness=False

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.centroider['centroid.sdss'].priority=0.0

'''if the peak's less thatn this insist on binning at least once '''
root.measurement.centroider['centroid.sdss'].peakMin=-1.0

'''fiddle factor for adjusting the binning '''
root.measurement.centroider['centroid.sdss'].wfac=1.5

'''maximum allowed binning '''
root.measurement.centroider['centroid.sdss'].binmax=16

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.centroider['centroid.naive'].priority=0.0

'''FIXME! NEVER DOCUMENTED! '''
root.measurement.centroider['centroid.naive'].background=0.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.centroider['centroid.gaussian'].priority=0.0

root.measurement.centroider.name='centroid.sdss'
'''prefix for all measurement fields '''
root.measurement.prefix=None

'''Largest aperture for which to use the slow, accurate, sinc aperture code '''
root.measurement.algorithms['flux.kron'].maxSincRadius=10.0

'''Number of times to iterate when setting the Kron radius '''
root.measurement.algorithms['flux.kron'].nIterForRadius=1

'''Use the Footprint size as part of initial estimate of Kron radius '''
root.measurement.algorithms['flux.kron'].useFootprintRadius=False

'''Minimum Kron radius (if == 0.0 use PSF's Kron radius) if enforceMinimumRadius. Also functions as fallback aperture radius if set. '''
root.measurement.algorithms['flux.kron'].minimumRadius=0.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['flux.kron'].priority=2.0

'''Multiplier of rms size for aperture used to initially estimate the Kron radius '''
root.measurement.algorithms['flux.kron'].nSigmaForRadius=6.0

'''If true check that the Kron radius exceeds some minimum '''
root.measurement.algorithms['flux.kron'].enforceMinimumRadius=True

'''if true, use existing shape and centroid measurements instead of fitting '''
root.measurement.algorithms['flux.kron'].fixed=False

'''Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K '''
root.measurement.algorithms['flux.kron'].smoothingSigma=-1.0

'''Number of Kron radii for Kron flux '''
root.measurement.algorithms['flux.kron'].nRadiusForFlux=2.5

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['flux.naive'].priority=2.0

'''FIXME! NEVER DOCUMENTED! '''
root.measurement.algorithms['flux.naive'].radius=12.0

'''Shapelet order of inner expansion (0 == Gaussian) '''
root.measurement.algorithms['multishapelet.psf'].innerOrder=2

'''Initial radius of inner component in pixels '''
root.measurement.algorithms['multishapelet.psf'].initialRadius=1.5

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['multishapelet.psf'].priority=2.0

'''Minimum inner radius in pixels. '''
root.measurement.algorithms['multishapelet.psf'].minRadius=0.1

'''outer radius divided by inner radius (fixed) '''
root.measurement.algorithms['multishapelet.psf'].radiusRatio=2.0

'''Minimum axis ratio for ellipse (b/a). '''
root.measurement.algorithms['multishapelet.psf'].minAxisRatio=0.1

'''outer Gaussian peak height divided by inner Gaussian peak height; held fixed in double-Gaussian ellipse fit, then allowed to vary when shapelets coefficients are fit and ellipses are held fixed. '''
root.measurement.algorithms['multishapelet.psf'].peakRatio=0.1

'''Shapelet order of outer expansion (0 == Gaussian) '''
root.measurement.algorithms['multishapelet.psf'].outerOrder=1

'''Use fast approximate exponential (good to ~1E-4) '''
root.measurement.algorithms['multishapelet.psf'].useApproximateExp=False

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['flux.peakLikelihood'].priority=2.0

'''Name of warping kernel (e.g. "lanczos4") used to compute the peak '''
root.measurement.algorithms['flux.peakLikelihood'].warpingKernelName='lanczos4'

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['classification.extendedness'].priority=5.0

'''correction factor for psfFlux error '''
root.measurement.algorithms['classification.extendedness'].psfErrFactor=0.0

'''correction factor for modelFlux error '''
root.measurement.algorithms['classification.extendedness'].modelErrFactor=0.0

'''critical ratio of model to psf flux '''
root.measurement.algorithms['classification.extendedness'].fluxRatio=0.95

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['flags.pixel'].priority=0.0

'''List of mask planes for which to search entire footprint '''
root.measurement.algorithms['flags.pixel'].any=[]

'''List of mask planes for which to search center of footprint '''
root.measurement.algorithms['flags.pixel'].center=[]

'''Root name of the FitProfileAlgorithm dev comoonent fields. '''
root.measurement.algorithms['multishapelet.combo'].devName='multishapelet.dev'

'''Number of pixels to grow the footprint by. '''
root.measurement.algorithms['multishapelet.combo'].growFootprint=5

'''Number of half-light radii used to determine the pixels to fit '''
root.measurement.algorithms['multishapelet.combo'].radiusInputFactor=4.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['multishapelet.combo'].priority=2.6

'''Mask planes that indicate pixels that should be ignored in the fit. '''
root.measurement.algorithms['multishapelet.combo'].badMaskPlanes=['EDGE', 'SAT']

'''Root name of the FitProfileAlgorithm exp component fields. '''
root.measurement.algorithms['multishapelet.combo'].expName='multishapelet.exp'

'''If true, individually weigh pixels using the variance image. '''
root.measurement.algorithms['multishapelet.combo'].usePixelWeights=False

'''Root name of the FitPsfAlgorithm fields. '''
root.measurement.algorithms['multishapelet.combo'].psfName='multishapelet.psf'

'''Use fast approximate exponential (good to ~1E-4) '''
root.measurement.algorithms['multishapelet.combo'].useApproximateExp=False

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['shape.hsm.moments'].priority=1.0

'''Mask planes used to reject bad pixels. '''
root.measurement.algorithms['shape.hsm.moments'].badMaskPlanes=[]

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['flux.sinc'].priority=2.0

'''major axis of inner boundary (pixels) '''
root.measurement.algorithms['flux.sinc'].radius1=0.0

'''major axis of outer boundary (pixels) '''
root.measurement.algorithms['flux.sinc'].radius2=12.0

'''measured from x anti-clockwise; radians '''
root.measurement.algorithms['flux.sinc'].angle=0.0

'''1 - b/a '''
root.measurement.algorithms['flux.sinc'].ellipticity=0.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['jacobian'].priority=3.0

'''Nominal pixel size (arcsec) '''
root.measurement.algorithms['jacobian'].pixelScale=0.168

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['shape.hsm.regauss'].priority=1.0

'''Mask planes used to reject bad pixels. '''
root.measurement.algorithms['shape.hsm.regauss'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

'''Field name for number of deblend children '''
root.measurement.algorithms['shape.hsm.regauss'].deblendNChild='deblend.nchild'

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['skycoord'].priority=5.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['flux.psf'].priority=2.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['countInputs'].priority=2.0

'''Name of a registered multi-Gaussian profile. '''
root.measurement.algorithms['multishapelet.dev'].profile='tractor-devaucouleur'

'''Minimum half-light radius in units of PSF inner radius for initial parameters. '''
root.measurement.algorithms['multishapelet.dev'].minInitialRadius=None

'''Number of pixels to grow the footprint by. '''
root.measurement.algorithms['multishapelet.dev'].growFootprint=None

'''Number of half-light radii used to determine the pixels to fit '''
root.measurement.algorithms['multishapelet.dev'].radiusInputFactor=None

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['multishapelet.dev'].priority=None

'''Mask planes that indicate pixels that should be ignored in the fit. '''
root.measurement.algorithms['multishapelet.dev'].badMaskPlanes=None

'''Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try. '''
root.measurement.algorithms['multishapelet.dev'].maxBadPixelFraction=None

'''Attempt to approximately deconvolve the canonical shape before using it to set the initial parameters. '''
root.measurement.algorithms['multishapelet.dev'].deconvolveShape=None

'''Minimum axis ratio for ellipse (b/a). '''
root.measurement.algorithms['multishapelet.dev'].minAxisRatio=None

'''If true, individually weigh pixels using the variance image. '''
root.measurement.algorithms['multishapelet.dev'].usePixelWeights=None

'''Use fast approximate exponential (good to ~1E-4) '''
root.measurement.algorithms['multishapelet.dev'].useApproximateExp=None

'''Root name of the FitPsfAlgorithm fields. '''
root.measurement.algorithms['multishapelet.dev'].psfName=None

'''Minimum half-light radius in units of PSF inner radius. '''
root.measurement.algorithms['multishapelet.dev'].minRadius=None

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['flux.aperture'].priority=2.0

'''Maximum number of radial annuli to measure '''
root.measurement.algorithms['flux.aperture'].nApertureMax=10

'''vector of radii for apertures (in pixels) '''
root.measurement.algorithms['flux.aperture'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

'''Largest aperture for which to use the slow, accurate, sinc aperture code '''
root.measurement.algorithms['flux.aperture'].maxSincRadius=10.0

'''Name of a registered multi-Gaussian profile. '''
root.measurement.algorithms['multishapelet.exp'].profile='tractor-exponential'

'''Minimum half-light radius in units of PSF inner radius for initial parameters. '''
root.measurement.algorithms['multishapelet.exp'].minInitialRadius=None

'''Number of pixels to grow the footprint by. '''
root.measurement.algorithms['multishapelet.exp'].growFootprint=None

'''Number of half-light radii used to determine the pixels to fit '''
root.measurement.algorithms['multishapelet.exp'].radiusInputFactor=None

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['multishapelet.exp'].priority=None

'''Mask planes that indicate pixels that should be ignored in the fit. '''
root.measurement.algorithms['multishapelet.exp'].badMaskPlanes=None

'''Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try. '''
root.measurement.algorithms['multishapelet.exp'].maxBadPixelFraction=None

'''Attempt to approximately deconvolve the canonical shape before using it to set the initial parameters. '''
root.measurement.algorithms['multishapelet.exp'].deconvolveShape=None

'''Minimum axis ratio for ellipse (b/a). '''
root.measurement.algorithms['multishapelet.exp'].minAxisRatio=None

'''If true, individually weigh pixels using the variance image. '''
root.measurement.algorithms['multishapelet.exp'].usePixelWeights=None

'''Use fast approximate exponential (good to ~1E-4) '''
root.measurement.algorithms['multishapelet.exp'].useApproximateExp=None

'''Root name of the FitPsfAlgorithm fields. '''
root.measurement.algorithms['multishapelet.exp'].psfName=None

'''Minimum half-light radius in units of PSF inner radius. '''
root.measurement.algorithms['multishapelet.exp'].minRadius=None

'''suffix of shape field flag to check if fixed is true '''
root.measurement.algorithms['flux.gaussian'].shapeFlag='.flags'

'''Convergence tolerance for FWHM '''
root.measurement.algorithms['flux.gaussian'].tol2=9.999999747378752e-05

'''Convergence tolerance for e1,e2 '''
root.measurement.algorithms['flux.gaussian'].tol1=9.999999747378752e-06

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['flux.gaussian'].priority=2.0

'''name of shape field to use if fixed is true '''
root.measurement.algorithms['flux.gaussian'].shape='shape.sdss'

'''name of centroid field to use if fixed is true '''
root.measurement.algorithms['flux.gaussian'].centroid='shape.sdss.centroid'

'''FIXME! NEVER DOCUMENTED! '''
root.measurement.algorithms['flux.gaussian'].background=0.0

'''Maximum number of iterations '''
root.measurement.algorithms['flux.gaussian'].maxIter=100

'''if true, use existing shape and centroid measurements instead of fitting '''
root.measurement.algorithms['flux.gaussian'].fixed=True

'''FIXME! NEVER DOCUMENTED! '''
root.measurement.algorithms['flux.gaussian'].shiftmax=10.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['shape.hsm.psfMoments'].priority=1.0

'''Convergence tolerance for FWHM '''
root.measurement.algorithms['shape.sdss'].tol2=9.999999747378752e-05

'''Convergence tolerance for e1,e2 '''
root.measurement.algorithms['shape.sdss'].tol1=9.999999747378752e-06

'''Whether to also compute the shape of the PSF model '''
root.measurement.algorithms['shape.sdss'].doMeasurePsf=True

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['shape.sdss'].priority=1.0

'''Additional value to add to background '''
root.measurement.algorithms['shape.sdss'].background=0.0

'''Maximum number of iterations '''
root.measurement.algorithms['shape.sdss'].maxIter=100

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['flux.scaled'].priority=2.0

'''scaling factor of PSF FWHM for aperture radius '''
root.measurement.algorithms['flux.scaled'].scale=3.14

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['centroid.sdss'].priority=0.0

'''if the peak's less thatn this insist on binning at least once '''
root.measurement.algorithms['centroid.sdss'].peakMin=-1.0

'''fiddle factor for adjusting the binning '''
root.measurement.algorithms['centroid.sdss'].wfac=1.5

'''maximum allowed binning '''
root.measurement.algorithms['centroid.sdss'].binmax=16

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['centroid.record'].priority=5.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['shape.hsm.linear'].priority=1.0

'''Mask planes used to reject bad pixels. '''
root.measurement.algorithms['shape.hsm.linear'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

'''Field name for number of deblend children '''
root.measurement.algorithms['shape.hsm.linear'].deblendNChild=''

'''Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution) '''
root.measurement.algorithms['cmodel'].minInitialRadius=0.1

'''Use this multiple of the initial fit ellipse then grow by the PSF width to determine the maximum final fit region size. '''
root.measurement.algorithms['cmodel'].region.nFitRadiiMax=3.0

'''Use this multiple of the Kron ellipse to set the fit region (for the final fit region, subject to the nFitRadiiMin and nFitRadiiMax constraints). '''
root.measurement.algorithms['cmodel'].region.nKronRadii=1.5

'''Mask planes that indicate pixels that should be ignored in the fit. '''
root.measurement.algorithms['cmodel'].region.badMaskPlanes=['EDGE', 'SAT', 'BAD', 'NO_DATA']

'''If the Kron radius is less than this multiple of the PSF width, ignore it and fall back to a PSF-oriented ellipse scaled to match the area of the footprint or this radius (whichever is larger). '''
root.measurement.algorithms['cmodel'].region.nPsfSigmaMin=4.0

'''Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try. '''
root.measurement.algorithms['cmodel'].region.maxBadPixelFraction=0.1

'''Use this multiple of the initial fit ellipse then grow by the PSF width to determine the minimum final fit region size. '''
root.measurement.algorithms['cmodel'].region.nFitRadiiMin=1.0

'''Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region '''
root.measurement.algorithms['cmodel'].region.nPsfSigmaGrow=2.0

'''Abort if the fit region grows beyond this many pixels. '''
root.measurement.algorithms['cmodel'].region.maxArea=100000

'''If the maximum of the gradient falls below this threshold, consider the algorithm converged '''
root.measurement.algorithms['cmodel'].initial.optimizer.gradientThreshold=0.01

'''steps with reduction radio less than this will decrease the trust radius '''
root.measurement.algorithms['cmodel'].initial.optimizer.trustRegionShrinkReductionRatio=0.25

'''when increase the trust region size, multiply the radius by this factor '''
root.measurement.algorithms['cmodel'].initial.optimizer.trustRegionGrowFactor=2.0

'''steps with length this fraction of the trust radius may increase the trust radius '''
root.measurement.algorithms['cmodel'].initial.optimizer.trustRegionGrowStepFraction=0.8

'''whether to save all iterations for debugging purposes '''
root.measurement.algorithms['cmodel'].initial.optimizer.doSaveIterations=False

'''steps with reduction radio greater than this may increase the trust radius '''
root.measurement.algorithms['cmodel'].initial.optimizer.trustRegionGrowReductionRatio=0.75

'''when reducing the trust region size, multiply the radius by this factor '''
root.measurement.algorithms['cmodel'].initial.optimizer.trustRegionShrinkFactor=0.3333333333333333

'''steps with reduction ratio greater than this are accepted '''
root.measurement.algorithms['cmodel'].initial.optimizer.stepAcceptThreshold=0.0

'''the initial trust region will be set to this value '''
root.measurement.algorithms['cmodel'].initial.optimizer.trustRegionInitialSize=1.0

'''If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method '''
root.measurement.algorithms['cmodel'].initial.optimizer.noSR1Term=False

'''relative step size used for numerical derivatives (added to other steps) '''
root.measurement.algorithms['cmodel'].initial.optimizer.numDiffRelStep=0.0

'''step size (in units of trust radius) used for numerical derivatives (added to relative step) '''
root.measurement.algorithms['cmodel'].initial.optimizer.numDiffTrustRadiusStep=0.1

'''maximum number of steps '''
root.measurement.algorithms['cmodel'].initial.optimizer.maxOuterIterations=500

'''Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold '''
root.measurement.algorithms['cmodel'].initial.optimizer.skipSR1UpdateThreshold=1e-08

'''maximum number of iterations (i.e. function evaluations and trust region subproblems) per step '''
root.measurement.algorithms['cmodel'].initial.optimizer.maxInnerIterations=20

'''If the trust radius falls below this threshold, consider the algorithm converged '''
root.measurement.algorithms['cmodel'].initial.optimizer.minTrustRadiusThreshold=0.01

'''value passed as the tolerance to solveTrustRegion '''
root.measurement.algorithms['cmodel'].initial.optimizer.trustRegionSolverTolerance=1e-08

'''absolute step size used for numerical derivatives (added to other steps) '''
root.measurement.algorithms['cmodel'].initial.optimizer.numDiffAbsStep=0.0

'''Number of degrees of freedom for the Student's T distribution on ln(radius) '''
root.measurement.algorithms['cmodel'].initial.empiricalPriorConfig.logRadiusNu=50.0

'''Width of the Student's T distribution in ln(radius). '''
root.measurement.algorithms['cmodel'].initial.empiricalPriorConfig.logRadiusSigma=0.45

'''Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T. '''
root.measurement.algorithms['cmodel'].initial.empiricalPriorConfig.logRadiusMu=-1.0

'''Width of exponential ellipticity distribution (conformal shear units) '''
root.measurement.algorithms['cmodel'].initial.empiricalPriorConfig.ellipticitySigma=0.3

'''Softened core width for ellipticity distribution (conformal shear units '''
root.measurement.algorithms['cmodel'].initial.empiricalPriorConfig.ellipticityCore=0.001

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.measurement.algorithms['cmodel'].initial.empiricalPriorConfig.logRadiusMinInner=-6.0

'''Minimum ln(radius) '''
root.measurement.algorithms['cmodel'].initial.empiricalPriorConfig.logRadiusMinOuter=-6.001

'''One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None '''
root.measurement.algorithms['cmodel'].initial.priorSource='EMPIRICAL'

'''Number of Gaussian used to approximate the profile '''
root.measurement.algorithms['cmodel'].initial.nComponents=3

'''ln(radius) at which the softened cutoff begins towards the maximum '''
root.measurement.algorithms['cmodel'].initial.linearPriorConfig.logRadiusMaxInner=3.0

'''The ratio P(logRadiusMinInner)/P(logRadiusMaxInner) '''
root.measurement.algorithms['cmodel'].initial.linearPriorConfig.logRadiusMinMaxRatio=1.0

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.measurement.algorithms['cmodel'].initial.linearPriorConfig.logRadiusMinInner=-6.0

'''Maximum ln(radius) '''
root.measurement.algorithms['cmodel'].initial.linearPriorConfig.logRadiusMaxOuter=3.001

'''Ellipticity magnitude (conformal shear units) at which the softened cutoff begins '''
root.measurement.algorithms['cmodel'].initial.linearPriorConfig.ellipticityMaxInner=2.0

'''Minimum ln(radius) '''
root.measurement.algorithms['cmodel'].initial.linearPriorConfig.logRadiusMinOuter=-6.001

'''Maximum ellipticity magnitude (conformal shear units) '''
root.measurement.algorithms['cmodel'].initial.linearPriorConfig.ellipticityMaxOuter=2.001

'''Name of the Prior that defines the model to fit (a filename in $MEAS_MULTIFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting. '''
root.measurement.algorithms['cmodel'].initial.priorName=''

'''Whether to record the time spent in this stage '''
root.measurement.algorithms['cmodel'].initial.doRecordTime=True

'''Name of the shapelet.RadialProfile that defines the model to fit '''
root.measurement.algorithms['cmodel'].initial.profileName='lux'

'''Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances) '''
root.measurement.algorithms['cmodel'].initial.usePixelWeights=True

'''Maximum radius used in approximating profile with Gaussians (0=default for this profile) '''
root.measurement.algorithms['cmodel'].initial.maxRadius=0

'''Whether to record the steps the optimizer takes (or just the number, if running as a plugin) '''
root.measurement.algorithms['cmodel'].initial.doRecordHistory=True

'''If the maximum of the gradient falls below this threshold, consider the algorithm converged '''
root.measurement.algorithms['cmodel'].dev.optimizer.gradientThreshold=1e-05

'''steps with reduction radio less than this will decrease the trust radius '''
root.measurement.algorithms['cmodel'].dev.optimizer.trustRegionShrinkReductionRatio=0.25

'''when increase the trust region size, multiply the radius by this factor '''
root.measurement.algorithms['cmodel'].dev.optimizer.trustRegionGrowFactor=2.0

'''steps with length this fraction of the trust radius may increase the trust radius '''
root.measurement.algorithms['cmodel'].dev.optimizer.trustRegionGrowStepFraction=0.8

'''whether to save all iterations for debugging purposes '''
root.measurement.algorithms['cmodel'].dev.optimizer.doSaveIterations=False

'''steps with reduction radio greater than this may increase the trust radius '''
root.measurement.algorithms['cmodel'].dev.optimizer.trustRegionGrowReductionRatio=0.75

'''when reducing the trust region size, multiply the radius by this factor '''
root.measurement.algorithms['cmodel'].dev.optimizer.trustRegionShrinkFactor=0.3333333333333333

'''steps with reduction ratio greater than this are accepted '''
root.measurement.algorithms['cmodel'].dev.optimizer.stepAcceptThreshold=0.0

'''the initial trust region will be set to this value '''
root.measurement.algorithms['cmodel'].dev.optimizer.trustRegionInitialSize=1.0

'''If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method '''
root.measurement.algorithms['cmodel'].dev.optimizer.noSR1Term=False

'''relative step size used for numerical derivatives (added to other steps) '''
root.measurement.algorithms['cmodel'].dev.optimizer.numDiffRelStep=0.0

'''step size (in units of trust radius) used for numerical derivatives (added to relative step) '''
root.measurement.algorithms['cmodel'].dev.optimizer.numDiffTrustRadiusStep=0.1

'''maximum number of steps '''
root.measurement.algorithms['cmodel'].dev.optimizer.maxOuterIterations=500

'''Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold '''
root.measurement.algorithms['cmodel'].dev.optimizer.skipSR1UpdateThreshold=1e-08

'''maximum number of iterations (i.e. function evaluations and trust region subproblems) per step '''
root.measurement.algorithms['cmodel'].dev.optimizer.maxInnerIterations=20

'''If the trust radius falls below this threshold, consider the algorithm converged '''
root.measurement.algorithms['cmodel'].dev.optimizer.minTrustRadiusThreshold=1e-05

'''value passed as the tolerance to solveTrustRegion '''
root.measurement.algorithms['cmodel'].dev.optimizer.trustRegionSolverTolerance=1e-08

'''absolute step size used for numerical derivatives (added to other steps) '''
root.measurement.algorithms['cmodel'].dev.optimizer.numDiffAbsStep=0.0

'''Number of degrees of freedom for the Student's T distribution on ln(radius) '''
root.measurement.algorithms['cmodel'].dev.empiricalPriorConfig.logRadiusNu=50.0

'''Width of the Student's T distribution in ln(radius). '''
root.measurement.algorithms['cmodel'].dev.empiricalPriorConfig.logRadiusSigma=0.45

'''Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T. '''
root.measurement.algorithms['cmodel'].dev.empiricalPriorConfig.logRadiusMu=-1.0

'''Width of exponential ellipticity distribution (conformal shear units) '''
root.measurement.algorithms['cmodel'].dev.empiricalPriorConfig.ellipticitySigma=0.3

'''Softened core width for ellipticity distribution (conformal shear units '''
root.measurement.algorithms['cmodel'].dev.empiricalPriorConfig.ellipticityCore=0.001

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.measurement.algorithms['cmodel'].dev.empiricalPriorConfig.logRadiusMinInner=-6.0

'''Minimum ln(radius) '''
root.measurement.algorithms['cmodel'].dev.empiricalPriorConfig.logRadiusMinOuter=-6.001

'''One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None '''
root.measurement.algorithms['cmodel'].dev.priorSource='EMPIRICAL'

'''Number of Gaussian used to approximate the profile '''
root.measurement.algorithms['cmodel'].dev.nComponents=8

'''ln(radius) at which the softened cutoff begins towards the maximum '''
root.measurement.algorithms['cmodel'].dev.linearPriorConfig.logRadiusMaxInner=3.0

'''The ratio P(logRadiusMinInner)/P(logRadiusMaxInner) '''
root.measurement.algorithms['cmodel'].dev.linearPriorConfig.logRadiusMinMaxRatio=1.0

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.measurement.algorithms['cmodel'].dev.linearPriorConfig.logRadiusMinInner=-6.0

'''Maximum ln(radius) '''
root.measurement.algorithms['cmodel'].dev.linearPriorConfig.logRadiusMaxOuter=3.001

'''Ellipticity magnitude (conformal shear units) at which the softened cutoff begins '''
root.measurement.algorithms['cmodel'].dev.linearPriorConfig.ellipticityMaxInner=2.0

'''Minimum ln(radius) '''
root.measurement.algorithms['cmodel'].dev.linearPriorConfig.logRadiusMinOuter=-6.001

'''Maximum ellipticity magnitude (conformal shear units) '''
root.measurement.algorithms['cmodel'].dev.linearPriorConfig.ellipticityMaxOuter=2.001

'''Name of the Prior that defines the model to fit (a filename in $MEAS_MULTIFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting. '''
root.measurement.algorithms['cmodel'].dev.priorName=''

'''Whether to record the time spent in this stage '''
root.measurement.algorithms['cmodel'].dev.doRecordTime=True

'''Name of the shapelet.RadialProfile that defines the model to fit '''
root.measurement.algorithms['cmodel'].dev.profileName='luv'

'''Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances) '''
root.measurement.algorithms['cmodel'].dev.usePixelWeights=False

'''Maximum radius used in approximating profile with Gaussians (0=default for this profile) '''
root.measurement.algorithms['cmodel'].dev.maxRadius=0

'''Whether to record the steps the optimizer takes (or just the number, if running as a plugin) '''
root.measurement.algorithms['cmodel'].dev.doRecordHistory=True

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['cmodel'].priority=2.5

'''If the maximum of the gradient falls below this threshold, consider the algorithm converged '''
root.measurement.algorithms['cmodel'].exp.optimizer.gradientThreshold=1e-05

'''steps with reduction radio less than this will decrease the trust radius '''
root.measurement.algorithms['cmodel'].exp.optimizer.trustRegionShrinkReductionRatio=0.25

'''when increase the trust region size, multiply the radius by this factor '''
root.measurement.algorithms['cmodel'].exp.optimizer.trustRegionGrowFactor=2.0

'''steps with length this fraction of the trust radius may increase the trust radius '''
root.measurement.algorithms['cmodel'].exp.optimizer.trustRegionGrowStepFraction=0.8

'''whether to save all iterations for debugging purposes '''
root.measurement.algorithms['cmodel'].exp.optimizer.doSaveIterations=False

'''steps with reduction radio greater than this may increase the trust radius '''
root.measurement.algorithms['cmodel'].exp.optimizer.trustRegionGrowReductionRatio=0.75

'''when reducing the trust region size, multiply the radius by this factor '''
root.measurement.algorithms['cmodel'].exp.optimizer.trustRegionShrinkFactor=0.3333333333333333

'''steps with reduction ratio greater than this are accepted '''
root.measurement.algorithms['cmodel'].exp.optimizer.stepAcceptThreshold=0.0

'''the initial trust region will be set to this value '''
root.measurement.algorithms['cmodel'].exp.optimizer.trustRegionInitialSize=1.0

'''If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method '''
root.measurement.algorithms['cmodel'].exp.optimizer.noSR1Term=False

'''relative step size used for numerical derivatives (added to other steps) '''
root.measurement.algorithms['cmodel'].exp.optimizer.numDiffRelStep=0.0

'''step size (in units of trust radius) used for numerical derivatives (added to relative step) '''
root.measurement.algorithms['cmodel'].exp.optimizer.numDiffTrustRadiusStep=0.1

'''maximum number of steps '''
root.measurement.algorithms['cmodel'].exp.optimizer.maxOuterIterations=250

'''Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold '''
root.measurement.algorithms['cmodel'].exp.optimizer.skipSR1UpdateThreshold=1e-08

'''maximum number of iterations (i.e. function evaluations and trust region subproblems) per step '''
root.measurement.algorithms['cmodel'].exp.optimizer.maxInnerIterations=20

'''If the trust radius falls below this threshold, consider the algorithm converged '''
root.measurement.algorithms['cmodel'].exp.optimizer.minTrustRadiusThreshold=1e-05

'''value passed as the tolerance to solveTrustRegion '''
root.measurement.algorithms['cmodel'].exp.optimizer.trustRegionSolverTolerance=1e-08

'''absolute step size used for numerical derivatives (added to other steps) '''
root.measurement.algorithms['cmodel'].exp.optimizer.numDiffAbsStep=0.0

'''Number of degrees of freedom for the Student's T distribution on ln(radius) '''
root.measurement.algorithms['cmodel'].exp.empiricalPriorConfig.logRadiusNu=50.0

'''Width of the Student's T distribution in ln(radius). '''
root.measurement.algorithms['cmodel'].exp.empiricalPriorConfig.logRadiusSigma=0.45

'''Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T. '''
root.measurement.algorithms['cmodel'].exp.empiricalPriorConfig.logRadiusMu=-1.0

'''Width of exponential ellipticity distribution (conformal shear units) '''
root.measurement.algorithms['cmodel'].exp.empiricalPriorConfig.ellipticitySigma=0.3

'''Softened core width for ellipticity distribution (conformal shear units '''
root.measurement.algorithms['cmodel'].exp.empiricalPriorConfig.ellipticityCore=0.001

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.measurement.algorithms['cmodel'].exp.empiricalPriorConfig.logRadiusMinInner=-6.0

'''Minimum ln(radius) '''
root.measurement.algorithms['cmodel'].exp.empiricalPriorConfig.logRadiusMinOuter=-6.001

'''One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None '''
root.measurement.algorithms['cmodel'].exp.priorSource='EMPIRICAL'

'''Number of Gaussian used to approximate the profile '''
root.measurement.algorithms['cmodel'].exp.nComponents=6

'''ln(radius) at which the softened cutoff begins towards the maximum '''
root.measurement.algorithms['cmodel'].exp.linearPriorConfig.logRadiusMaxInner=3.0

'''The ratio P(logRadiusMinInner)/P(logRadiusMaxInner) '''
root.measurement.algorithms['cmodel'].exp.linearPriorConfig.logRadiusMinMaxRatio=1.0

'''ln(radius) at which the softened cutoff begins towards the minimum '''
root.measurement.algorithms['cmodel'].exp.linearPriorConfig.logRadiusMinInner=-6.0

'''Maximum ln(radius) '''
root.measurement.algorithms['cmodel'].exp.linearPriorConfig.logRadiusMaxOuter=3.001

'''Ellipticity magnitude (conformal shear units) at which the softened cutoff begins '''
root.measurement.algorithms['cmodel'].exp.linearPriorConfig.ellipticityMaxInner=2.0

'''Minimum ln(radius) '''
root.measurement.algorithms['cmodel'].exp.linearPriorConfig.logRadiusMinOuter=-6.001

'''Maximum ellipticity magnitude (conformal shear units) '''
root.measurement.algorithms['cmodel'].exp.linearPriorConfig.ellipticityMaxOuter=2.001

'''Name of the Prior that defines the model to fit (a filename in $MEAS_MULTIFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting. '''
root.measurement.algorithms['cmodel'].exp.priorName=''

'''Whether to record the time spent in this stage '''
root.measurement.algorithms['cmodel'].exp.doRecordTime=True

'''Name of the shapelet.RadialProfile that defines the model to fit '''
root.measurement.algorithms['cmodel'].exp.profileName='lux'

'''Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances) '''
root.measurement.algorithms['cmodel'].exp.usePixelWeights=False

'''Maximum radius used in approximating profile with Gaussians (0=default for this profile) '''
root.measurement.algorithms['cmodel'].exp.maxRadius=0

'''Whether to record the steps the optimizer takes (or just the number, if running as a plugin) '''
root.measurement.algorithms['cmodel'].exp.doRecordHistory=True

'''If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this.  If <= 0.0, abort the fit early instead. '''
root.measurement.algorithms['cmodel'].fallbackInitialMomentsPsfFactor=1.5

'''Root name of the FitPsfAlgorithm fields. '''
root.measurement.algorithms['cmodel'].psfName='multishapelet.psf'

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['focalplane'].priority=3.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['centroid.naive'].priority=0.0

'''FIXME! NEVER DOCUMENTED! '''
root.measurement.algorithms['centroid.naive'].background=0.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['shape.hsm.bj'].priority=1.0

'''Mask planes used to reject bad pixels. '''
root.measurement.algorithms['shape.hsm.bj'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

'''Field name for number of deblend children '''
root.measurement.algorithms['shape.hsm.bj'].deblendNChild=''

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['centroid.gaussian'].priority=0.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['flux.aperture.elliptical'].priority=1.899999976158142

'''Maximum number of radial annuli to measure '''
root.measurement.algorithms['flux.aperture.elliptical'].nApertureMax=10

'''vector of radii for apertures (in pixels) '''
root.measurement.algorithms['flux.aperture.elliptical'].radii=[1.0, 1.5625, 2.44140625, 3.814697265625, 5.9604644775390625, 9.313225746154785, 14.551915228366852, 22.737367544323206, 35.52713678800501, 55.51115123125783]

'''Largest aperture for which to use the slow, accurate, sinc aperture code '''
root.measurement.algorithms['flux.aperture.elliptical'].maxSincRadius=10.0

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['variance'].priority=2.0

'''scale factor to apply to shape for aperture '''
root.measurement.algorithms['variance'].scale=5.0

'''mask planes to ignore '''
root.measurement.algorithms['variance'].mask=['DETECTED', 'DETECTED_NEGATIVE']

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['correctfluxes'].priority=3.0

'''List of flux fields that should not be corrected (otherwise all fields in getApCorrRegistry() will be) '''
root.measurement.algorithms['correctfluxes'].ignored=[]

'''Whether to propagate aperture correction uncertainties into flux uncertainties '''
root.measurement.algorithms['correctfluxes'].doPropagateErrors=False

'''Whether to set the general failure flag for a flux when it cannot be aperture-corrected '''
root.measurement.algorithms['correctfluxes'].doFlagApCorrFailures=True

'''Whether to save the per-source per-flux aperture corrections and their errors '''
root.measurement.algorithms['correctfluxes'].doRecordApCorr=True

'''Parameter that sets the sort order for algorithms; lower numbers go first. Typically, priority=0 for centroids, 1 for shapes, and 2 for fluxes. '''
root.measurement.algorithms['shape.hsm.ksb'].priority=1.0

'''Mask planes used to reject bad pixels. '''
root.measurement.algorithms['shape.hsm.ksb'].badMaskPlanes=['BAD', 'SAT', 'INTRP']

'''Field name for number of deblend children '''
root.measurement.algorithms['shape.hsm.ksb'].deblendNChild=''

root.measurement.algorithms.names=['flux.psf', 'flags.pixel', 'shape.hsm.moments', 'flux.aperture', 'flux.naive', 'focalplane', 'flux.gaussian', 'centroid.naive', 'flux.sinc', 'shape.sdss', 'jacobian', 'shape.hsm.regauss', 'flux.kron', 'correctfluxes', 'classification.extendedness', 'skycoord', 'shape.hsm.psfMoments']
'''The seed value to use for random number generation. '''
root.measurement.replaceWithNoise.noiseSeed=0

'''Add ann offset to the generated noise. '''
root.measurement.replaceWithNoise.noiseOffset=0.0

'''How do we choose the mean and variance of the Gaussian noise we generate?
Allowed values:
	variance	Mean = 0, variance = the image's variance
	meta	Mean = 0, variance = the "BGMEAN" metadata entry
	measure	Measure clipped mean and variance from the whole image
 '''
root.measurement.replaceWithNoise.noiseSource='measure'

'''the name of the aperture flux algorithm used for calibration '''
root.measurement.slots.calibFlux='flux.naive'

'''the name of the algorithm used to set the source aperture flux slot '''
root.measurement.slots.apFlux='flux.sinc'

'''the name of the algorithm used to set the source inst flux slot '''
root.measurement.slots.instFlux='flux.gaussian'

'''the name of the algorithm used to set source moments parameters '''
root.measurement.slots.shape='shape.hsm.moments'

'''the name of the centroiding algorithm used to set source x,y '''
root.measurement.slots.centroid='centroid.sdss'

'''the name of the algorithm used to set the source model flux slot '''
root.measurement.slots.modelFlux='flux.gaussian'

'''the name of the algorithm used to set the source psf flux slot '''
root.measurement.slots.psfFlux='flux.psf'

'''When measuring, replace other detected footprints with noise? '''
root.measurement.doReplaceWithNoise=True

'''Measure sources? '''
root.doMeasurement=False

'''What to do when a peak to be deblended is close to the edge of the image
Allowed values:
	None	Field is optional
	ramp	Ramp down flux at the image edge by the PSF
	noclip	Ignore the edge when building the symmetric template.
	clip	Clip the template at the edge AND the mirror of the edge.
 '''
root.deblend.edgeHandling='ramp'

'''Assign stray flux to deblend children.  Implies findStrayFlux. '''
root.deblend.assignStrayFlux=True

'''Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended. '''
root.deblend.maskLimits={'NO_DATA': 0.25}

'''How to split flux among peaks
Allowed values:
	trim	Shrink the parent footprint to pixels that are not assigned to children
	r-to-peak	~ 1/(1+R^2) to the peak
	None	Field is optional
	r-to-footprint	~ 1/(1+R^2) to the closest pixel in the footprint.  CAUTION: this can be computationally expensive on large footprints!
	nearest-footprint	Assign 100% to the nearest footprint (using L-1 norm aka Manhattan distance)
 '''
root.deblend.strayFluxRule='trim'

'''If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up '''
root.deblend.catchFailures=False

'''When the deblender should attribute stray flux to point sources
Allowed values:
	always	Always
	None	Field is optional
	never	Never; stray flux will not be attributed to any deblended child if the deblender thinks all peaks look like point sources
	necessary	When there is not an extended object in the footprint
 '''
root.deblend.strayFluxToPointSources='necessary'

'''Chi-squared per DOF cut for deciding a source is a PSF during deblending (shifted PSF model #2) '''
root.deblend.psfChisq2b=1.5

'''Mask planes to ignore when performing statistics '''
root.deblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

'''Maximum area for footprints before they are ignored as large; non-positive means no threshold applied '''
root.deblend.maxFootprintArea=10000

'''Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied '''
root.deblend.minFootprintAxisRatio=0.0

'''Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied '''
root.deblend.maxFootprintSize=0

'''Chi-squared per DOF cut for deciding a source is a PSF during deblending (un-shifted PSF model) '''
root.deblend.psfChisq1=1.5

'''Find stray flux---flux not claimed by any child in the deblender. '''
root.deblend.findStrayFlux=True

'''Footprints smaller in width or height than this value will be ignored; minimum of 2 due to PSF gradient calculation.
	Valid Range = [2,inf) '''
root.deblend.tinyFootprintSize=2

'''Chi-squared per DOF cut for deciding a source is PSF during deblending (shifted PSF model) '''
root.deblend.psfChisq2=1.5

'''Guarantee that all peaks produce a child source. '''
root.deblend.propagateAllPeaks=False

'''Mask name for footprints not deblended, or None '''
root.deblend.notDeblendedMask='NOT_DEBLENDED'

'''Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited) '''
root.deblend.maxNumberOfPeaks=0

'''When splitting stray flux, clip fractions below this value to zero. '''
root.deblend.clipStrayFluxFraction=0.001

'''Perform calibration? '''
root.doCalibrate=True

'''If True persist background model with background subtracted calexp.          If False persist calexp with the background included. '''
root.persistBackgroundModel=True

'''Write calibration results? '''
root.doWriteCalibrate=False

