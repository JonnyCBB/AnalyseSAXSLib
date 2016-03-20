import numpy as np
#import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import math
from scipy import ndimage
from scipy import signal
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy import exp
from skimage import restoration


class Beam():

    # #Class Attributes
    # beamFlux = 1e12 #photons per second. ###NEED TO CONFIRM THE REAL VALUE
    # beamEnergy = 12.66 #keV
    # doseType = "DWD" #the dose type parsed from the RADDOSE-3D log file.
    # beamPixelSize = [0.3027, 0.2995] #THIS SHOULDN'T BE A CLASS ATTRIBUTE [0.3027, 0.2995]
    # extraX = 0
    # extraY = 0

    #Class constructor
    def __init__(self, beamArray, pgmFileName):
        """Beam object constructor

        INPUTS:
            beamArray           -A 2D numpy array of integer values as a spatially resolved representation
                                    of relative intensities. The values should be between 0 and 255 to be
                                    compatible with the .pgm file format.
            pgmFileName         -A string containing the name of the pgm file to point to for this particular beam instance.
        """
        self.beamArray = beamArray
        self.pgmFileName = writePGMFile(self.beamArray,pgmFileName)
        print "SUCCESS :-) Beam has been successfully created"
        print "****************************************************"
        print

    #Note that class methods can access class attributes but not instance
    #attributes.
    @classmethod
    def initialiseBeamFromApMeas(cls, beamApMeasXFilename, beamApMeasYFilename, beamPostProcessingType, apDiameter, apStep, outputPGMFileName,
                                fitGauss, deconvolve, deblur):
        """Create a beam object from aperture scan measurements

        INPUTS:
            beamApMeasXFilename     -string with the location of the aperture scan measurements in the horizontal (x)
                                        direction.
            beamApMeasYFilename     -string with the location of the aperture scan measurements in the vertical (y)
                                        direction.
            beamPostProcessingType  -string giving the type of processing that should be carried out on the beam
                                        array once the deconvolution has taken place.
                                        Type "?BeamModule.beamPostProcessManip" for more details about how this parameter
                                        is used
            apDiameter              -The diameter of the aperture in microns
            apStep                  -The incremental position at which consecutive i_pin readings
                                        are taken (in microns).
            outputPGMFileName       -A string containing the name of the pgm file that will be written during the actual
                                        object initialisation.

        OUTPUTS:
            beamFromApMeas          -A beam object.
        """

        print "****************************************************"
        print "Creating a beam object from aperture measurements..."
        print
        beamArray = generateBeamFromApMeas(beamApMeasXFilename, beamApMeasYFilename, beamPostProcessingType, apDiameter, apStep,
                                          fitGauss, deconvolve, deblur)
        print "Beam array has been generated."
        print "Now centering the beam array..."
        beamCentroid = findCentroid(beamArray)
        beamArray = cropBeamArray(beamCentroid,beamArray)
        print "Beam centering successful."
        beamFromApMeas = cls(beamArray, outputPGMFileName)
        return beamFromApMeas

    @classmethod
    def initialiseBeamFromPNG(cls, pngImage, redWeightValue, greenWeightValue, blueWeightValue, outputPGMFileName, collimation, beamPixelSize):
        """Create a beam object from a PNG file.

        INPUTS:
            pngImage           -The path to a png file that contains an image of an X-ray beam as a string.
            redWeightValue     -A scalar (float) value giving the weight of the red pixels in the png image
                                for the conversion to grayscale.
            greenWeightValue   -A scalar (float) value giving the weight of the green pixels in the png image
                                for the conversion to grayscale.
            blueWeightValue    -A scalar (float) value giving the weight of the blue pixels in the png image
                                for the conversion to grayscale.
            outputPGMFileName  -A string containing the name of the pgm file that will be written during the actual
                                        object initialisation.

        OUTPUTS:
            beamFromPNG     -A beam object.
        """

        print "****************************************************"
        print "Creating a beam object from PNG image..."
        print
        beamArray = generateBeamFromPNG(pngImage, redWeightValue, greenWeightValue, blueWeightValue)
        print "Beam array has been generated."
        print "Now centering the beam array..."
        beamCentroid = findCentroid(beamArray)
        print "Beam centering successful."
        print "Cropping Beam Array so the beam is central..."
        beamArray = cropBeamArray(beamCentroid,beamArray)
        print "Finished cropping"
        print "Subtracting background from the beam image..."
        beamArray = subtractBackground(beamArray, collimation[0], collimation[1], beamPixelSize[0], beamPixelSize[1])
        print "Background substraction successful"
        beamFromPNG = cls(beamArray, outputPGMFileName)
        return beamFromPNG

    @classmethod
    def initialiseBeamFromPGM(cls, pgmImageFile):
        """Create a beam object from a PGM file

        INPUTS:
            pgmImageFile       -A string giving the location of the pgm file containing the beam image
            outputPGMFileName  -A string containing the name of the pgm file that will be written during the actual
                                        object initialisation.

        OUTPUTS:
            beamFromPGM     -A beam object.
        """

        # print "****************************************************"
        # print "Creating a beam object from PGM file..."
        # print
        # beamArray = generateBeamFromPGM(pgmImageFile)
        # beamCentroid = findCentroid(beamArray)
        # beamArray = cropBeamArray(beamCentroid,beamArray)
        # slits = np.zeros(2, dtype=np.float)
        # cameraSettings = (0,0,0,0,slits)
        # beamFromPGM = cls(beamArray, cls.beamFlux, cls.beamEnergy, cls.beamPixelSize, pgmImageFile, cameraSettings, cls.doseType)
        # return beamFromPGM

    # def formRADDOSE3DBeamInputString(self, pixelSize):
    #     """Method that forms a string containing the Beam parameters that is
    #     suitable to be transferred directly into a RADDOSE-3D input file.
    #
    #     INPUTS:
    #         pixelSize       -A 2-element list, array or tuple, containing the horizontal and vertical values for the
    #                             collimation.
    #
    #     OUTPUTS:
    #         beamParameters  -A string containing the Beam block inforamtion formatted to be written directly into a
    #                             RADDOSE-3D input file.
    #     """
    #     #Create strings for each of the lines required in the beam block of the RADDOSE-3D input file.
    #     beamLine        = "Beam"
    #     typeLine        = "Type ExperimentalPGM"
    #     fileLine        = "File {pgmFile}".format(pgmFile=self.pgmFileName)
    #     pixelLine       = "PixelSize {horz} {vert}".format(horz=str(pixelSize[0]), vert=str(pixelSize[1]))
    #     fluxLine        = "Flux {flux}".format(flux=self.beamFlux)
    #     energyLine      = "Energy {energy}".format(energy=self.beamEnergy)
    #     collimationLine = "Collimation Rectangular {vert} {horz}".format(horz=self.collimation[0], vert=self.collimation[1])
    #
    #     #Put all of the lines into one formatted string.
    #     beamParameters  = "{beam}\n{type}\n{file}\n{pixel}\n{flux}\n{energy}\n{collimation}\n".format(
    #     beam=beamLine, type=typeLine, file=fileLine, pixel=pixelLine, flux=fluxLine,
    #     energy=energyLine, collimation=collimationLine)
    #
    #     return beamParameters

def generateBeamFromApMeas(beamApMeasXFilename, beamApMeasYFilename, beamPostProcessingType, apDiameter, apStep,
                          fitConvolvedGaussian, deconvole, deblurBeam):
    """Create a beam array from aperture scan measurements

    INPUTS:
        beamApMeasXFilename     -string with the location of the aperture scan measurements in the horizontal (x)
                                    direction.
        beamApMeasYFilename     -string with the location of the aperture scan measurements in the vertical (y)
                                    direction.
        beamPostProcessingType  -string giving the type of processing that should be carried out on the beam
                                    array once the deconvolution has taken place.
        apertureDiameter        -The diameter of the aperture in microns
        apertureStep            -The incremental position at which consecutive i_pin readings
                                    are taken (in microns).

    OUTPUTS:
        beamArray               -A 2D numpy array of integer values as a spatially resolved representation
                                    of relative intensities. The values should be between 0 and 255 to be
                                    compatible with the .pgm file format.
    """
    #load the .dat file and store aperture data in a 2D numpy array.
    #Then split the 2D array into two 1D arrays, one to store the i_pin readings.
    #and one to store the aperture positions.
    #Do this for both the vertical and horizontal directions.
    apertureX = np.loadtxt(beamApMeasXFilename,skiprows=7)
    apertureXPosition = apertureX[:,0]
    apertureXMeasurement = apertureX[:,1]

    apertureY = np.loadtxt(beamApMeasYFilename,skiprows=7)
    apertureYPosition = apertureY[:,0]
    apertureYMeasurement = apertureY[:,1]

    #Get maximum i_pin reading
    maxIpin = max(apertureXMeasurement.max(), apertureYMeasurement.max())

    #If Gaussian fit is selected then fit a Gaussian profile to the beam.
    if "gaussfitconv" in fitConvolvedGaussian:
        #Fit Gaussian to the x direction data
        meanGuessX = (apertureXPosition[-1] + apertureXPosition[0])/2.0
        sigmaGuessX = 1
        gaussFitXparams, pcovariance = curve_fit(gauss,apertureXPosition,apertureXMeasurement,
                                           p0=[apertureXMeasurement.max(), meanGuessX, sigmaGuessX])

        #Fit Gaussian to the y direction data
        meanGuessY = (apertureYPosition[-1] + apertureYPosition[0])/2.0
        sigmaGuessY = 1
        gaussFitYparams, pcovariance = curve_fit(gauss,apertureYPosition,apertureYMeasurement,
                                           p0=[apertureYMeasurement.max(), meanGuessY, sigmaGuessY])

        #Create a 2d Gaussian
        X, Y = np.meshgrid(apertureXPosition, apertureYPosition)
        muX = gaussFitXparams[1]
        muY = gaussFitYparams[1]
        sigmaX = gaussFitXparams[2]
        sigmaY = gaussFitYparams[2]
        convolvedBeamArray = maxIpin * exp(-((X-muX)**2 / (2 * sigmaX**2) + (Y-muY)**2 / (2 * sigmaY**2)))

    #If the Gaussian fit wasn't selected then generate the beam by averaging the data.
    else:
        #Create temporary 2D beam arrays using both row and column wise scaling
        tempBeamArrayX = beamScalingRowWise(apertureYMeasurement,apertureXMeasurement)
        tempBeamArrayY = beamScalingColWise(apertureYMeasurement,apertureXMeasurement)

        #Take an average of the two temporary beam arrays to get a single convolved beam array
        convolvedBeamArray = (tempBeamArrayX + tempBeamArrayY) / 2

    ##############################################################
    # Create plot
    ##############################################################
#     rowNum = math.floor(convolvedBeamArray.shape[0]/2.0)
#     colNum = round(convolvedBeamArray.shape[1]/2.0)
#     fig = plt.figure()
#     plt.subplot(211)
#     plt.plot(apertureXPosition, convolvedBeamArray[rowNum,0:], 'b-', apertureXPosition, apertureXMeasurement, 'ro')

#     plt.subplot(212)
#     plt.plot(apertureYPosition, convolvedBeamArray[0:,colNum], 'b-', apertureYPosition, apertureYMeasurement, 'ro')
#     plt.show()

    if deconvole:
        #Create the aperture Point Spread Function.
        aperturePSF = createAperturePSF(apDiameter,apStep)

        #Deconvolve the beam image
        blurredBeamTuple = restoration.unsupervised_wiener(convolvedBeamArray, aperturePSF)
        blurredBeamArray = blurredBeamTuple[0]     #Get beam array

        if "smoothweiner" in deblurBeam:
            #Deblur beam with a Weiner filter
            initialBeamNoiseGuess = 1
            res = minimize(lambda beamNoise: deblurBeamObjectiveFunction(beamNoise,blurredBeamArray,aperturePSF,
                                                                        apertureXMeasurement,apertureYMeasurement),
                           initialBeamNoiseGuess, method='nelder-mead', options={'xtol': 1e-8, 'disp': True})

            #Deblur the image
            deconvolvedBeamArray = signal.wiener(blurredBeamArray,aperturePSF.shape,res.x[0])
            deblurredBeamArray = deconvolvedBeamArray
        elif "smoothgauss" in deblurBeam:
            #Fit gaussian to the blurred beam array
            params = fitgaussian(blurredBeamArray)
            fit = gauss2d(*params)
            deblurredBeamArray = fit(*np.indices(blurredBeamArray.shape))
        else:
            deblurredBeamArray = blurredBeamArray

        resultingBeamArray = deblurredBeamArray #Assign deblurred beam to another variable.
    else:
        resultingBeamArray = convolvedBeamArray

    #Apply post processing on the beam array
    processedBeamArray = beamPostProcessManip(resultingBeamArray,beamPostProcessingType)

    #Scale the beam array
    beamArray = scaleArray(processedBeamArray)

    ##############################################################
    # Create plot
    ##############################################################
#     rowNum = math.floor(beamArray.shape[0]/2.0)
#     colNum = round(beamArray.shape[1]/2.0)
#     fig = plt.figure()
#     plt.subplot(131)
#     plt.imshow(beamArray)

#     plt.subplot(132)
#     plt.plot(apertureXPosition[0:len(beamArray[rowNum,0:])], beamArray[rowNum,0:], 'b-')

#     plt.subplot(133)
#     plt.plot(apertureYPosition, beamArray[0:,colNum], 'b-')
#     plt.show()

    return beamArray

def generateBeamFromPNG(pngImage, redWeightValue, greenWeightValue, blueWeightValue):
    """Function that generates a beam array from a PNG file.

    INPUTS:
        pngImage           -The path to a png file that contains an image of an X-ray beam as a string.
        redWeightValue     -A scalar (float) value giving the weight of the red pixels in the png image
                            for the conversion to grayscale.
        greenWeightValue   -A scalar (float) value giving the weight of the green pixels in the png image
                            for the conversion to grayscale.
        blueWeightValue    -A scalar (float) value giving the weight of the blue pixels in the png image
                            for the conversion to grayscale.

    OUTPUTS:
        beamArray          -A 2D numpy array of integer values as a spatially resolved representation
                            of relative intensities. The values should be between 0 and 255 to be
                            compatible with the .pgm file format.
    """

    rgbBeamImage = mpimg.imread(pngImage)     #read the png file
    grayscaleBeamImage = rgb2Grayscale(rgbBeamImage, redWeightValue, greenWeightValue, blueWeightValue)     #convert from rgb values to grayscale
    beamArray = scaleArray(grayscaleBeamImage)     #scale the array whilst converting floats to integer values
    return beamArray     #return beam array

def generateBeamFromPGM(pgmImageFile):
    """Function that generates a beam array from a PGM file

    INPUTS:
        pgmImageFile     -A string giving the location of the pgm file containing the beam image

    OUTPUTS:
        beamArray        -2D numpy array of ints representing a spatial relative intensity distribution.
    """

    #Print an informative string to the console.
    print "Generating a 2D array from file: \"" + pgmImageFile + "\""

    #Local function variables
    maxValueLine = False     #boolean variable to determine if the line corresponding to the max pixel value has been reached
    pixelCounter = 0     #variable to count the number of pixel values
    beamListCounter = 0     #variable to count the number pixel elements in the pixel value list

    pgmfile = open(pgmImageFile,'r') #Open pgm file for reading

    #Loop through each line in the pgm file
    for line in pgmfile:
        #If the line begins with a 'P' (magic number line) or a '#' (comment line) then do nothing. Else...
        if (line[0] != 'P' and line[0] != '#'):
            #If the line has two values then these represent the pixel width and height
            if len(line.split()) == 2:
                beamArrayWidth = int(line.split()[0])     #extract the width
                beamArrayHeight = int(line.split()[1])     #extract the height
                beamList = np.zeros(beamArrayWidth * beamArrayHeight)     #calculate the total number of pixel values and preallocate an array to store them
                maxValueLine = True     #The next line is the max value line which we can ignore
            elif maxValueLine == True:     #If the current line is the max value line...
                maxValueLine = False    #...then ignore it but set the max value line checker to false for the subsequent lines
            else:
                #The rest of the lines represent pixel values so we want to store these
                beamList[pixelCounter] = int(line)
                pixelCounter += 1     #increment pixel counter

    pgmfile.close()     #close the file

    #preallocate the beam array
    beamArray = np.zeros((beamArrayHeight,beamArrayWidth), dtype=np.int)

    #Loop through each value in the beam array and insert the corresponding pixel value.
    for row in xrange(0,beamArrayHeight):
        for col in xrange(0,beamArrayWidth):
            beamArray[row,col] = beamList[beamListCounter]
            beamListCounter += 1

    print "Successful: Beam has been generated."

    return beamArray     #Return the beam array.

def beamPostProcessManip(beamArray,processType):
    """Function that applies some processing to the aperture measurements
    prior to using them to create a beam.

    INPUTS:
        beamArray      - The deconvolved beam array as a 2D numpy array of floats.
        ProcessType    - The type of processing performed input as a string. There are 3
                        types of processing options available:
                        1) "shift" - If the minimum beam array reading is negative than it will
                        add the absolute value of that reading to all values in the beamArray
                        array. This has the effect of shifting all the values up so every
                        value is above zero.
                        2)"threshold" - Any negative beam array reading is set to zero.
                        3) "" - Any other string input will not do any preprocessing.

    OUTPUTS:
        beamArray         - A 2D numpy array of floats containing the processed beam array values.
    """
    if processType == "positive":
        print 'Post processing type used on the beam array: "{}" '.format(processType)
        beamArray = beamArray[beamArray >= 0]     #Only use non-negative values (NOT RECOMMENDED)
    elif processType == "shift":
        print 'Post processing type used on the beam array: "{}" '.format(processType)
        minValue = beamArray.min()         #Find minimum value
        if minValue < 0:
            beamArray = math.fabs(minValue) + beamArray     #If minimum value is negative then shift all the values up so they are positive
    elif processType == "threshold":
        print 'Post processing type used on the beam array: "{}" '.format(processType)
        beamArray[beamArray < 0] = 0         #Set negative values to zero
    else:
        print 'No post processing has been performed on the beam array measurements.'
    return beamArray

def beamScalingRowWise(rows,cols):
    """Scales the i_pin measurements in the horizontal direction according to the normalised
    values in the vertical directon to generate a 2D beam array.

    INPUTS:
        rows         -i_pin readings in the vertical direction as a numpy array of floats
        cols         -i_pin readings in the horizontal direction as a numpy array of floats

    OUTPUTS:
        beamArray    -2D numpy array of floats representing an expected set of i_pins values
                        if the readings were taken across the 2D area
    """
    print 'Performing row-wise scaling of i_pin measurements to generate temporary beam array'
    beamArray = np.zeros((rows.size,cols.size), dtype=np.float64)     #Preallocate beam array
    maxValue = rows.max()             #Get the maximum i_pin value in the vertical direction

    #For each row
    for i in xrange(0,rows.size):
        beamArray[i,:] = (rows[i]/maxValue) * cols  #Scale the horizontal i_pin measurements by the corresponding normalised vertical i_pin reading
    return beamArray

def beamScalingColWise(rows,cols):
    """Scales the i_pin measurements in the vertical direction according to the relative
    values in the horizontal directon to generate a 2D beam array.

    INPUTS:
        rows         -i_pin readings in the vertical direction as a numpy array of floats
        cols         -i_pin readings in the horizontal direction as a numpy array of floats

    OUTPUTS:
        beamArray    -2D numpy array of floats representing an expected set of i_pins values
                        if the readings were taken across the 2D area
    """
    print 'Performing column-wise scaling of i_pin measurements to generate temporary beam array'
    beamArray = np.zeros((rows.size,cols.size), dtype=np.float64)     #Preallocate beam array
    maxValue = cols.max()          #Get the maximum i_pin value in the horizontal direction

    #for each column
    for i in xrange(0,cols.size):
        beamArray[:,i] = (cols[i]/maxValue) * rows  #Scale the vertical i_pin measurements by the corresponding normalised horizontal i_pin reading
    return beamArray

def createAperturePSF(apertureDiameter,apertureStep):
    """Creates a matrix that acts as a Point Spread Function (PSF) of the aperture used to
    generate the i_pin readings.

    INPUTS:
        apertureDiameter        -The diameter of the aperture in microns
        apertureStep            -The incremental position at which consecutive i_pin readings
                                    are taken (in microns).

    OUTPUTS:
        psf                     - A 2D numpy array of integers, either 1's or 0's. 1's represent
                                    the area covered by the aperture.
    """
    print 'Generating the Point Spread Function of the {} micron diameter aperture with a measurement step of {} microns'.format(apertureDiameter,apertureStep)
    #Calculate the aperture radius from the diameter
    apertureRadius = apertureDiameter / 2.0

    ######################################################################
    ### Calculate from the radius how big the psf matrix should be and set all
    ### Values to zeros.
    ###
    #The matrix consists of all surrounding points from the centre of the
    #aperture that are within the area of the diameter - i.e. these matrix
    #elements could potentially contribute to the signal reading
    if apertureRadius%2 == 0:
        psf = np.zeros((apertureRadius + 1, apertureRadius + 1), dtype=np.int)
    else:
        psf = np.zeros((math.floor(apertureRadius), math.floor(apertureRadius)), dtype=np.int)

    ######################################################################
    #Get dimensions of the matrix (it should be a square matrix but we'll
    #take both dimensions anyway)
    dimensionsOfKernel = psf.shape
    #Get the index of the central point of the aperture psf
    indexOfPSFCentre = [dimensionsOfKernel[0] / 2, dimensionsOfKernel[1] / 2]
    #Set the central position in the psf to 1
    psf[indexOfPSFCentre[0], indexOfPSFCentre[1]] = 1

    ######################################################################
    ### Find the indices of the points that lie within the area of the aperture

    #Calculate the Euclidean distance of each point from the centre of the aperture
    distanceMatrix = ndimage.distance_transform_edt(psf==0,sampling=[apertureStep,apertureStep])
    #Find the points that lie within the area of the aperture
    booleanMatrix = distanceMatrix <= apertureRadius

    #Loop through matrix and set points that lie within the aperture to equal 1
    #and points that lie outside the aperture to equal 0
    for i in xrange(0,dimensionsOfKernel[0]):
        for j in xrange(0,dimensionsOfKernel[1]):
            if booleanMatrix[i,j] == True:
                psf[i,j] = 1

    return psf     #Return the point spread function

##########################################################
def simulateApertureScans(beamArray,psf):
    """Function that uses a beam array and a point spread function representing the aperture
    to simulate the aperture scans carried out on the I02 beamline at Diamond Light Source.

    INPUTS:
        beamArray         -A 2D numpy array of floats that represents the theoretical deconvolved i_pin
                            readings at each point in the space.
        psf               -A 2D numpy array of floats that represents the aperture used to measure the
                            beam intensity.

    OUTPUTS:
        simApScanX        -1D numpy array of floats containing i_pin readings from the simulated
                            aperture scan in the horizontal direction.
        simApScanY        -1D numpy array of floats containing i_pin readings from the simulated
                            aperture scan in the vertical direction.
    """

    #Calculate total number of elements to be added in each dimension of beam matrix
    matrixBuffer = int(2*math.floor(psf.shape[0]/2))

    #Preallocate the buffered matrix. It's buffered with zeros around the outside
    bufferedBeamMatrix = np.zeros((beamArray.shape[0] + matrixBuffer, beamArray.shape[1] + matrixBuffer), dtype=np.float)

    #Check if the beam matrix has been cropped in the second dimension by the deconvolution process
    if bufferedBeamMatrix.shape[0] == bufferedBeamMatrix.shape[1]:
        beamArrayCrop = 0
    else:
        beamArrayCrop = bufferedBeamMatrix.shape[0] - bufferedBeamMatrix.shape[1]

    #Fill bufferedBeamMatrix with elements from the actual beam array.
    for i in xrange(matrixBuffer / 2, bufferedBeamMatrix.shape[0] - (matrixBuffer / 2)):
        for j in xrange(matrixBuffer / 2, bufferedBeamMatrix.shape[1] - (matrixBuffer / 2)):
            bufferedBeamMatrix[i,j] = beamArray[i - (matrixBuffer / 2), j - (matrixBuffer / 2)]

    #Preallocate arrays to contain simulated aperture scan measurements
    simApScanX = np.zeros(beamArray.shape[1], dtype=np.float)
    simApScanY = np.zeros(beamArray.shape[0], dtype=np.float)

    #Get the central element for horizontal (X) and the vertical (Y) scans
    apScanRow = int(math.floor(beamArray.shape[0] / 2) + matrixBuffer / 2)
    apScanCol = int(math.floor((beamArray.shape[1] + beamArrayCrop) / 2) + matrixBuffer / 2)

    #Simulate horizontal aperture scan
    for col in xrange(matrixBuffer / 2, bufferedBeamMatrix.shape[1] - (matrixBuffer / 2)):
        for i in xrange(0,psf.shape[0]):
            for j in xrange(0,psf.shape[1]):
                a = i - matrixBuffer / 2
                b = j - matrixBuffer / 2
                simApScanX[col - matrixBuffer] += bufferedBeamMatrix[apScanRow + a, col + b] * psf[i,j]

    #Simulate vertical aperture scan
    for row in xrange(matrixBuffer / 2, bufferedBeamMatrix.shape[0] - (matrixBuffer / 2)):
        for i in xrange(0,psf.shape[0]):
            for j in xrange(0,psf.shape[1]):
                a = i - matrixBuffer / 2
                b = j - matrixBuffer / 2
                simApScanY[row - matrixBuffer] += bufferedBeamMatrix[row + a, apScanCol + b] * psf[i,j]

    return simApScanX, simApScanY

def rootMeanSquaredDeviation(xPredicted,xMeasured):
    """Find the root mean squared deviation between two 1D numpy arrays.
        Note: both arrays must be the same size.

    INPUTS:
        xPredicted        -1D numpy array of floats. These are the predicted values
        xMeasured         -1D numpy array of floats. These are the measured values

    OUTPUTS:
        rmsd              -A scalar float value representing the root mean squared deviation
    """
    rmsd = np.sqrt(np.sum(np.square(xPredicted - xMeasured)))

    return rmsd


def deblurBeamObjectiveFunction(noiseRatio,beamArray,psf,actualApMeasurementX,actualApMeasurementY):
    """Objective function used to deblur the deconvolved beam image.
    An objective funtion is a function whos output is required to be optimal (in this case our optimal value
    is the minimal one) by some optimisation routine. One, or many arguments can be altered by the
    optimisation routine. In this case that parameter is the noiseRatio.

    INPUTS:
        noiseRatio            -A scalar value that represents the noise-power term in the wiener deconvolution
                                function
        beamArray             -A 2D numpy array of floats that represents the theoretical deconvolved i_pin
                                readings at each point in the space.
        psf                   -A 2D numpy array of floats that represents the aperture used to measure the
                                beam intensity.
        actualApMeasurementX  -Measured i_pin readings in the horizontal direction as a 1D numpy array of floats
        actualApMeasurementY  -Measured i_pin readings in the vertical direction as a 1D numpy array of floats

    OUTPUTS:
        totalRMSD             -The sum of the root mean squared deviations of the theoretical and measured
                                i_pin readings.
    """
    #Deblur the image
    deconvolvedBeamArray = signal.wiener(beamArray,psf.shape,noiseRatio)

    #simulate the aperture scans
    simulatedApMeasurementsX,simulatedApMeasurementsY = simulateApertureScans(deconvolvedBeamArray,psf)

    #Calculate the root mean squared deviations (rmsd)
    rmsdX = rootMeanSquaredDeviation(simulatedApMeasurementsX, actualApMeasurementX)
    rmsdY = rootMeanSquaredDeviation(simulatedApMeasurementsY, actualApMeasurementY)
    ######### I used to need the line below (as opposed to the rmsdX above),
    ######### however since changing distribution to Anacdonda the code below
    ######### throws errors because there is a mismatch in the dimensions.
    #rmsdX = rootMeanSquaredDeviation(simulatedApMeasurementsX, actualApMeasurementX[0:-1])

    #Add the rmsds
    totalRMSD = rmsdX + rmsdY

    return totalRMSD

def writePGMFile(beamArray,fileName):
    """Function that writes a plain text PGM file to the specified location.
    Note: The magic number at the header of the file gives the type of PGM.
    In this case the magic number is P2.

    INPUTS:
        beamArray         -A 2D numpy array of integer values as a spatially resolved representation
                            of relative intensities. The values should be between 0 and 255 to be
                            compatible with the .pgm file format.
        fileName          -Name (and location if the full file path is given) of the PGM to be created
                            from the beam array.

    OUTPUTS:
        fileName          -The location of the pgm file specified by the user.

        Note: the other output is the actual pgm file.
    """
    #Create header lines for the pgm file
    magicNumberLine = "P2\n"
    commentLine = "# CREATOR: Diamond Light Source Beamline I02\n"
    arrayHeight, arrayWidth = beamArray.shape     #Get array dimensions
    maxPixelValue = 255     #Get max value of array

    #Print informative string
    print "Writing a PGM file: \"" + fileName + "\" to current directory..."

    #Open a pgm file for writing
    with open(fileName,"w") as pgmFile:
        #Write header lines to pgm file
        pgmFile.write(magicNumberLine)
        pgmFile.write(commentLine)
        pgmFile.write(str(arrayWidth) + " " + str(arrayHeight) + "\n")
        pgmFile.write(str(maxPixelValue) + "\n")
        #Loop through the beam array and insert each element value into the pgm file.
        for row in xrange(0,arrayHeight):
            for col in xrange(0,arrayWidth):
                pgmFile.write(str(beamArray[row,col]) + "\n")

    print "Finished writing PGM file."

    return fileName

def rgb2Grayscale(imageFile,redWeight,greenWeight,blueWeight):
    """Function to convert rgb images to grayscale

    INPUTS:
        imageFile              -An N x M x 3 (3D) array of floats representing the RGB values of an image.
        redWeight              -A scalar (float) value giving the weight of the red pixels in the png image
                                for the conversion to grayscale.
        greenWeight            -A scalar (float) value giving the weight of the green pixels in the png image
                                for the conversion to grayscale.
        blueWeight             -A scalar (float) value giving the weight of the blue pixels in the png image
                                for the conversion to grayscale.

    OUTPUTS:
        grayscaleimage         -An N x M (2D) array of floats representing the weighted average grayscale values
                                of the image.
    """

    #Extract the red, green and blues pixels from the image
    redPixels, greenPixels, bluePixels = imageFile[:,:,0], imageFile[:,:,1], imageFile[:,:,2]
    #Calculate the weighted average
    grayscaleImage = (redWeight * redPixels) + (greenWeight * greenPixels) + (blueWeight * bluePixels)

    return grayscaleImage     #return the grayscale image

def scaleArray(array):
    """ Function to scale the values of the array to range between 0 and 255.
    This function also converts all values to integer values so it's compatible
    with the pgm image file format.

    INPUTS:
        array           -A numpy array of floats/ints

    OUTPUTS:
        scaledArray     -A scaled version of the input array. Values in the array are
                            converted to integers.
    """
    scalingValue = 255/array.max()
    scaledArrayAsFloats = np.around(array * scalingValue)
    scaledArray = scaledArrayAsFloats.astype(int)
    return scaledArray

def parseBeamImageFilename(imageFilename):
    """Parse the image file name to get the information about the camera settings used to capture the image

    INPUTS:
        imageFilename         -String pointing to the png image file.

    OUTPUTS:
        (Tuple of values)     -The output is a tuple of values giving various camera settings set during
                                image capture. These are:
                                zoom
                                gain
                                exposureTime
                                transmission
                                slits
    """

    splitFilename = imageFilename.split("_")     #Split the filename by underscores
    slits = np.zeros(2, dtype=np.float)          #Preallocate array to store the slit dimensions
    #Loop through each split entry of the filename
    for entry in splitFilename:
        if entry != "beam":     #ignore the "beam" entry
            if ".png" in entry:             #check if the entry contains the file extension. If so then remove it.
                splitEntry = entry.split(".")
                entry = splitEntry[0]
            if entry.isdigit():             #check if the entry only contains numbers. If so then it's the vertical slit distance
                slits[1] = float(entry)
            else:                           #else check the different camera settings and assign them to the correct variable.
                firstLetter = entry[0]
                if firstLetter == 'z':
                    zoom = float(entry[1:])
                elif firstLetter == 'g':
                    gain = float(entry[1:])
                elif firstLetter == 'e':
                    exposureTime = float(entry[1:])
                elif firstLetter == 't':
                    transmission = float(entry[1:])
                elif firstLetter == 's':
                    slits[0] = float(entry[1:])

    return(zoom,gain,exposureTime,transmission,slits)     #Return tuple

def findCentroid(array):
    """Find centroid of a numpy array.

    INPUTS:
        array         -a 2D numpy array of floats. In this module this array represents the beam.

    OUTPUTS:
        xcen, ycen    -a 2 element tuple containing the x and y coordinates of the centroid of the array.
    """
    h, w = array.shape
    ygrid, xgrid  = np.mgrid[0:h:1, 0:w:1]
    xcen, ycen = xgrid[array == 255].mean(), ygrid[array == 255].mean()
    return xcen, ycen

def cropBeamArray(arrayElementCoordinates, array):
    """Crop the beam array so that the given array element coordinates are central in the resulting beam
    array.
    Note: the leftIndex variable should have it's named swapped with topIndex and the same
    goes for rightIncdex and bottomIndex but I just can't be bothered to be honest.

    INPUTS:
        arrayElementCoordinates     -A tuple, list or array of floats corresponding to the pixel coordinates
                                        which will be central in the resulting cropped image.
        array                       -A 2D numpy array of floats corresponding to an array of the original beam.

    OUTPUTS:
        croppedArray                -A 2D numpy array of floats corresponding to a cropped array of the original
                                        beam.

    """
    xCoord, yCoord = arrayElementCoordinates[0], arrayElementCoordinates[1]
    arrayHeight, arrayWidth = array.shape

    if math.fabs(yCoord - arrayHeight/2.0) < 1:
        leftIndex = 0
        rightIndex = arrayHeight
    else:
        if arrayHeight/2.0 < yCoord:
            leftIndex = math.floor(yCoord - arrayHeight/2.0)
        else:
            leftIndex = 0
            rightIndex = math.ceil(2 * yCoord)

    if math.fabs(xCoord - arrayWidth/2.0) < 1:
        topIndex = 0
        bottomIndex = arrayWidth
    else:
        if arrayWidth/2.0 < xCoord:
            topIndex = math.floor(xCoord - arrayWidth/2.0)
        else:
            topIndex = 0
            bottomIndex = math.ceil(2 * xCoord)

    if leftIndex != 0 and topIndex != 0:
        croppedArray = array[leftIndex:,topIndex:]
    elif leftIndex != 0:
        croppedArray = array[leftIndex:,topIndex:bottomIndex]
    elif topIndex != 0:
        croppedArray = array[leftIndex:rightIndex,topIndex:]
    else:
        croppedArray = array[leftIndex:rightIndex,topIndex:bottomIndex]

    return croppedArray

def gauss(x, a, mu, sigma):
    return a * exp(-(x-mu)**2 / (2 * sigma**2))

def gauss2d(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gauss2d(*p)(*np.indices(data.shape)) - data)
    p, success = leastsq(errorfunction, params)
    return p

def createWindowAroundCentroid(array, xDistance, yDistance, pixelSizeX, pixelSizeY):
    """Create a window that surrounds the centroid of the array
    """
    xCen, yCen = findCentroid(array)
    pixelsInXDir = xDistance/(2.0 * pixelSizeX)
    pixelsInYDir = yDistance/(2.0 * pixelSizeY)
    xMin = xCen - pixelsInXDir
    xMax = xCen + pixelsInXDir
    yMin = yCen - pixelsInYDir
    yMax = yCen + pixelsInYDir
    return xMin, xMax, yMin, yMax

def subtractBackground(array, xDistance, yDistance, pixelSizeX, pixelSizeY):
    """Function that subtracts an average background from the beam array
    based on the size of the slits.
    """
    window = createWindowAroundCentroid(array, xDistance, yDistance, pixelSizeX, pixelSizeY)
    xMin = np.floor(window[0])
    xMax = np.ceil(window[1])
    yMin = np.floor(window[2])
    yMax = np.ceil(window[3])

    beamArray = array
    backgroundList = []
    #Loop through each element in the beam array and append background elements
    #in a list
    for i in xrange(0,beamArray.shape[0]):
        for j in xrange(0,beamArray.shape[1]):
            if j < xMin or j > xMax or i < yMin or i > yMax:
                backgroundList.append(beamArray[i,j])

    averageBackground = np.mean(backgroundList) #Average the background
    beamMinusBackground = beamArray - averageBackground #Substract the average background
    beamMinusBackground[beamMinusBackground < 0] = 0 #Set -ve values to zero
    beamMinusBackground = np.around(beamMinusBackground) #round values to nearest integer
    beamMinusBackground = beamMinusBackground.astype(int) #Store values as integers
    return beamMinusBackground
