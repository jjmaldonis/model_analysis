
import sys
import math
import numpy as np
from model import Model
import scipy
from scipy import *
import scipy.signal
import scipy.ndimage
import scipy.spatial
rot3d = __import__('3d_rot')
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion


def dist2d(x1,y1,x2,y2):
    return math.sqrt((x1-x2)**2+(y1-y2)**2)


def detect_peaks(image):
    """
    Takes an image and detect the peaks usingthe local maximum filter.
    Returns a boolean mask of the peaks (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """

    # define an 8-connected neighborhood
    neighborhood = generate_binary_structure(2,2)

    #apply the local maximum filter; all pixel of maximal value 
    #in their neighborhood are set to 1
    local_max = maximum_filter(image, footprint=neighborhood)==image
    #local_max is a mask that contains the peaks we are 
    #looking for, but also the background.
    #In order to isolate the peaks we must remove the background from the mask.

    #we create the mask of the background
    background = (image==0)

    #a little technicality: we must erode the background in order to 
    #successfully subtract it form local_max, otherwise a line will 
    #appear along the background border (artifact of the local maximum filter)
    eroded_background = binary_erosion(background, structure=neighborhood, border_value=1)

    #we obtain the final mask, containing only peaks, 
    #by removing the background from the local_max mask
    detected_peaks = local_max - eroded_background

    return detected_peaks


def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = scipy.mgrid[-size:size+1, -sizey:sizey+1]
    g = scipy.exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None) :
    """ smooths the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = scipy.signal.convolve2d(im,g, mode='valid')
    return(improc)


def rxrydist(atomi,atomj,m):
    x = (atomj.coord[0] - atomi.coord[0])
    y = (atomj.coord[1] - atomi.coord[1])
    while(x > m.lx/2): x = -m.lx + x
    while(y > m.ly/2): y = -m.ly + y
    while(x < -m.lx/2): x = m.lx + x
    while(y < -m.ly/2): y = m.ly + y
    return x,y


def rdf_2d(m,dr):
    mean = 0.0
    std = 0.0
    # dr is the bin size in atom coord units (A probably)
    # it should be large enough so that the intensity isn't 1 in every bin
    # but small enough to not overlook important information (ie peaks)

    # Algorithm:
    # Calculation the distance between all pairs of atoms, and calculate
    # the rx and ry components of them. This will results in two huge
    # lists: rx and ry. Then create a square matrix of size boxlen/dr,
    # and go through rx and ry simultaneously, incrementing the matrix
    # index if an rx,ry pair falls in that spot. Contour plot the matrix.

    size = int(np.ceil(m.lx/dr))
    #print("Initializing matrix: {0}x{1}".format(size,size))
    hist2d = np.zeros((size,size),dtype=np.int)
    rx = []
    ry = []
    bindict = {}

    #print("Calculating all distances...")
    for i,atomi in enumerate(m.atoms):
        for j,atomj in enumerate(m.atoms):
            if(j != i):
                x,y = rxrydist(atomi,atomj,m)
                xbin = int(math.floor((x+m.lx/2.0)/dr))
                ybin = int(math.floor((y+m.ly/2.0)/dr))
                #try:
                #    bindict[(xbin,ybin)].append( [i,j] )
                #except KeyError:
                #    bindict[(xbin,ybin)] = [[i,j]]
                #rx.append(x)
                #ry.append(y)
                hist2d[xbin][ybin] += 1
    # This print statement prints the raw histogram image.
    #print(hist2d.tolist())
    orig_hist = copy(hist2d)

    # Blur the image to get thing smoother for analysis.
    # hist2d is no longer usable.
    hist2d = blur_image(hist2d,10)
    #print(hist2d)

    # Note: scipy.ndimage.measurements.find_objects did not work well.

    # Use the stdev and mean to find out what value constitutes a "peak".
    std = np.std(hist2d)
    mean = np.mean(hist2d)
    sm = 3*std+mean
    lsm = np.array( [[0 if x < sm else x for x in list] for list in hist2d] )
    # This print statement prints the blurred image where the background is removed (set to 0).
    #print(lsm)

    # Detected peaks is a T/F mask. This actually finds where the peaks are.
    detected_peaks = detect_peaks(lsm)
    # Set the center of each peak to 0 for viewing purposes, hist2d is no longer usable.
    for i in range(0,len(detected_peaks)):
        for j in range(0,len(detected_peaks[i])):
            if(detected_peaks[i][j]): hist2d[i][j] = -1.0  # hist2d is blurred at this point as well
            if(detected_peaks[i][j]): orig_hist[i][j] = -1.0
    # This print line prints out the image matrix with each center black (ie 0).
    #print(hist2d.tolist())

    # Here we find where the peaks occur, using detected_peaks.
    peak_indexes = [[i*dr,j*dr] for i,ilist in enumerate(detected_peaks) for j,val in enumerate(detected_peaks[i]) if detected_peaks[i][j] == True]
    #print(peak_indexes)

    ## Calculate all the distances between peaks.
    peak_dists = scipy.spatial.distance.pdist( peak_indexes, 'euclidean')
    #peak_dists = [ x for x in peak_dists if x < 4.0 ]
    peak_dists.sort()
    #print(peak_dists)
    for x in peak_dists: print(x)

    return (orig_hist,hist2d)

    # Create a histogram of the distances between peaks.
    # These are the arrays you plot vs. each other in excel.
    # Take the weighted average before the cutoff to get the
    # correct plane spacing.
    # Not sure how to automate this.
    #peak_dist_hist = np.histogram(peak_dists,size*10.0)
    #print(peak_dist_hist[1].tolist())
    #print(peak_dist_hist[0].tolist())


    ## Get atoms that are in spots:
    ## Algorithm (OLD! Revise this):
    ## Create another matrix of size sizeXsize
    ## that will be used as a holder to identify all positions
    ## in the histogram that have intensities higher than
    ## sm = 3*stddev + mean. ie for every entry in mat that is
    ## greater than sm, increment the index in this new matrix
    ## as well as the indexes around that index. Then, go back
    ## through this matrix and set every index > 1 to True, and
    ## and everything else to False. These are the 'pixels' that
    ## we want to extract atoms from. Using GetAtomsInBin(...),
    ## put all the atoms in the True bins into a new model.
    ## Save that model!

    ## This says the spot size radius is 10 pixels.
    ## Need to look at image in Igor to figure this out.
    #tol = 10
    #tfmat = np.zeros((size,size))
    #for i in range(0,len(hist2d)):
    #    for j in range(0,len(hist2d[i])):
    #        kmin = i-tol
    #        kmax = i+tol
    #        if(kmin < 0): kmin = 0
    #        if(kmax > size-1): kmax = size-1
    #        for k in range(kmin,kmax):
    #            nmin = i-tol
    #            nmax = i+tol
    #            if(nmin < 0): nmin = 0
    #            if(nmax > size-1): nmax = size-1
    #            for n in range(nmin,nmax):
    #                if(hist2d[k][n] >= sm): tfmat[k][n] += 1
    #tfmat = [ [1 if x > 1 else 0 for x in ilist] for ilist in tfmat]

    #atoms = []
    #for i in range(0,size):
    #    for j in range(0,size):
    #        if(hist2d[i][j]):
    #            print(bindict[(i,j)])
    #            for k in range(0,len(bindict[(i,j)])):
    #                atoms.append(bindict[(i,j)][k][0])
    #                atoms.append(bindict[(i,j)][k][1])
    #atoms = list(set(atoms))
    #atoms = [m.atoms[i] for i in atoms]
    #new_model = Model("2D Hist extracted atoms\n", m.lx,m.ly,m.lz, atoms)

    #new_model.write_real_xyz('temp.xyz')
    #new_model.write_real_xyz()



def main():
    # Input file should be a (optionally flattened) model file rotated
    # in the correct orientation to see planes without streaking of 
    # the atom columns.
    # You should send the output to a file.
    # In that file, using vim, do :%s/\], \[/\r/g
    # which puts each matrix row on a new line,
    # then delete the leading [[ and ]] on the first and last lines.
    # This can be loaded into an excel file, from which the matrix
    # can then be pasted into a 2d wave of the correct column size
    # in Igor. The rows will be pasted in automatticaly. Then in 
    # the data browser, you can right click and click "New Image"!
    # From there, you can save the image by going into the Igor help
    # and searching for "save image". Click on the first of the two
    # links and then click on the button to #include <Image Saver>.
    # Go to Data > Save Waves > Save Image.
    # Set: Bit Depth = 16 bit, File Type = TIFF, and Path = _none_.
    # OR: ImageSave/D=16/T="TIFF" w as "w.tif" where 'w' is the name
    # of the wave you pasted into above. Now you can open ImageJ, set
    # the pixel size to 0.1 Ang (or whatever it is below), and use the
    # distance tool to measure the peak distances.
    m = Model(sys.argv[1])
    hist2d,blurred_hist2d = rdf_2d(m,0.1)
    temp = np.zeros((len(blurred_hist2d),len(blurred_hist2d[0])),dtype=np.int)
    for i,row in enumerate(blurred_hist2d):
        for j,x in enumerate(row):
            if(x == -1): temp[i][j] = 1
    np.savetxt('temp1.txt',temp)

    m = Model(sys.argv[2])
    hist2d,blurred_hist2d = rdf_2d(m,0.1)
    temp = np.zeros((len(blurred_hist2d),len(blurred_hist2d[0])),dtype=np.int)
    for i,row in enumerate(blurred_hist2d):
        for j,x in enumerate(row):
            if(x == -1): temp[i][j] = 2
    np.savetxt('temp2.txt',temp)

    #print(hist2d.tolist())
    #print(peak_dist_hist[0].tolist()) # histogram
    #print(peak_dist_hist[1].tolist()) # bin edges
    #np.savetxt('2d_rot_matrix.txt',hist2d.tolist())
    #np.savetxt('2d_rot_matrix_blurred.txt',blurred_hist2d.tolist())

    ## A first attempt at rotating the model and finding interesting things.
    #for t1 in [36]*7:
    #    ra = rot3d.calc_rot_array(t1,0,0)
    #    rot3d.rot(m,ra)
    #    rdf_2d(m,0.1)

if __name__ == '__main__':
    main()

