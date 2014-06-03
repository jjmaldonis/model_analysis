import numpy as np
import itertools
import scipy.spatial as spatial
import math
from collections import defaultdict


def main():
    # These are for t3 - full model - 0 rot
    #xx = np.array([-1.59661,-6.49117,1.31493,8.93178,-11.3095,-3.87922,3.95049,-8.6776,-1.14368,6.60032,-6.09572,1.5423,-11.1344,9.28251,-3.42097,4.29183,-8.44396,12.0538,-0.744052,7.02859,-13.6018,-5.75384,1.91543,-10.7154,9.777,-3.0811,4.6375,-8.09474,12.1967,-0.366933,7.26269,-12.995,-5.3918,2.3143,9.94702,-10.3263,-2.72239,4.97288,-7.71958,12.6045,-0.0339517,7.6468,-12.6997,-5.04385,2.65599,10.2525,-10.0009,-2.38506,5.3219,-7.32355,12.8661,0.301,8.02849,-12.2529,-4.69847,3.0119,10.6601,-9.85641,-1.97964,5.67557,13.5107,-7.10072,0.674167,8.36635,-12.1423,-4.38631,3.35134,11.1154,-9.3475,-1.61294,5.98824,-6.6541,1.07579,8.59859,-4.02672,3.78621,11.2388,-9.23234,-1.38635,6.49489,1.46454])
    #yy = np.array([-13.0135,-12.0045,-11.6794,-11.2696,-11.3734,-11.0096,-10.4323,-10.078,-9.71647,-9.38735,-8.98603,-8.5237,-8.07549,-8.11897,-7.71436,-7.32815,-6.90611,-6.94489,-6.5612,-6.15396,-6.18243,-5.73658,-5.33358,-4.88745,-4.80902,-4.53373,-4.111,-3.67341,-3.66839,-3.27151,-2.8533,-2.8538,-2.46533,-2.0914,-1.69791,-1.61275,-1.2895,-0.842847,-0.458095,-0.438408,-0.0354641,0.383534,0.385201,0.76331,1.2109,1.55113,1.65869,2.0126,2.39573,2.77632,2.75875,3.20249,3.60239,3.62001,4.03653,4.46168,4.85265,4.73056,5.27541,5.66213,6.13531,6.08266,6.48706,6.81082,6.94545,7.26832,7.64378,8.03617,8.08076,8.4479,8.78229,9.27676,9.74843,9.99324,10.3115,10.9273,11.2747,11.2236,11.5886,12.0599,12.7769])

    # These are for xtal.t3.allclusters.0rot.xyz
    #xx=np.array([-3.96954,-1.23944,-6.13183,1.52377,-3.50204,-8.37638,-0.844936,-5.76161,1.83726,-3.06934,4.59663,-8.04486,-0.365664,7.33512,-5.40676,2.29207,-2.6976,4.95053,-7.70217,-0.0376925,7.62148,-5.00266,2.63548,-2.36109,5.34406,-7.4045,0.302294,7.89081,-4.64835,2.99459,-1.93214,5.69076,0.752213,8.29594,3.40975,-1.67746,6.03777,1.21923,3.87908])
    #yy=np.array([-11.0532,-9.81174,-9.03512,-8.48442,-7.79604,-6.98633,-6.62881,-5.83473,-5.42329,-4.54319,-4.1367,-3.7085,-3.30234,-2.84677,-2.48796,-2.10459,-1.31195,-0.853427,-0.487057,-0.0365463,0.410394,0.799499,1.23164,2.02635,2.42725,2.78372,3.22103,3.53854,4.08515,4.45838,5.35181,5.74884,6.57751,7.03055,7.7151,8.33448,8.95794,10.0016,10.965])

    # These are for xtal.t1.cluster0.big.1rot.xyz
    #xx = np.array([ -13.3448, -11.1207, -9.06439, -12.0859, -6.88223, -9.81976, -4.44827, -7.67956, -2.47592, 2.81164, -5.2456, -0.0419648, 5.16167, -2.89557, 2.43396, 7.42777, 4.36434, 9.73583, 6.84026, 12.0019, 8.98047, 11.0787, 13.3028])
    #yy = np.array([ -9.31618, -7.63759, -6.37865, -5.45542, -4.90988, -3.8188, -3.31522, -2.35003, -1.67859, -1.09108, -0.671437, -0.0419648, 0.587507, 1.04912, 1.63663, 2.2661, 3.23129, 3.73487, 4.86792, 5.37149, 6.29472, 7.55366, 9.23225])

    # For JWH t3 xtal.all cluster
    #xx = np.array([ -5.07778, -0.419651, 2.05629, -9.69395, -7.30193, -4.86796, -2.51791, -0.167861, 2.22415, 4.61616, 6.96621, -9.52608, 9.31626, -7.218, -4.78403, -2.43398, -0.0419651, 2.39201, 4.74206, 7.17604, 9.48412, -9.35822, -7.05014, -4.7001, -2.30808, 0.0839303, 2.47594, 4.78403, 7.218, 9.61002, -2.14022, 0.377686, 4.99385])
    #yy = np.array([-6.16887, -5.5394, -5.11975, -4.19651, -3.7349, -3.35721, -3.02149, -2.72773, -2.35005, -2.05629, -1.8045, -1.38485, -1.42681, -1.00716, -0.755372, -0.377686, -0.0419651, 0.293756, 0.671442, 0.965198, 1.30092, 1.34288, 1.72057, 1.97236, 2.26612, 2.6438, 2.97952, 3.31525, 3.65097, 4.11258, 5.03582, 5.45547, 6.08494])

    # For al 3x3x3
    xx = np.array([ -2.73873, 0, 2.68893, -2.73873, 2.68893, 0, -2.73873, 0, 2.68893])
    yy = np.array([ -2.73873, -2.73873, -2.73873, 0, 0, 0, 2.68893, 2.68893, 2.68893])

    dists = spatial.distance.cdist([ [xx[i],yy[i]] for i in xrange(len(xx)) ],[[0,0]])
    dists = [ x[0] for x in dists ]
    sorted_dists = sorted(dists)
    #print(sorted_dists)
    center = dists.index(sorted_dists[0])
    #print(center,dists[center])
    dists = spatial.distance.cdist([ [xx[i],yy[i]] for i in xrange(len(xx)) ],[[xx[center],yy[center]]])
    dists = [ x[0] for x in dists ]
    sorted_dists = sorted(dists)
    print('Distances from center spot to other spots:')
    print(sorted_dists[1:])

    # y = m*x + b
    #A = np.vstack([xx, np.ones(len(xx))]).T
    #out = np.linalg.lstsq(A, yy)
    #m,b = out[0]
    #res = out[3]
    #print m,b,res

    symnum = 2 # number of degrees of symmetry about the 0 point
    starting_points = []
    i = 1
    while( len(starting_points) < symnum ):
        starting_points.append(dists.index(sorted_dists[i]))
        #print(starting_points[len(starting_points)-1],dists[starting_points[len(starting_points)-1]])
        i = i + 2
    print('Center: {0} @ ({1},{2})'.format(center,xx[center],yy[center]))
    #starting_points[2] = 13
    print('Starting points: {0}'.format(starting_points))

    # Go thru all the orientations (symnum) and find the corresponding lines
    # This only gives the lines thru the 0 point, I find the rest later
    point_array = []
    lines = []
    for point_index in starting_points:
        testx = np.array([ xx[center], xx[point_index] ])
        testy = np.array([ yy[center], yy[point_index] ])
        #print("test",testx,testy)
        #print(testx,testy)
        A = np.vstack([testx, np.ones(len(testx))]).T
        out = np.linalg.lstsq(A, testy)
        m,b = out[0]
        res = out[3]
        #print m,b,res
        # Search thru all the other points, looking for ones close to this line
        points = [center,point_index]
        for i in range(len(xx)):
            if(i != center and i != point_index):
                p = [xx[i],yy[i]] # point
                d = abs((m*xx[i] - yy[i] + b)/math.sqrt(m**2+1)) #Distance from p to line
                if( d < 1.0 ):
                    #print(p,d)
                    points.append(i) # accept if within 1.5 A
                    # Update line
                    np.array(list(testx).append(xx[i]))
                    np.array(list(testy).append(yy[i]))
                    A = np.vstack([testx, np.ones(len(testx))]).T
                    out = np.linalg.lstsq(A, testy)
                    m,b = out[0]
        #print(points)
        #if(len(points)): return
        point_array.append(sorted(points[:]))
        lines.append((m,b))
    point_array[1].append(7)
    print('\nOriginal lines and atoms in them: START')
    for i,l in enumerate(point_array):
        print(lines[i])
        for p in l:
            print('Index: {0}, coord: ({1},{2})'.format(p,xx[p],yy[p]))
        if(i < len(point_array)-1): print('')
    print('END\n')

    # Now that I have the line directions and all the points on those lines
    # that go thru the center point, the algorithm is to go thru each point
    # on a line and calculate the points on a DIFFERENT line that goes thru
    # that point.
    for j,line in enumerate(lines[0:symnum]): # for each line going thru the center
        for line2 in lines[0:symnum]:
            m = line2[0]
            #b = line[1]
            for point in point_array[j]: # for all points along the line
                if( point != center):
                    b = -m*xx[point] + yy[point]
                    #print([xx[point],yy[point]],line,b)
                    # Search thru all the other points, looking for ones close to this line
                    testx = np.array([ xx[point] ])
                    testy = np.array([ yy[point] ])
                    points = [point]
                    for i in range(len(xx)):
                        if(i not in point_array[j] ):
                            p = [xx[i],yy[i]] # point
                            d = abs((m*xx[i] - yy[i] + b)/math.sqrt(m**2+1)) #Distance from p to line
                            if( d < 1.0 ):
                                #print(p,d)
                                points.append(i) # accept if within 1.5 A
                                ## Update line
                                #np.array(list(testx).append(xx[i]))
                                #np.array(list(testy).append(yy[i]))
                                #A = np.vstack([testx, np.ones(len(testx))]).T
                                #out = np.linalg.lstsq(A, testy)
                                #m,b = out[0]
                    #print(points)
                    #if(len(points)): return
                    if(len(points) > 1 and sorted(points) not in point_array):
                        point_array.append(sorted(points[:]))
                        lines.append((m,b))
    #print('')
    #for i,l in enumerate(point_array):
    #    print(lines[i])
    #    for p in l:
    #        print('{0} {1}'.format(xx[p],yy[p]))
    #    print('')
    #print(len(lines))

    linecounts = defaultdict(int)
    # Count the number of lines in each direction:
    for line in lines:
        linecounts[line[0]] += 1
    linecounts = [ linecounts[x] for x in sorted(list(linecounts)) ]
    # linecounts is now in sorted order

    # Update all lines:
    for i,l in enumerate(point_array):
        testx = np.array([ xx[j] for j in l ])
        testy = np.array([ yy[j] for j in l ])
        A = np.vstack([testx, np.ones(len(testx))]).T
        out = np.linalg.lstsq(A, testy)
        m,b = out[0]
        lines[i] = (m,b)
        

    # Print waves
    # Need to copy into excel, delete the 0's
    # then past into igor and append on top of
    # the unblurred image
    print('Waves START')
    maxlen = 0
    minlen = 10000
    for line in point_array:
        if(len(line) > maxlen): maxlen = len(line)
        if(len(line) < minlen): minlen = len(line)
    for i in xrange(maxlen):
        string = ''
        for j,line in enumerate(point_array):
            try:
                m = lines[j][0]
                b = lines[j][1]
                y = m*xx[line[i]]+b
                string += str(xx[line[i]]) + ' ' + str(y) + ' ' # This line prints the line's values
                #string += str(xx[line[i]]) + ' ' + str(yy[line[i]]) + ' ' # This line prints the actual values
            except:
                string += '0 0 '
        string = string.strip()
        print(string)
    print('Waves END')

    ##for line in lines:
    #for line in sorted(lines, key=lambda tup: tup[0]):
    #    print(line)
    #print(len(lines))

    lines.sort(key=lambda tup: tup[0])
    plane_spacings = []
    start = 0
    for i in range(len(linecounts)): # ie for linecounts[i] in [11,9,11]
        avgm = 0.0
        # Go thru the selection of lines with similar slopes and avg them
        for j in range(linecounts[i]):
            avgm += lines[start+j][0]
        avgm /= linecounts[i]
        # Using the avg slope, set mperp. Also just set bperp to zero.
        mperp = -1.0/avgm
        bperp = 0.0
        #points = [] # DEBUGGING
        #for p in point_array[0]:
        #    x = xx[p]
        #    y = avgm*x
        #    points.append([x,y])
        #for point in points:
        #    print('{0} {1}'.format(point[0],point[1]))
        #print('')
        # Now go back thru all the lines and calculate the intersections
        # of the line with the perpendicular one we just found
        points = []
        for j in range(linecounts[i]):
            m = lines[start+j][0]
            b = lines[start+j][1]
            x = (bperp-b)/(m-mperp)
            y = m*x+b
            points.append([x,y])
            #print(lines[start+j],start)
        #for point in points:
        #    print('{0} {1}'.format(point[0],point[1]))
        #print('')
        # Calculate the distances between all the points, the first linecounts[i]-1
        # of these will be the plane spacings
        ds = spatial.distance.pdist(points)
        ds.sort()
        #print(ds)
        #print(ds[:linecounts[i]-1])
        plane_spacings.append((mperp,bperp))
        for d in ds[:linecounts[i]-1]:
            plane_spacings.append(d)
        start += linecounts[i]
            
    #plane_spacings.sort()
    print("Plane spacings!:")
    for d in plane_spacings:
        print(d)






if __name__ == '__main__':
    main()
