import base
import numpy as np
import sys
import subprocess
import pickle
import os
import tempfile

PROCESSDIR = os.path.join(os.path.dirname(__file__), "process_data/")

def get_pos_midpoint(r1, r2, box):
    """
    return the midpoint of two vectors r1 and r2
    
    use the minimum image: i.e. make sure we get a sensible answer if the                                                                     
    positions r1 and r2 correspond to nucleotides which were put at opposite                                                                  
    ends of the box due to pbc's                                                                                                            
    """
    assert (isinstance (box, np.ndarray) and len(box) == 3)
    return r1 - min_distance(r1, r2, box)/2

def array_index(mylist, myarray):
    """
    Return the index in a list mylist of a numpy array myarray
    """
    return map(lambda x: (myarray == x).all(), mylist).index(True)

def min_distance (r1, r2, box):
    """
    return the minimum image distance in going from r1 to r2, in a box of size box

    stolen from base.py Nucleotide.distance()
    """
    assert (isinstance (box, np.ndarray) and len(box) == 3)
    dr = r2 - r1
    dr -= box * np.rint (dr / box)
    return dr

def vecs2spline(vecs, per):
    import scipy.interpolate
    # interpolate vecs by interpolating each cartesian co-ordinate in turn
    xx = [vec[0] for vec in vecs]
    yy = [vec[1] for vec in vecs]
    zz = [vec[2] for vec in vecs]
    # NB s = 0 forces interpolation through all data points
    spline_xx = scipy.interpolate.splrep(range(len(xx)), xx, k = 3, s = 0, per = per)
    spline_yy = scipy.interpolate.splrep(range(len(yy)), yy, k = 3, s = 0, per = per)
    spline_zz = scipy.interpolate.splrep(range(len(zz)), zz, k = 3, s = 0, per = per)
    return [spline_xx, spline_yy, spline_zz], (0, len(xx)-1)

def get_base_spline(strand, reverse = False):
    """
    return a cartesian spline that represents a fit through the bases for the strand 'strand'

    args:
    strand: base.Strand object
    """
    import scipy

    base_pos = []
    for nuc in strand._nucleotides:
        base_pos.append(nuc.get_pos_base())

    if reverse:
        base_pos.reverse()

    if strand._circular:
        if reverse:
            base_pos.append(strand._nucleotides[-1].get_pos_base())
        else:
            base_pos.append(strand._nucleotides[0].get_pos_base())

    # interpolate bbms by interpolating each cartesian co-ordinate in turn
    xx = [vec[0] for vec in base_pos]
    yy = [vec[1] for vec in base_pos]
    zz = [vec[2] for vec in base_pos]
    # NB s = 0 forces interpolation through all data points
    spline_xx = scipy.interpolate.splrep(range(len(xx)), xx, k = 3, s = 0, per = strand._circular)
    spline_yy = scipy.interpolate.splrep(range(len(yy)), yy, k = 3, s = 0, per = strand._circular)
    spline_zz = scipy.interpolate.splrep(range(len(zz)), zz, k = 3, s = 0, per = strand._circular)
    return [spline_xx, spline_yy, spline_zz], (0, len(xx)-1)

def get_bb_spline(strand, reverse = False):
    """
    return a cartesian spline that represents a fit through the backbone of a duplex

    args:
    strand: base.Strand object
    """
    import scipy

    bb_pos = []
    for nuc in strand._nucleotides:
        bb_pos.append(nuc.get_pos_back())
    if reverse:
        bb_pos.reverse()
    if strand._circular:
        if reverse:
            bb_pos.append(strand._nucleotides[-1].get_pos_back())
        else:
            bb_pos.append(strand._nucleotides[0].get_pos_back())
 
    # interpolate bbms by interpolating each cartesian co-ordinate in turn
    xx = [vec[0] for vec in bb_pos]
    yy = [vec[1] for vec in bb_pos]
    zz = [vec[2] for vec in bb_pos]
    # NB s = 0 forces interpolation through all data points
    spline_xx = scipy.interpolate.splrep(range(len(xx)), xx, k = 3, s = 0, per = strand._circular)
    spline_yy = scipy.interpolate.splrep(range(len(yy)), yy, k = 3, s = 0, per = strand._circular)
    spline_zz = scipy.interpolate.splrep(range(len(zz)), zz, k = 3, s = 0, per = strand._circular)
    return [spline_xx, spline_yy, spline_zz], (0, len(xx)-1)

def get_sayar_twist(spline1, spline2, smin, smax, npoints = 1000, circular = False, integral_type = "simple"):
    """
    return the twist for a given pair of spline fits, one through the backbone of each duplex

    from Sayar et al. 2010 Phys. Rev. E

    Just need to integrate along the contour parameter s that is common to both splines. We need the normalised tangent vector to the spline formed by the midpoint of the two input splines t(s), the normalised normal vector formed by the vectors between the splines u(s), and the derivative of the normalised normal vector between the splines d/ds (u(s)). NB, the normal u(s) vector should be orthogonal to the tangent vector t(s); we ensure this by using only the component orthogonal to t(s).

    Using integral_type = 'simple' and npoints = 200, it will give a correct twist, or at least one that gives a conserved linking number when combined with get_sayar_writhe

    args:
    spline1: list of 3 splines corresponding to x, y and z spline through strand 1's backbone
    spline2: list of 3 splines corresponding to x, y and z spline through strand 2's backbone -- NB the splines should run in the same direction, i.e. one must reverse one of the splines if they come from get_base_spline (e.g. use get_base_spline(reverse = True))
    smin: minimum value for s, which parameterises the splines
    smax: maximum value for s, which parameterises the splines
    npoints: number of points for the discrete integration
    """

    import scipy.interpolate
    import scipy.integrate
    
    s1xx, s1yy, s1zz = spline1[0]
    s2xx, s2yy, s2zz = spline2[0]

    # bpi is the base pair index parameter that common to both splines
    bpi = np.linspace(smin, smax, npoints)

    # find the midpoint between the input splines, as a function of base pair index
    mxx = (scipy.interpolate.splev(bpi, s1xx) + scipy.interpolate.splev(bpi, s2xx)) / 2
    myy = (scipy.interpolate.splev(bpi, s1yy) + scipy.interpolate.splev(bpi, s2yy)) / 2
    mzz = (scipy.interpolate.splev(bpi, s1zz) + scipy.interpolate.splev(bpi, s2zz)) / 2
    

#    mxx = scipy.interpolate.splev(bpi, s1xx) 
#    myy = scipy.interpolate.splev(bpi, s1yy) 
#    mzz = scipy.interpolate.splev(bpi, s1zz) 

    # contour_len[ii] is contour length along the midpoint curve of point ii
    delta_s = [np.sqrt((mxx[ii+1]-mxx[ii])**2+(myy[ii+1]-myy[ii])**2+(mzz[ii+1]-mzz[ii])**2) for ii in range(len(bpi)-1)]
    contour_len = np.cumsum(delta_s)
    contour_len = np.insert(contour_len, 0, 0)

    # ss is a linear sequence from first contour length element (which is 0) to last contour length element inclusive
    ss = np.linspace(contour_len[0], contour_len[-1], npoints)

    # get the midpoint spline as a function of contour length
    msxx = scipy.interpolate.splrep(contour_len, mxx, k = 3, s = 0, per = circular)
    msyy = scipy.interpolate.splrep(contour_len, myy, k = 3, s = 0, per = circular)
    mszz = scipy.interpolate.splrep(contour_len, mzz, k = 3, s = 0, per = circular)

    # find the tangent of the midpoint spline. 
    # the tangent t(s) is d/ds [r(s)], where r(s) = (mxx(s), myy(s), mzz(s)). So the tangent is t(s) = d/ds [r(s)] = (d/ds [mxx(s)], d/ds [myy(s)], d/ds [mzz(s)])
    # get discrete array of normalised tangent vectors; __call__(xxx, 1) returns the first derivative
    # the tangent vector is a unit vector
    dmxx = scipy.interpolate.splev(ss, msxx, 1)
    dmyy = scipy.interpolate.splev(ss, msyy, 1)
    dmzz = scipy.interpolate.splev(ss, mszz, 1)
    
    def tsmooth(s):
        return np.array([dmxx[s], dmyy[s], dmzz[s]])
    
    tt = list(range(len(ss)))
    for ii in range(len(ss)):
        tt[ii] = np.array([dmxx[ii], dmyy[ii], dmzz[ii]])
  
    
    # get the normal vector via n(s) = dt(s)/ds = d^2/ds^2[r(s)]    
    ddmxx = scipy.interpolate.splev(ss, msxx, 2)
    ddmyy = scipy.interpolate.splev(ss, msyy, 2)
    ddmzz = scipy.interpolate.splev(ss, mszz, 2)
    nn = list(range(len(ss)))
    for ii in range(len(ss)):
        nn[ii] = np.array([ddmxx[ii], ddmyy[ii], ddmzz[ii]])
        
    def nsmooth(s):
        return np.array(ddmxx[s], ddmxx[s], ddmxx[s])
        
        
    # get the binormal vector via cross product
    bb = list(range(len(ss)))
    for ii in range(len(ss)):
        bb[ii] = np.cross(tt[ii],nn[ii])
    

    # we also need the 'normal' vector u(s) which points between the base pairs. (or between the spline fits through the backbones in this case)
    # n.b. these uxx, uyy, uzz are not normalised
    uxx_bpi = scipy.interpolate.splev(bpi, s2xx) - scipy.interpolate.splev(bpi, s1xx)
    uyy_bpi = scipy.interpolate.splev(bpi, s2yy) - scipy.interpolate.splev(bpi, s1yy)
    uzz_bpi = scipy.interpolate.splev(bpi, s2zz) - scipy.interpolate.splev(bpi, s1zz)

    # get the normal vector spline as a function of contour length
    suxx = scipy.interpolate.splrep(contour_len, uxx_bpi, k = 3, s = 0, per = circular)
    suyy = scipy.interpolate.splrep(contour_len, uyy_bpi, k = 3, s = 0, per = circular)
    suzz = scipy.interpolate.splrep(contour_len, uzz_bpi, k = 3, s = 0, per = circular)

    # evaluate the normal vector spline as a function of contour length
    uxx = scipy.interpolate.splev(ss, suxx)
    uyy = scipy.interpolate.splev(ss, suyy)
    uzz = scipy.interpolate.splev(ss, suzz)

    uu = list(range(len(ss)))
    for ii in list(range(len(ss))):
        uu[ii] = np.array([uxx[ii], uyy[ii], uzz[ii]])
        uu[ii] = uu[ii] - np.dot(tt[ii], uu[ii]) * tt[ii]
        # the normal vector should be normalised
        uu[ii] = norm(uu[ii])
        
    def usmooth(s):
        usmooth = np.array([uxx[s], uyy[s], uzz[s]])
        usmooth = usmooth - np.dot(tt[s], uu[s]) * tt[s]
        return usmooth

    # and finally we need the derivatives of that vector u(s). It takes a bit of work to get a spline of the normalised version of u from the unnormalised one
    nuxx = [vec[0] for vec in uu]
    nuyy = [vec[1] for vec in uu]
    nuzz = [vec[2] for vec in uu]
    nusxx = scipy.interpolate.splrep(ss, nuxx, k = 3, s = 0, per = circular)
    nusyy = scipy.interpolate.splrep(ss, nuyy, k = 3, s = 0, per = circular)
    nuszz = scipy.interpolate.splrep(ss, nuzz, k = 3, s = 0, per = circular)
    duxx = scipy.interpolate.splev(ss, nusxx, 1)
    duyy = scipy.interpolate.splev(ss, nusyy, 1)
    duzz = scipy.interpolate.splev(ss, nuszz, 1)
    duu = list(range(len(ss)))
    for ii in list(range(len(ss))):
        duu[ii] = np.array([duxx[ii], duyy[ii], duzz[ii]])
        
    def dusmooth(s):
        return np.array([duxx[s], duyy[s], duzz[s]])
    
    def twist_integrand(s):
        return np.dot(tsmooth(s), np.cross(usmooth(s), dusmooth(s)))
    
    ds = float(contour_len[-1] - contour_len[0]) / (npoints - 1)
    # do the integration w.r.t. s
    if circular:
        srange = range(len(ss)-1)
    else:
        srange = range(len(ss))
    if integral_type == "simple":
        integral = 0
        for ii in srange:
            triple_scalar_product = np.dot(tt[ii], np.cross(uu[ii], duu[ii]))
            integral += triple_scalar_product * ds
    elif integral_type == "quad":
        assert False, "not currently supported; shouldn't be difficult to implement if wanted"
        integral, err = scipy.integrate.quad(twist_integrand, ss[0], ss[-1], args = (msxx, msyy, mszz, nusxx, nusyy, nuszz), limit = 500)
        print(sys.stderr, "error estimate:", err)
        
    twist = integral/(2 * np.pi)

    return twist

def get_sayar_writhe(splines1, smin, smax, splines2 = False, npoints = 1000, debug = False, circular = False, integral_type = "simple"):
    """
    return the writhe for a 3D spline fit through a set of duplex midpoints

    from Sayar et al. 2010 Phys. Rev. E

    Using integral_type = 'simple' and npoints = 200, it will give a correct writhe, or at least one that gives a conserved linking number when combined with get_sayar_twist

    args:
    splines1: list of 3 splines corresponding to either (if not splines2) a 3D spline through the duplex or (if splines2) strand 1's bases
    smin: minimum value for s, which parameterises the splines
    smax: maximum value for s, which parameterises the splines
    splines2: optionally, (see splines1) list of 3 splines corresponding to a 3D spline through strand2's bases
    npoints: number of points for the discrete integration
    debug: print a load of debugging information
    """
    import scipy.integrate

    # bpi is the base pair index parameter that common to both duplex's splines
    bpi = np.linspace(smin, smax, npoints)

    ## get midpoint splines sxx, syy, szz
    if not splines2:
        # splines1 is the midpoint 3D spline as a function of base pair index
        sxx_bpi, syy_bpi, szz_bpi = splines1[0]
        xx_bpi = scipy.interpolate.splev(bpi, sxx_bpi)
        yy_bpi = scipy.interpolate.splev(bpi, syy_bpi)
        zz_bpi = scipy.interpolate.splev(bpi, szz_bpi)
    else:
        # take splines1 and splines2 to be the splines through the bases of each backbone; in that case we need to find the midpoint here first
        s1xx_bpi, s1yy_bpi, s1zz_bpi = splines1[0]
        s2xx_bpi, s2yy_bpi, s2zz_bpi = splines2[0]

        # find the midpoint as a function of base pair index between the input splines 
        xx_bpi = (scipy.interpolate.splev(bpi, s1xx_bpi) + scipy.interpolate.splev(bpi, s2xx_bpi)) / 2
        yy_bpi = (scipy.interpolate.splev(bpi, s1yy_bpi) + scipy.interpolate.splev(bpi, s2yy_bpi)) / 2
        zz_bpi = (scipy.interpolate.splev(bpi, s1zz_bpi) + scipy.interpolate.splev(bpi, s2zz_bpi)) / 2

    # contour_len[ii] is contour length along the midpoint curve of point ii
    delta_s = [np.sqrt((xx_bpi[ii+1]-xx_bpi[ii])**2+(yy_bpi[ii+1]-yy_bpi[ii])**2+(zz_bpi[ii+1]-zz_bpi[ii])**2) for ii in range(len(bpi)-1)]
    contour_len = np.cumsum(delta_s)
    contour_len = np.insert(contour_len, 0, 0)

    # ss is a linear sequence from first contour length element (which is 0) to last contour length element inclusive
    ss = np.linspace(contour_len[0], contour_len[-1], npoints)

    sxx = scipy.interpolate.splrep(contour_len, xx_bpi, k = 3, s = 0, per = circular)
    syy = scipy.interpolate.splrep(contour_len, yy_bpi, k = 3, s = 0, per = circular)
    szz = scipy.interpolate.splrep(contour_len, zz_bpi, k = 3, s = 0, per = circular)
    xx = scipy.interpolate.splev(ss, sxx)
    yy = scipy.interpolate.splev(ss, syy)
    zz = scipy.interpolate.splev(ss, szz)

    # find the tangent of the midpoint spline. 
    # the tangent t(s) is d/ds [r(s)], where r(s) = (mxx(s), myy(s), mzz(s)). So the tangent is t(s) = d/ds [r(s)] = (d/ds [mxx(s)], d/ds [myy(s)], d/ds [mzz(s)])
    # get discrete array of tangent vectors; __call__(xxx, 1) returns the first derivative
    dxx = scipy.interpolate.splev(ss, sxx, 1)
    dyy = scipy.interpolate.splev(ss, syy, 1)
    dzz = scipy.interpolate.splev(ss, szz, 1)
    tt = list(range(len(ss)))
    for ii in list(range(len(ss))):
        tt[ii] = np.array([dxx[ii], dyy[ii], dzz[ii]])
    
    # do the double integration w.r.t. s and s'
    if integral_type == "simple":
        integral = 0
        if circular:
            srange = list(range(len(ss)-1))
            ds = float(contour_len[-1] - contour_len[0]) / (npoints - 1)
        else:
            srange = list(range(len(ss)))
            ds = float(contour_len[-1] - contour_len[0]) / npoints
        for ii in srange:
            for jj in srange:
                # skip ii=jj and use symmetry in {ii, jj}
                if ii > jj:
                    diff = np.array([xx[ii]-xx[jj], yy[ii] - yy[jj], zz[ii] - zz[jj]])
                    diff_mag = np.sqrt(np.dot(diff, diff))
                    diff_frac = diff / (diff_mag ** 3)
                    triple_scalar_product = np.dot(np.cross(tt[ii], tt[jj]), diff_frac)
                    integral += triple_scalar_product * ds * ds
        # multiply by 2 because we exploited the symmetry in {ii, jj} to skip half of the integral
        integral *= 2
    elif  integral_type == "dblquad":
        # contour_len[0] to contour[-1] SHOULD be a complete integral on the closed curve; i.e. sxx(contour_len[0]) = sxx(contour_len[-1]) etc.
        val, err = scipy.integrate.dblquad(writhe_integrand, ss[0], ss[-1], lambda x: ss[0], lambda x: ss[-1], args = (sxx, syy, szz, ss[-1]))
        print(sys.stderr, err)
        integral = val
    elif integral_type == "chopped dblquad":
        integral = 0
        for ss_coarse in np.linspace(ss[0], ss[-1], 10):
            for ss_coarse_prime in np.linspace(ss[0], ss[-1], 10):
                val, err = scipy.integrate.dblquad(writhe_integrand, ss_coarse, ss_coarse + float(ss[-1]-ss[0])/9, lambda x: ss_coarse_prime, lambda x: ss_coarse_prime + float(ss[-1]-ss[0])/9, args = (sxx, syy, szz, contour_len[-1]))
                print(err)
                integral += val
    elif integral_type == "quad":
        integral, err = scipy.integrate.quad(writhe_integrand2, ss[0], ss[-1], args = (sxx, syy, szz, ss[0], ss[-1]), limit = 100, epsabs = 1e-5, epsrel = 0)
    elif integral_type == "simps":
        srange = range(len(ss))
        integrand = [[] for ii in srange]
        for ii in srange:
            for jj in srange:
                # skip ii=jj
                if ii == jj:
                    triple_scalar_product = 0
                else:
                    diff = np.array([xx[ii]-xx[jj], yy[ii] - yy[jj], zz[ii] - zz[jj]])
                    diff_mag = np.sqrt(np.dot(diff, diff))
                    diff_frac = diff / (diff_mag ** 3)
                    triple_scalar_product = np.dot(np.cross(tt[ii], tt[jj]), diff_frac)
                integrand[ii].append(triple_scalar_product)
        integral = scipy.integrate.simps(scipy.integrate.simps(integrand, ss), ss)
    else:
        assert False
        
    writhe = float(integral) / (4*np.pi)

    return writhe
    
def get_vhelix_vis(fname):
    """
    read the 'virtual helix visibility' from a text file and ignore any virtual helices set to invisible
    the syntax is the same as for the base.py visibility
    modified from the base.py strand visibility parsing
    """
    actions = {'vis' : True, 'inv' : False}
    visibility_list = []
    path = fname

    try:
        inp = open (path, 'r')
    except:
        base.Logger.log ("Origami visibility file `" + path + "' not found. Assuming default visibility", base.Logger.INFO)
        return False

    base.Logger.log("Using origami visibility file " +path, base.Logger.INFO)
    # read the visibility from a file; if we got here the file was opened
    lines = []
    for line in inp.readlines():
        line = line.strip().lower()
        # remove everything that comes after '#'
        line = line.split('#')[0]
        if len(line) > 0: lines.append(line)

    for line in lines:
        if '=' in line:
            sp = line.split("=")
            for p in sp[0]:
                one, two, three = [p.strip(), '=', sp[1]]
        else:
            one, two, three = line, "", ""

        if two != '=':
            base.Logger.log ("Lines in visibility must begin with one of inv=, vis= and default=. Skipping this line: --" + line + "--", base.Logger.WARNING)
            continue

        if one == 'default':
            if three not in ['inv', 'vis']:
                base.Logger.log ("Wrong default in visibility file. Assuming visible as default", base.Logger.WARNING)
                three = 'vis'
            if three == 'inv':
                mode = 'inv'
            else:
                mode = 'vis'
        else:
            # filter removes all the empty strings
            arr = [a.strip() for a in filter(None, three.split(','))]
            for a in arr:
                try:
                    ind = int(a)
                except:
                    base.Logger.log ("Could not cast '%s' to int. Assuming 0" % a, base.Logger.WARNING)
                    ind = 0
                try:
                    visibility_list.append(ind)
                except:
                    base.Logger.log ("vhelix %i does not exist in system, cannot assign visibility. Ignoring" % ind, base.Logger.WARNING)

    return mode, visibility_list

def norm(vec):
    return vec / np.sqrt(np.dot(vec,vec))

def parse_scanfile(infile):
    try:
        f = open(infile, "r")
    except IOError:
        base.Logger.log("could not open file %s" % infile, base.Logger.CRITICAL)
        sys.exit()
    for line in f.readlines():
        if line.startswith('begin_vb'):
            begin_vb = int(line.split()[2])
        elif line.startswith('end_vb'):
            end_vb = int(line.split()[2])
        elif line.startswith('vh_list'):
            vh_list = [int(x) for x in (line.split()[2]).split(",")]
        elif line.startswith('trim'):
            trim = int(line.split()[2])
    try:
        return begin_vb, end_vb, vh_list, trim
    except NameError:
        base.Logger.log("error while reading scan file %s, dying now" % infile, base.Logger.CRITICAL)
        sys.exit()

def get_scaffold_index(system):
    # find the index of the scaffold strand in the system
    strand_lengths = [strand.get_length() for strand in system._strands]
    return strand_lengths.index(max(strand_lengths))

def get_bb_midpoint(system, strand, n_index, interaction_list):
    base.Logger.log("origami_utils.get_bb_midpoint: deprecated function, use origami_utils.Origami.get_bb_midpoint", base.Logger.WARNING)
    # get midpoint vector between 2 hybridised bases
    r1 = strand._nucleotides[n_index].get_pos_base()
    r2 = system._nucleotides[interaction_list[strand._nucleotides[n_index].index]].get_pos_base()
    vec = (r1+r2)/2
    return vec

def get_nucleotide(vhelix, vbase, vhelix_indices):
    # find the system nucleotide index of a nucleotide given a position on the origami
    if vhelix % 2 == 0:
        dir = 1
    else:
        dir = -1
    return vhelix_indices[vhelix] + vbase * dir

def parse_vh_data(filename, origami):
    """
    get vhelices data from file - format is either <auto,[number of vhelices]> or, if a region of an origami is to be analysed, <[region width],[list of starting nucleotide index for the region for each vhelix]>
    """
    scaf_index = get_scaffold_index(origami._sys)
    vhelix_def_file = open(filename, "r")
    data = [x for x in vhelix_def_file.readline().replace("\n","").replace(" ","").split(",")]
    if data[0] == "auto":
        origami.num_vh = int(data[1])
        origami.width = origami._sys._strands[scaf_index].get_length() / origami.num_vh
        origami.vhelix_indices = []
        start_nuc_ind = -1
        for i in range(origami.num_vh):
            if i % 2 == 0:
                start_nuc_ind += 1
            else:
                start_nuc_ind += origami.width*2 - 1
            origami.vhelix_indices.append(start_nuc_ind)
    else:
        origami.width = int(data[0])
        origami.vhelix_indices = [int(x) for x in data[1:]]
        origami.num_vh = len(origami.vhelix_indices)
    base.Logger.log("using file data.vhd, %d virtual helices found" % origami.num_vh, base.Logger.INFO)

def print_arrow_debug_line(begin, end, file):
    if file:
        file.write("draw arrow {%f %f %f} {%f %f %f}\n" % (begin[0], begin[1], begin[2], end[0], end[1], end[2]))
    return 0

def open_arrow_debug_file(filename, type="w"):
    f_arrow = open(filename, type)
    f_arrow.write(
        "color Display Background white\n" +
        "set mName [mol new]\n" +
        "proc vmd_draw_arrow {mol start end} {\n" +
        "# an arrow is made of a cylinder and a cone\n"
        "set middle [vecadd $start [vecscale 0.65 [vecsub $end $start]]]\n" +
        "graphics $mol cylinder $start $middle radius 0.05\n" +
        "graphics $mol cone $middle $end radius 0.15\n" +
        "}\n")
    return f_arrow

def print_box_debug_line(vecs, file):
    if file:
        if len(vecs) == 4:
            file.write("draw box {%f %f %f} {%f %f %f} {%f %f %f} {%f %f %f}\n" % (vecs[0][0], vecs[0][1], vecs[0][2], vecs[1][0], vecs[1][1], vecs[1][2], vecs[2][0], vecs[2][1], vecs[2][2], vecs[3][0], vecs[3][1], vecs[3][2]))
        else:
            base.Logger.log("drawing boxes only works for 4 vertices at the moment", base.Logger.WARNING)

def open_box_debug_file(filename):
    f_boxes = open(filename, "w")
    f_boxes.write(
        "color Display Background white\n" +
        "set mName [mol new]\n" +
        "proc vmd_draw_box {mol vert1 vert2 vert3 vert4} {\n" +
        "# a 'box' is a plane made of 2 triangles here\n" +
        "graphics $mol triangle $vert1 $vert2 $vert3\n" +
        "graphics $mol triangle $vert1 $vert4 $vert3\n" +
        "}\n")
    return f_boxes

def angle_sense(v1, v2, axis):
    # return the angle between two vectors, using a 3rd vector to determine the sense of the angle
    v1 = norm(v1)
    v2 = norm(v2)
    axis = norm(axis)
    angle = np.arccos(np.dot(v1,v2))
    if np.dot(norm(np.cross(v1,v2)),axis) < 0:
        angle *= -1 # attempt to check the 'sense' of the angle w.r.t. an axis
    return angle

def dihedral_angle_sense(v1, v2, axis, sanity_check = False):
    # return the dihedral angle between two vectors, using a 3rd vector to determine the sense of the angle and to define the axis normal to the plane we want the vectors in
    v1n = norm(v1)
    v2n = norm(v2)
    axis = norm(axis)
    if sanity_check and (abs(np.dot(v1n,axis)) > 0.6 or abs(np.dot(v2n,axis)) > 0.6):
        return False
    # in plane
    v1p = v1n - np.dot(v1n,axis)*axis
    v2p = v2n - np.dot(v2n,axis)*axis
    v1p = norm(v1p)
    v2p = norm(v2p)
    angle = np.arccos(np.dot(v1p,v2p))
    if np.dot(norm(np.cross(v1p,v2p)),axis) < 0:
        angle *= -1 # attempt to check the 'sense' of the angle w.r.t. an axis
    return angle

class vhelix_vbase_to_nucleotide(object):
    # at the moment squares with skips in have entries in the dicts but with the nucleotide list empty (rather than having no entry) - I'm not sure whether or not this is desirable. It's probably ok
    def __init__(self):
        self._scaf = {}
        self._stap = {}
        self.nuc_count = 0 # record the nucleotide count, updated only after a whole strand is added
        self.strand_count = 0

    def add_scaf(self, vh, vb, strand, nuc):
        self._scaf[(vh, vb)] = (strand, nuc)

    def add_stap(self, vh, vb, strand, nuc):
        self._stap[(vh, vb)] = (strand, nuc)

    # these methods use a reference vhvb2n object to make the final vhvb2n object
    def add_scaf_strand(self, add_strand, reference, continue_join = False):
        count = 0
        size = len(self._scaf)
        for (vh, vb), [strand_ind, nuc] in reference._scaf.iteritems():
            if strand_ind == add_strand:
                self.add_scaf(vh, vb, self.strand_count, [x + self.nuc_count for x in nuc])
                count += len(nuc)
        self.nuc_count += count
        if len(self._scaf) == size:
            return 1
        else:
            if continue_join == False:
                self.strand_count += 1
            return 0

    def add_stap_strand(self, add_strand, reference, continue_join = False):
        count = 0
        size = len(self._stap)
        for (vh, vb), [strand_ind, nuc] in reference._stap.iteritems():
            if strand_ind == add_strand:
                self.add_stap(vh, vb, self.strand_count, [x + self.nuc_count for x in nuc])
                count += len(nuc)
        self.nuc_count += count
        if len(self._stap) == size:
            return 1
        else:
            if continue_join == False:
                self.strand_count += 1
            return 0

    def add_strand(self, add_strand, reference, continue_join = False):
        if self.add_scaf_strand(add_strand, reference, continue_join) and self.add_stap_strand(add_strand, reference, continue_join):
            #base.Logger.log("while adding strand %s to vhelix_vbase_to_nucleotide object; either strand already present or strand not found in reference object" % add_strand, base.Logger.WARNING)
            return 1
        else:
            return 0

