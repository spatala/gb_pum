import ovito.data as ovd
from ovito.data import CutoffNeighborFinder
import ovito.modifiers as ovm
import ovito.io as oio
import numpy as np

def analyze_gb_atoms(dir, dump_name):
    fname = dir+dump_name;
    data = compute_ovito_data(fname);
    non_p = identify_pbc(data);

    pts1 = np.array(data.particles['Position'][...]);

    lat_par = 4.05
    r_atom = lat_par/(2*np.sqrt(2))

    ##############################################################################
    ### Code to find the indices of grain boundary atoms
    gb_csm_cut = 1;
    csm1 = np.array(data.particles['c_csym'][...]);
    gb1_inds = np.where(np.array(csm1) > gb_csm_cut)[0];

    ##############################################################################
    ### Get all atoms within cut-off distance of GB atoms
    rCut = 3*lat_par;
    finder1 = CutoffNeighborFinder(rCut, data);

    nn_tot_inds = []
    for index in gb1_inds:
        nn_inds = []
        for neigh in finder1.find(index):
            nn_inds.append(neigh.index)
        nn_tot_inds.append(nn_inds)

    flat_list = [item for sublist in nn_tot_inds for item in sublist]
    nn_tot_inds = np.unique(np.array(flat_list))

    nn_gb_inds = 0*gb1_inds - 1;
    ct1 = 0;
    for gb1 in gb1_inds:
       nn_gb_inds[ct1]=np.where((nn_tot_inds - gb1) == 0)[0][0]
       ct1 = ct1 + 1;

    ##############################################################################
    #### Trim simulation data-points and reindex GB atoms
    pts = np.array(data.particles['Position'][...]);
    pts1 = pts[nn_tot_inds,:];
    gb_inds = np.copy(nn_gb_inds)

    pbc = data.cell.pbc
    pbc = np.asarray(pbc) + 0
    non_p = np.where(pbc == 0)[0][0]

    arr = np.array([0, 1, 2])
    arr = np.delete(arr, non_p)


    sim_cell = data.cell[...]
    sim_nonp_vec = np.array(sim_cell[:, non_p])
    sim_1vec = np.array(sim_cell[:, arr[0]])
    sim_2vec = np.array(sim_cell[:, arr[1]])
    sim_orig = np.array(sim_cell[:, 3])

    p1_vec = np.array([sim_1vec[arr[0]], sim_1vec[arr[1]]])
    p2_vec = np.array([sim_2vec[arr[0]], sim_2vec[arr[1]]])
    [n1, n2] = num_rep_2d(p1_vec, p2_vec, rCut)

    pts_w_imgs, inds_arr = create_imgs(pts1, n1, n2, sim_1vec, sim_2vec, non_p)

    #####################################################################################
    import scipy.io as sio
    to_plot = {}
    to_plot['pts'] = pts_w_imgs
    box_cell = np.array(data.cell[...])
    to_plot['box_cell'] = box_cell
    to_plot['gb_inds'] = gb_inds
    to_plot['uc_inds'] = inds_arr
    mat_name = dump_name + '_pts_box.mat'
    sio.savemat(mat_name, to_plot)
    #####################################################################################

def num_rep_2d(xvec, yvec, rCut):
    """
    Function finds the number of replications necessary such that thecircle of radius rCut at the
    center of the primitive-cell lies completely inside the super-cell.

    Parameters
    ------------
    xvec :
        The basis vector in x direction in x-z plane
    yvec :
        The basis vector in z direction in x-z plane
    rCut
        Cut-off radius for computing Delaunay triangulations

    Returns
    ------------
    [int(m_x), int(m_y)] :
        int(m_x) is the number of replications in x direction, int(m_y)
        is the number of replication in z direction.

    """
    c_vec_norm = np.linalg.norm(np.cross(xvec, yvec))
    d_y = c_vec_norm/(np.linalg.norm(yvec))
    d_x = c_vec_norm/(np.linalg.norm(xvec))
    m_x = np.ceil(rCut/d_y)
    m_y = np.ceil(rCut/d_x)

    return [int(m_x), int(m_y)]

def compute_ovito_data(filename0):
    """
    Computes the attributes of ovito

    Parameters
    ------------
    filename0 : string
        The name of the input file.

    Returns
    --------
    data : class
        all the attributes of data
    """
    pipeline = oio.import_file(filename0, sort_particles=True)
    dmod = ovm.PolyhedralTemplateMatchingModifier(rmsd_cutoff=.1)
    pipeline.modifiers.append(dmod)
    data = pipeline.compute()
    return data


def identify_pbc(data):
    """
    Function finds the non-periodic direction

    Parameters
    ------------
    data : class
        all the attributes of data

    Returns
    --------
    non_pbc : int
        The non-periodic direction. 0 , 1 or 2 which corresponds to
        x, y and z direction, respectively.
    """
    pbc = data.cell.pbc
    pbc = np.asarray(pbc) + 0
    non_pbc = np.where(pbc == 0)[0][0]
    return non_pbc


def create_imgs(pts1, n1, n2, sim_1vec, sim_2vec, non_p):
    """
    Creates the replicates of the main cell in X and Z direction.

    Parameters
    -------------
    pts1 :
        Indices of the atoms which Y value is in range [GBRegion[0] - rCut, GBRegion[1] + rCut].
    n1 :
        Number of replications in 1st periodic direction
    n2 :
        Number of replications in 2nd periodic direction
    sim_1vec :
        The simulation cell basis vector in 1st periodic direction
    sim_2vec :
        The simulation cell basis vector in 2nd periodic direction
    non_pbc : int
        The non-periodic direction. 0 , 1 or 2 which corresponds to
        x, y and z direction, respectively.

    Returns
    ----------
    pts_w_imgs :
        The position of atoms after replicating the box n_x and n_z times in X and Z direction.
    inds_array :
        The atom indices of the initial unit cell with no replicates.
    """
    num1 = np.shape(pts1)[0]
    pts_w_imgs = np.zeros((num1*(2*n1+1)*(2*n2+1), 3))
    inds_array = np.zeros((num1*(2*n1+1)*(2*n2+1), ))
    tinds1 = np.arange(0, num1)

    # The first set of atoms correspond to the main
    # cell.
    ct1 = 0
    ind_st = num1*ct1
    ind_stop = num1*(ct1+1)-1
    pts_w_imgs[ind_st:ind_stop+1, :] = pts1
    inds_array[ind_st:ind_stop+1] = tinds1
    ct1 = ct1 + 1

    # Array for translating the main cell
    n1_val = np.linspace(-n1, n1, 2*n1+1)
    n2_val = np.linspace(-n2, n2, 2*n2+1)
    mval = np.meshgrid(n1_val, n2_val)
    m1 = np.ndarray.flatten(mval[0])
    m2 = np.ndarray.flatten(mval[1])
    i1 = np.where((m1 == 0) & (m2 == 0))[0][0]
    m1 = np.delete(m1, i1)
    m2 = np.delete(m2, i1)
    p1_trans = np.tile(sim_1vec, (num1, 1))
    p2_trans = np.tile(sim_2vec, (num1, 1))

    # Creating the images
    for ct2 in range(np.size(m1)):
        mp1 = m1[ct2]
        mp2 = m2[ct2]
        pts_trans = pts1 + mp1*p1_trans + mp2*p2_trans

        ind_st = num1*ct1
        ind_stop = num1*(ct1+1)-1
        pts_w_imgs[ind_st:ind_stop+1, :] = pts_trans
        inds_array[ind_st:ind_stop+1] = tinds1
        ct1 = ct1 + 1

    return pts_w_imgs, inds_array.astype(int)
