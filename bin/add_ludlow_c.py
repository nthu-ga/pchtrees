import numpy as np
import tables as tb
import os
import sys
import tempfile
import scipy.interpolate as spi

import logging
log = logging.getLogger('add_so_mass')
log.setLevel(logging.INFO)

modpath = '/data/apcooper/sfw/colossus/'
if modpath not in sys.path:
    sys.path.append(modpath)
from colossus.cosmology import cosmology
from colossus.halo      import profile_nfw, mass_so, concentration

modpath = '/data/apcooper/sfw/coco_trees/py'
if modpath not in sys.path:
    sys.path.append(modpath)
import coco_trees.trees as trees

def append_hdf5_data(filename,dataset_path,dataset_name,data,
                     comment="",createparents=False,
                     complevel=6,
                     create_file=False,overwrite=False):
    """
    Appends a dataset to an hdf5 file.
    """
    mode = "a"
    if not os.path.exists(filename):
        if create_file:
            mode = "w"
        else:
            raise Exception("HDF5 file {} does not exist".format(filename))

    with tb.open_file(filename,mode) as t:
        if overwrite:
            try:
                log.debug('Removing node...')
                t.remove_node("{}/{}".format(dataset_path,dataset_name))
                log.debug('Removed node...')
            except tb.exceptions.NoSuchNodeError:
                pass

        filters = tb.Filters(complevel=complevel, complib='zlib')
        t.create_carray(dataset_path,
                       dataset_name,
                       obj=data,
                       title=comment,
                       filters=filters,
                       createparents=createparents)
    return

def compile_mcoll(root,f=0.02,truncate=True,convert=True,
                  mdef_from='vir', mres=None,fiducial_c=10,verbose=False):
    """
    Walk the tree about the given node (of mass M0, redshift z0) to obtain the
    'collapsed mass history', defined as a the sum of the mass of all progenitors
    of the node extant at a given redshift z > z0 with individual masses greater
    than f * M0.

    The definition of mass should be provided as mdef_from; if this is '200', the
    corresponding M200c masses should already have been added to the tree files,
    in which case they should have been read and included with each node automatically
    when building the tree.
    """
    nsnap = root.snapshot

    # Note that PCH Trees masses are h^-1
    mcoll = np.zeros(nsnap+1,dtype=float)

    # Specify the definition of the root node mass
    root_mass = root.m200c if mdef_from == '200c' else root.mass

    # Optionally impose a different mass resolution limit on the tree nodes
    # (for debugging)
    limiting_mass = f*root_mass if mres is None else np.maximum(f*root_mass,mres)

    if verbose:
        log.debug('Limiting mass (log10):', np.log10(limiting_mass))

    # Walk the full tree above this node
    walker = trees.TreeHostIterator(root)
    for n in walker:
        node_mass = n.m200c if mdef_from == '200c' else n.mass
        if node_mass > limiting_mass:
            mcoll[n.snapshot] += node_mass

    if truncate:
        # Only return the collapsed mass history where information
        # from the tree exists (i.e. where mcoll is nonzero, accounting
        # for 'noise' at early times by requiring a contiguous run)
        for i in range(root.snapshot+1):
            if mcoll[i] > 0:
                break
        return mcoll[i:]
    else:
        return mcoll

def iterate_for_c(mcoll, mcoll_z, z0=None, m0=None,
                  C_FACTOR = 400.0,
                  iter_max=20, eps=0.01,
                  verbose=False, c_trial=10.0):
    """
    """
    RHO_FACTOR = 3/(4*np.pi)
    converged  = False

    colo_cosmo = cosmology.getCurrent()

    minimum_mcoll = mcoll.min()

    # Interpolates a mass to a redshift
    interpolate_log_mcoll = spi.CubicSpline(np.log10(mcoll),mcoll_z)

    if z0 is None:
        z0 = mcoll_z[-1]
    if m0 is None:
        m0 = mcoll[-1]

    z_trial = list()
    error   = list()

    # The concentration-Mcoll relation from López-Cano et al. is based on a
    # 200m definition of halo mass. C = 493.

    # In Ludlow 2016, it is based on a 200c definition. C = 400.

    # PCH Trees masses are h^-1 and use the virial definition of halo
    # mass.

    # As long as we can calculate the mass within the scale radius, which has the
    # same meaning for all mass definitions, we should be OK. We can search for
    # the characteristic time, zc, for the native mass definition of the halo. At
    # the end, we can report the concentration according to some other mass definition.

    native_mdef = 'vir'
    nfw = profile_nfw.NFWProfile(M = m0, c = c_trial, z = z0, mdef = native_mdef)
    r_vir_native = nfw.RDelta(z0, native_mdef)

    # Starting guess
    r_s_trial = r_vir_native/c_trial

    rho_0 = nfw.par['rhos']

    if verbose:
        print(f"Starting: r_s = {r_s_trial:4.3f}")

    # Iterate to convergence
    for j in range(0,iter_max):
        # When does the collapsed mass equal the mass within r_s?
        # Note the enclosed mass does not depend on the definition of total mass.
        m_target = profile_nfw.NFWProfile.M(rho_0, r_s_trial, 1)

        if m_target < minimum_mcoll:
            zc = -1
            break

        zc = interpolate_log_mcoll(np.log10(m_target))

        if verbose:
            print(f"Step {j:d}: Target mass = {np.log10(m_target):4.3f}")
            print(f"Step {j:d}:          zc = {zc:4.3f}")

        # For debugging
        z_trial.append(zc)

        # Calculate < rho_(-2) >
        # This needs to be h^2 Msol/kpc^3 because the mass is h^-1 Msol etc.
        # (Colossus returns rho_crit in Msol/kpc^3)
        rho_scale = C_FACTOR*colo_cosmo.rho_c(zc)

        # Calculate updated scale factor
        # < rho_s > = M_s / [ (4/3)*pi*r_s**3 ]
        # => r_s = { [3 / (4pi)] M_s / < rho_s> }**(1/3)
        r_s = (RHO_FACTOR*(m_target/rho_scale))**(1/3)

        # Determine error in r_s
        new_error = np.abs(r_s  - r_s_trial)
        error.append(new_error)

        if verbose:
            print(f"Step {j:d}:       error = {new_error:4.3f}")

        if new_error < eps:
            # Converged, no need to update r_s_trial
            converged = True
            break

        # Have not converged, update trial r_s
        r_s_trial = r_s

        # Update corresponding trial concentration
        # (the reference rvir is fixed)
        c_trial = r_vir_native/r_s_trial

        rho_0 = m0 / r_s_trial**3 / 4.0 / np.pi / profile_nfw.NFWProfile.mu(c_trial)

        if verbose:
            print(f"Step {j:d}: r_s = {r_s_trial:4.3f}, c= {c_trial:4.3f}")

    retdict = dict()
    retdict['converged'] = converged
    retdict['n_iter']    = j
    retdict['n_mcoll']   = len(mcoll)
    retdict['zc']        = zc
    retdict['ztrial']    = np.array(z_trial)
    retdict['error']     = np.array(error)

    return c_trial, retdict

def main(filename):
    """
    """
    log.info(f'Reading trees from: {filename:s}')

    M = trees.PCHTreeFile(filename)

    with tb.open_file(filename,'r') as f:

        # FIXME: Confirm we have the required arrays

        # FIXME: these parameters could be read by the tree file reader.
        # Assume flat so don't read Lambda0
        Om0    = f.root.Parameters._v_attrs['cosmo_omega0'][0]
        Ob0    = f.root.Parameters._v_attrs['cosmo_omegab'][0]
        sigma8 = f.root.Parameters._v_attrs['pspec_sigma8'][0]
        ns     = f.root.Parameters._v_attrs['pspec_nspec'][0]

        # Set the colossus cosmology to that of the tree file (h = 1)
        colo_cosmo_data = {'flat': True, 'H0': 100.0, 'Om0': Om0, 'Ob0': Ob0, 'sigma8': sigma8, 'ns': ns}
        colo_cosmo = cosmology.setCosmology('pchtrees_input_cosmology', **colo_cosmo_data)

    MIN_MCOLL_LENGTH = 3
    LUDLOW_LITTLE_F  = 0.02
    LUDLOW_C_FACTOR  = 400

    # Arrays to hold concentrations for those halos we recompute here.
    # The default is -1.
    # FIXME some of these arrays could be float16
    Cvir_output  = np.zeros(M.nnodes,dtype=np.float32) - 1

    STATUS_NOT_ATTEMPTED = -1
    STATUS_OK            = 0
    STATUS_LIMITED_MCOLL = 1
    STATUS_LOW_MRES      = 2
    STATUS_NOT_CONVERGED = 3
    Cvir_status = np.repeat(STATUS_NOT_ATTEMPTED, M.nnodes)

    print(f'Have {(int(M.tree_table_ntrees)):d} trees to process')
    for itree in range(M.tree_table_ntrees):
        T = M.get_tree(itree)

        # FIXME: also add first order progenitors
        nodes_to_process = list()
        # Add concentrations to the main branch only
        MBW = T.host_main_branch_walker()
        for n in MBW:
            nodes_to_process.append(n)

        print(f'Tree {itree:d}: {len(nodes_to_process):d} nodes to process')

        for n in nodes_to_process:
            # Compute collapsed mass history in m200c
            # Requires the tree-reader to have already set m200c
            if np.any(M.mres[:n.snapshot]) > n.m200c*LUDLOW_LITTLE_F:
                # Mass resolution is too low at snapshots above this node.
                Cvir_status[n.offset_in_tree_file] = STATUS_LOW_MRES
                continue
            else:
                mcoll = compile_mcoll(n,f=LUDLOW_LITTLE_F,
                                      mdef_from='200c',convert=False)
            # FIXME: explain the following indexing
            mcoll_z = M.output_redshift[n.snapshot-len(mcoll)+1:n.snapshot+1] # Some thought needed

            # We give up trying to compute the concentration if:
            # - the collapsed mass history is too short
            # - the mass resolution is not high enough to resolve f*M0
            # - the iterations to compute the concentration do not converge

            # This trial value has already been computed by `add_so_mass` using some
            # fiducial mean relation.
            c_trial = n.c200c
            if len(mcoll) > MIN_MCOLL_LENGTH:
                c, diag = iterate_for_c(mcoll, mcoll_z, c_trial=c_trial,
                                        C_FACTOR = LUDLOW_C_FACTOR,
                                        verbose=False)
                if diag['converged']:
                    Cvir_output[n.offset_in_tree_file] = c
                    Cvir_status[n.offset_in_tree_file] = STATUS_OK
                else:
                    Cvir_status[n.offset_in_tree_file] = STATUS_NOT_CONVERGED

    print(f'Appending output to {filename:s}')
    append_hdf5_data(filename,'/TreeHalos','Group_C_Virial_L16_History',  Cvir_output,
                    comment='NFW concentration from L16 collapsed mass history',
                    overwrite = True)
    append_hdf5_data(filename,'/TreeHalos','Group_C_Virial_L16_History_Status',  Cvir_status,
                    comment='Status flag for computation of Group_C_Virial_L16_History',
                    overwrite = True)

    # FIXME diagnostics and logging
    print('Done!')

    # STORE PARAMETERS
    # Should stamp a ledger with the version of this script and the runtime, and a checksum?
    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')

    args = parser.parse_args()
    assert(os.path.exists(args.filename))
    main(args.filename)