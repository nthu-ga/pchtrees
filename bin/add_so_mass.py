import numpy as np
import tables as tb
import os
import sys
import tempfile

import logging
log = logging.getLogger('add_so_mass')
log.setLevel(logging.INFO)

modpath = '/data/apcooper/sfw/colossus/'
if modpath not in sys.path:
    sys.path.append(modpath)
from colossus.cosmology import cosmology
from colossus.halo      import profile_nfw, mass_so, concentration

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


def convert_SO_masses(m_input, z_input, c_input=10.0, mdef_input='vir',mdef_output='200c'):
    """
    """
    # When computing concentrations from accretion histories later, we will
    # treat these fiducial masses as fixed, even though in principle they
    # depend on the choice of fiducial concentraiton. The inconsistency here
    # should be of order 1% in the total masses.

    # Use colossus routines to find NFW parameters for each halo.
    # These are 'fiducial' in the sense that they assume a concentration.
    # The 'virial' definition is the B&N fitting function.
    rvir = mass_so.M_to_R(m_input, z_input, mdef_input)
    rs   = rvir / c_input
    rhos = m_input / rs**3 / 4.0 / np.pi / profile_nfw.NFWProfile.mu(c_input)

    # Compute the density threshold for the output mass definition
    density_threshold = mass_so.densityThreshold(z_input, mdef_output)

    # Compute the radius enclosing the new density threshold:
    x = profile_nfw.NFWProfile.xDelta(rhos, density_threshold)
    R_out = x * rs

    # And the mass enclosed by that raidus:
    M_out = 4.0 / 3.0 * np.pi * density_threshold * R_out**3

    if np.isscalar(c_input):
        c_output = np.repeat(c_input,len(M_out))
    else:
        c_output = c_input

    return M_out, R_out, c_output, rvir

def main(filename):
    """
    """
    log.info(f'Reading trees from: {filename:s}')
    with tb.open_file(filename,'r') as f:
        redshift_list = f.root.OutputTimes.Redshift.read()[::-1]

        mvir_input = f.get_node('/TreeHalos','Group_M_Virial').read()
        snap_input = f.get_node('/TreeHalos','SnapNum').read()
        zred_input = redshift_list[snap_input]

        # Assume flat so don't read Lambda0
        Om0    = f.root.Parameters._v_attrs['cosmo_omega0'][0]
        Ob0    = f.root.Parameters._v_attrs['cosmo_omegab'][0]
        sigma8 = f.root.Parameters._v_attrs['pspec_sigma8'][0]
        ns     = f.root.Parameters._v_attrs['pspec_nspec'][0]

        ps_k  = f.get_node('/Powerspec','k').read()
        ps_pk = f.get_node('/Powerspec','Pk').read()

    # Set the colossus cosmology to that of the tree file (h = 1)
    colo_cosmo_data = {'flat': True, 'H0': 100.0, 'Om0': Om0, 'Ob0': Ob0, 'sigma8': sigma8, 'ns': ns}
    colo_cosmo = cosmology.setCosmology('pchtrees_input_cosmology', **colo_cosmo_data)

    # Write the tabulated PS to a temp file for the colossus routines
    temp_ps_file = tempfile.NamedTemporaryFile()
    np.savetxt(temp_ps_file.name, np.column_stack([np.log10(ps_k),np.log10(ps_pk)]))

    print()
    print(cosmology.current_cosmo)
    print()

    # Compute SO masses assuming fiducial concentration (virial mass definition).
    Cvir_input, c_err = concentration.modelIshiyama21(mvir_input,zred_input,'vir',
                                                   ps_args=dict(path=temp_ps_file.name,model='ps_func'))


    print('Concentration errors (1):', np.sum(~c_err), '/', len(c_err))

    M200c_out, R200c_out, _, Rvir_out = convert_SO_masses(mvir_input, zred_input,
                                                          c_input=Cvir_input,
                                                          mdef_input='vir',mdef_output='200c')

    C200c_out, c_err = concentration.modelIshiyama21(mvir_input,zred_input,'200c',
                                                   ps_args=dict(path=temp_ps_file.name,model='ps_func'))

    print('Concentration errors (2):', np.sum(~c_err), '/', len(c_err))

    append_hdf5_data(filename,'/TreeHalos','Group_M_Crit200',     M200c_out,
                    comment='M200c mass, h^-1 Msol',
                    overwrite=True)
    append_hdf5_data(filename,'/TreeHalos','Group_R_Crit200',     R200c_out,
                    comment='R200c radius, h^-1 kpc',
                    overwrite=True)
    append_hdf5_data(filename,'/TreeHalos','Group_R_Virial',      Rvir_out,
                    comment='Virial radius (Bryan & Norman), h^-1 kpc',
                    overwrite=True)
    append_hdf5_data(filename,'/TreeHalos','Group_C_Crit200_fid', C200c_out,
                    comment='Fiducial NFW concentration R_s/R_200c (Ishiyama21)',
                    overwrite = True)
    append_hdf5_data(filename,'/TreeHalos','Group_C_Virial_fid',  Cvir_input,
                    comment='Fiducial NFW concentration R_s/R_vir (Ishiyama21)',
                    overwrite = True)

    # Should stamp a ledger with the version of this script and the runtime, and a checksum?
    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')

    args = parser.parse_args()
    assert(os.path.exists(args.filename))
    main(args.filename)