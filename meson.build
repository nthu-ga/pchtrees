project('pchtrees', 'fortran',
  version : '0.1',
  default_options : ['warning_level=3'])

tomlf_dep = dependency(
  'toml-f',
  version: '>=0.2.0',
  fallback: ['toml-f', 'tomlf_dep'],
  default_options: ['default_library=static'],
)

#h5fortran_proj = subproject('h5fortran')
#h5fortran_dep  = h5fortran_proj.get_variable('h5fortran_dep')

hdf5 = dependency('hdf5', language: 'fortran', method: 'config-tool')

executable('pchtrees',
		'parameter_file.F90',
		'defined_types.F90',
		'hyperbolic.F90',
		'indexxx.F90',
		'interp.F90',
		'locate.F90',
		'make_tree.F90',
		'memory.F90',
		'memory_modules.F90',
		'modified_merger_tree.F90',
		'num_pars.F90',
		'ran3.F90',
		'sigmacdm_spline.F90',
		'spline.F90',
		'split_PCH.F90',
		'transfer_function.F90',
		'tree_routines.F90',
		'unresolved_mass.F90',
		'run_statistics.F90',
		'trees.F90',
		'kind_numbers.F90',
		'cosmological_parameters.F90',
		'time_parameters.F90',
		'power_spectrum_parameters.F90',
		'runtime_parameters.F90',
		'deltcrit.F90',
		'halo_mass_function.F90',
		'io.F90',
		'commandline.F90',
		# dependencies : [tomlf_dep, hdf5, h5fortran_dep],
		dependencies : [tomlf_dep, hdf5],
		install : true)
