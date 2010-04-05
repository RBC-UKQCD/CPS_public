sample_dir is a sample directory structure that can be used to run the contraction_production code and the associated scripts that loop over configurations.  You can copy this to the host directory and rename it to whatever you wish.

You need to specify the configuration numbers you wish to analyze.  The "add_configs.sh" script does this for you if give it a list of the configuration numbers in a text file (more details in the comments at the top of that script).

IMPORTANT: Make sure to put a copy of the binary in the binary_v2_8 subdirectory.

In addition, there will be some other things you will likely have to change in the scripts and the vml files to suit your particular setup.  Here is a summary:

pbs_script_v2_8/qbatch.pbs
-job name (line 16)
-PBS queue/machine name (line 19)
-topology of the machine (line 42)
-your e-mail address to get PBS messages (line 25, 45)

pbs_script_v2_8/user.qcsh
-the options that follow the "while" statement (line 16-34);  for lines 19-34 you can just replace "/host2/lightman/contraction_production_r08-11_run2" with the path of the directory that you change sample_dir into, and "/host/lightman/contraction_production_r08-11_run2" with the path that the machine would see for this directory
-the range of configuration numbers that exist and can be used (line 52)
-the command that copies this configuration number from its original location (line 53)
-the last four arguments of the qrun command are set to 0 so that propagators and gauge rotated lattices aren't saved (line 63); if you want to save propagators those four options correspond to saving light quark propagators, saving strange quark propagators, saving gauge rotated lattices, and overwriting the propagators saved from previous configurations

vml_files/do_arg.vml and vml_files/evo_arg.vml
-these are standard.  do_arg.vml needs to be changed, for example, to specify the lattice size you are using.

vml_files/threept_arg.vml
-the location to write out the contraction data (line 2-4); you can just change "/host/lightman/contraction_production_r08-11_run2/" as you did in the user.qcsh file
-many other options, flags, and inputs that are described in the comments of the threept_arg.h file found in the code directory (contraction_production_v2_8/), for example the light quark masses (line 24-46), the strange quark masses (line 47-69), and the times at which to place the kaon source (line 93-115)
