#include<config.h>
CPS_START_NAMESPACE
Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page   1
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

* Error: Duplicate definition of global symbol.
	"_fopen" in section ".text" in module "sysfunc" in file
	"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
* Information: Ignoring redefinition in module:
	"stdio" in file "/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
* Error: Duplicate definition of global symbol.
	"_fprintf" in section ".text" in module "sysfunc" in file
	"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
* Information: Ignoring redefinition in module:
	"stdio" in file "/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
* Error: Duplicate definition of global symbol.
	"_printf" in section ".text" in module "sysfunc" in file
	"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
* Information: Ignoring redefinition in module:
	"stdio" in file "/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"

Allocation to Output Section "vector_overlay" in module *unnamed* in file
"/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/qcdsp.outtof"

  Offset   Length		Input Section

  809F95       2C	".text" in module "vector_util_asm_cram" in file
			"../../lib/util.olb"

  Total Allocation = 2C (hex)

Allocation to Output Section "wfm0_overlay" in module *unnamed* in file
"/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/qcdsp.outtof"

  Offset   Length		Input Section

  809800       E0	"wfm0" in module "wfm_trick_spproj" in file
			"../../lib/d_op_wilson_opt_nos.olb"
  8098E0        4	"wfm0" in module "wfm_spproj" in file
			"../../lib/d_op_wilson_opt_nos.olb"
  8098E4        C	"wfm0" in module "wfm_nga_reg" in file
			"../../lib/d_op_wilson_opt_nos.olb"
  8098F0      10B	"wfm0" in module "wfm_mat_trick" in file
			"../../lib/d_op_wilson_opt_nos.olb"
  8099FB       AF	"wfm0" in module "wfm_cmat_spproj" in file
			"../../lib/d_op_wilson_opt_nos.olb"
  809AAA       32	"wfm0" in module "wfm_buffers" in file
			"../../lib/d_op_wilson_opt_nos.olb"

  Total Allocation = 2DC (hex)

Allocation to Output Section "wfm1_overlay" in module *unnamed* in file

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page   2
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

"/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/qcdsp.outtof"

  Offset   Length		Input Section

  809C00        7	"wfm1" in module "wilson_dslash" in file
			"../../lib/d_op_wilson_opt_nos.olb"
  809C07        C	"wfm1" in module "wfm_trick_spproj" in file
			"../../lib/d_op_wilson_opt_nos.olb"
  809C13       4B	"wfm1" in module "wfm_trick_segment" in file
			"../../lib/d_op_wilson_opt_nos.olb"
  809C5E       69	"wfm1" in module "wfm_spproj" in file
			"../../lib/d_op_wilson_opt_nos.olb"
  809CC7        4	"wfm1" in module "wfm_dslash" in file
			"../../lib/d_op_wilson_opt_nos.olb"
  809CCB       DA	"wfm1" in module "wfm_buffers" in file
			"../../lib/d_op_wilson_opt_nos.olb"

  Total Allocation = 1A5 (hex)

Allocation to Output Section "cfm_overlay" in module *unnamed* in file
"/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/qcdsp.outtof"

  Offset   Length		Input Section

  809E00       C1	".text" in module "clover_mat_mlt_asm" in file
			"../../lib/d_op_clover.olb"

  Total Allocation = C1 (hex)

Allocation to Output Section "stag_ds_overlay" in module *unnamed* in file
"/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/qcdsp.outtof"

  Offset   Length		Input Section

  809800      282	".text" in module "dirac_serial" in file
			"../../lib/d_op_stag_opt.olb"

  Total Allocation = 282 (hex)
* Warning: No sections matched ALLOCATE Command specification:
	MODULE = "assign_cram0" KIND = CODE
* Warning: No sections matched ALLOCATE Command specification:
	MODULE = "assign_cram1"
* Warning: No sections matched ALLOCATE Command specification:
	MODULE = "cmhb_kern"
* Warning: No sections matched ALLOCATE Command specification:
	MODULE = "cmhb_kern_sup"
* Warning: No sections matched ALLOCATE Command specification:

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page   3
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

	MODULE = "grand"
* Warning: No sections matched ALLOCATE Command specification:
	MODULE = "mtxfast"

Allocation to Output Section ".bss" in module *unnamed* in file
"/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/qcdsp.outtof"

  Offset   Length		Input Section

    1000      42E	"GLOBAL" in module "main" in file
	    "/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/main.tof"
    142E       3C	"TL.Init" in module "main" in file
	    "/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/main.tof"
    146A        9	"GLOBAL" in module "alg_hmd_r" in file
			"../../lib/alg.olb"
    1473      1A0	"TL.Init" in module "alg_hmd_r" in file
			"../../lib/alg.olb"
    1613        9	"GLOBAL" in module "alg_hmd" in file "../../lib/alg.olb"
    161C       37	"TL.Init" in module "alg_hmd" in file
			"../../lib/alg.olb"
    1653        9	"GLOBAL" in module "alg_hmc_phi" in file
			"../../lib/alg.olb"
    165C      2FA	"TL.Init" in module "alg_hmc_phi" in file
			"../../lib/alg.olb"
    1956        9	"GLOBAL" in module "alg_base" in file
			"../../lib/alg.olb"
    195F       29	"TL.Init" in module "alg_base" in file
			"../../lib/alg.olb"
    1988        9	"GLOBAL" in module "verbose" in file
			"../../lib/util.olb"
    1991      219	"TL.Init" in module "verbose" in file
			"../../lib/util.olb"
    1BAA        1	".data" in module "su3_util_asm" in file
			"../../lib/util.olb"
    1BAB        1	".const" in module "su3_util_asm" in file
			"../../lib/util.olb"
    1BAC        1	".const" in module "sproj_tr_asm" in file
			"../../lib/util.olb"
    1BAD       16	"TL.Init" in module "smalloc" in file
			"../../lib/util.olb"
    1BC3        4	"OWN" in module "rcomplex" in file "../../lib/util.olb"
    1BC7       1E	".data" in module "random_asm" in file
			"../../lib/util.olb"
    1BE5       39	"GLOBAL" in module "random" in file "../../lib/util.olb"
    1C1E       16	"TL.Init" in module "pmalloc" in file
			"../../lib/util.olb"
    1C34       B1	"GLOBAL" in module "lattice_ws" in file

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page   4
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

			"../../lib/util.olb"
    1CE5       2C	"TL.Init" in module "lattice_ws" in file
			"../../lib/util.olb"
    1D11       B4	"GLOBAL" in module "lattice_wc" in file
			"../../lib/util.olb"
    1DC5       32	"TL.Init" in module "lattice_wc" in file
			"../../lib/util.olb"
    1DF7       68	"GLOBAL" in module "lattice_base" in file
			"../../lib/util.olb"
    1E5F       37	"OWN" in module "lattice_base" in file
			"../../lib/util.olb"
    1E96      524	"TL.Init" in module "lattice_base" in file
			"../../lib/util.olb"
    23BA       10	"GLOBAL" in module "gjp" in file "../../lib/util.olb"
    23CA      2E7	"TL.Init" in module "gjp" in file "../../lib/util.olb"
    26B1       6C	"GLOBAL" in module "g_wilson" in file
			"../../lib/util.olb"
    271D        7	"OWN" in module "g_wilson" in file "../../lib/util.olb"
    2724       70	"TL.Init" in module "g_wilson" in file
			"../../lib/util.olb"
    2794       63	"GLOBAL" in module "f_wilson_t" in file
			"../../lib/util.olb"
    27F7       69	"TL.Init" in module "f_wilson_t" in file
			"../../lib/util.olb"
    2860       60	"GLOBAL" in module "f_stag_t" in file
			"../../lib/util.olb"
    28C0       26	"TL.Init" in module "f_stag_t" in file
			"../../lib/util.olb"
    28E6       9C	"GLOBAL" in module "f_stag" in file "../../lib/util.olb"
    2982       37	"OWN" in module "f_stag" in file "../../lib/util.olb"
    29B9      14A	"TL.Init" in module "f_stag" in file
			"../../lib/util.olb"
    2B03       9F	"GLOBAL" in module "f_clover" in file
			"../../lib/util.olb"
    2BA2       10	"OWN" in module "f_clover" in file "../../lib/util.olb"
    2BB2      193	"TL.Init" in module "f_clover" in file
			"../../lib/util.olb"
    2D45       12	"GLOBAL" in module "error" in file "../../lib/util.olb"
    2D57      173	"TL.Init" in module "error" in file "../../lib/util.olb"
    2ECA        1	"GLOBAL" in module "convert_func" in file
			"../../lib/util.olb"
    2ECB        1	"OWN" in module "convert_func" in file
			"../../lib/util.olb"
    2ECC      119	"TL.Init" in module "convert_func" in file
			"../../lib/util.olb"
    2FE5        2	"GLOBAL" in module "convert" in file
			"../../lib/util.olb"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page   5
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

    2FE7        1	"OWN" in module "convert" in file "../../lib/util.olb"
    2FE8      1AD	"TL.Init" in module "convert" in file
			"../../lib/util.olb"
    3195        8	".data" in module "dirac_serial" in file
			"../../lib/d_op_stag_opt.olb"
    319D       A0	"GLOBAL" in module "stag_scu" in file
			"../../lib/d_op_stag_opt.olb"
    323D       23	"OWN" in module "stag_scu" in file
			"../../lib/d_op_stag_opt.olb"
    3260       CE	"GLOBAL" in module "dirac_init" in file
			"../../lib/d_op_stag_opt.olb"
    332E       2D	"GLOBAL" in module "d_op_stag" in file
			"../../lib/d_op_stag_opt.olb"
    335B       7A	"TL.Init" in module "d_op_stag" in file
			"../../lib/d_op_stag_opt.olb"
    33D5       E5	"TL.Init" in module "eigen_stag" in file
			"../../lib/d_op_stag_types.olb"
    34BA       2D	"GLOBAL" in module "d_op_stag_types" in file
			"../../lib/d_op_stag_types.olb"
    34E7       52	"TL.Init" in module "d_op_stag_types" in file
			"../../lib/d_op_stag_types.olb"
    3539        5	".data" in module "clover_mat_mlt_asm" in file
			"../../lib/d_op_clover.olb"
    353E       19	"OWN" in module "d_op_clover_supp" in file
			"../../lib/d_op_clover.olb"
    3557       37	"TL.Init" in module "d_op_clover_supp" in file
			"../../lib/d_op_clover.olb"
    358E       39	"GLOBAL" in module "d_op_clover" in file
			"../../lib/d_op_clover.olb"
    35C7       DA	"TL.Init" in module "d_op_clover" in file
			"../../lib/d_op_clover.olb"
    36A1       26	"OWN" in module "clover" in file
			"../../lib/d_op_clover.olb"
    36C7       20	"TL.Init" in module "clover" in file
			"../../lib/d_op_clover.olb"
    36E7       1B	"TL.Init" in module "lapack" in file
			"../../lib/lapack.olb"
    3702        8	"GLOBAL" in module "wilson_init" in file
			"../../lib/d_op_wilson_opt_nos.olb"
    370A       32	"TL.Init" in module "wilson_init" in file
			"../../lib/d_op_wilson_opt_nos.olb"
    373C       31	"TL.Init" in module "wilson_end" in file
			"../../lib/d_op_wilson_opt_nos.olb"
    376D       16	"OWN" in module "wfm_sublatt_pointers" in file
			"../../lib/d_op_wilson_opt_nos.olb"
    3783       1C	"TL.Init" in module "wfm_sublatt_pointers" in file
			"../../lib/d_op_wilson_opt_nos.olb"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page   6
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

    379F      15A	"TL.Init" in module "ritz" in file
			"../../lib/d_op_base.olb"
    38F9       5B	"TL.Init" in module "jacobi" in file
			"../../lib/d_op_base.olb"
    3954       DF	"TL.Init" in module "inv_cg" in file
			"../../lib/d_op_base.olb"
    3A33       2E	"GLOBAL" in module "dirac_op_base" in file
			"../../lib/d_op_base.olb"
    3A61        8	"OWN" in module "dirac_op_base" in file
			"../../lib/d_op_base.olb"
    3A69       81	"TL.Init" in module "dirac_op_base" in file
			"../../lib/d_op_base.olb"
    3AEA      15D	"TL.Init" in module "eigen_wilson" in file
			"../../lib/d_op_wilson_types.olb"
    3C47       39	"GLOBAL" in module "d_op_wilson_types" in file
			"../../lib/d_op_wilson_types.olb"
    3C80       E8	"TL.Init" in module "d_op_wilson_types" in file
			"../../lib/d_op_wilson_types.olb"
    3D68        9	"OWN" in module "glb_sum_five" in file
			"../../lib/glb_cpp.olb"
    3D71        9	"OWN" in module "glb_sum" in file
			"../../lib/glb_cpp.olb"
    3D7A      1C7	"TL.Init" in module "p2v" in file "../../lib/mem.olb"
    3F41        5	"OWN" in module "cbuf" in file "../../lib/cbuf.olb"
    3F46       50	".bss" in module "tclowio" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
    3F96        1	"iodata" in module "tclowio" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
    3F97        2	"OWN" in module "mypf" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
    3F99        7	"TL.Init" in module "mypf" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
    3FA0       4E	"GLOBAL" in module "cvtf64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
    3FEE        C	"OWN" in module "cvtdbl" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
    3FFA        1	".bss" in module "qcdsp_tcroot" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    3FFB        1	"OWN" in module "new" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    3FFC        2	"OWN" in module "memory" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    3FFE       3A	"GLOBAL" in module "math64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    4038       21	"OWN" in module "inifin" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    4059        2	".bss" in module "roundu32f64" in file

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page   7
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    405B        2	".bss" in module "mpyf40" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    405D        1	".bss" in module "flooru32f64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    405E        1	".bss" in module "floori32f64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    405F       3F	".bss" in module "f64trig" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    409E       15	".bss" in module "f64log" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    40B3       44	".bss" in module "f64itrig" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    40F7       35	".bss" in module "f64ihyp" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    412C       2D	".bss" in module "f64hyp" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    4159       68	".bss" in module "f64exp" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    41C1        F	".const" in module "f64exp" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    41D0       48	".bss" in module "f64const" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    4218        1	".bss" in module "ceilu32f64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    4219        2	".bss" in module "addf64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    421B        1	".bss" in module "addf40" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    421C        1	"errno" in module "trig" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
    421D    10000	".stack" in file "../../mem/include/link_p2v.lcf"
   1421D    30000	".sysmem" in file "../../mem/include/link_p2v.lcf"

  Total Allocation = 4321D (hex)

Allocation to Output Section ".data" in module *unnamed* in file
"/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/qcdsp.outtof"

  Offset   Length		Input Section

   4421D        0	"DEFALT" in module "alg_hmd_r" in file
			"../../lib/alg.olb"
   4421D        0	"DEFALT" in module "unitarize" in file
			"../../lib/util.olb"
   4421D        0	"DEFALT" in module "ritz" in file
			"../../lib/d_op_base.olb"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page   8
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

   4421D        0	"DEFALT" in module "jacobi" in file
			"../../lib/d_op_base.olb"
   4421D        0	"DEFALT" in module "inv_cg" in file
			"../../lib/d_op_base.olb"
   4421D        0	"DEFALT" in module "eigen_wilson" in file
			"../../lib/d_op_wilson_types.olb"
   4421D        1	"DEFALT" in module "tclowio" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   4421E        0	"DEFALT" in module "cvtdbl" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   4421E        D	"DEFALT" in module "qcdsp_tcroot" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   4422B        4	"DEFALT" in module "trig" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   4422F        3	"DEFALT" in module "trig" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   44232        4	"DEFALT" in module "fp_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   44236        4	".cinit" in module "roundu32f64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   4423A        4	".cinit" in module "mpyf40" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   4423E        3	".cinit" in module "flooru32f64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   44241        3	".cinit" in module "floori32f64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   44244       5B	".cinit" in module "f64trig" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   4429F       21	".cinit" in module "f64log" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   442C0       58	".cinit" in module "f64itrig" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   44318       3B	".cinit" in module "f64ihyp" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   44353       49	".cinit" in module "f64hyp" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   4439C       8C	".cinit" in module "f64exp" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   44428       90	".cinit" in module "f64const" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   444B8        3	".cinit" in module "ceilu32f64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   444BB        4	".cinit" in module "addf64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   444BF        3	".cinit" in module "addf40" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   444C2     2A49	".cinit" in file "../../mem/include/link_p2v.lcf"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page   9
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

   46F0B        F	"TL.Const" in module "main" in file
	    "/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/main.tof"
   46F1A        1	"TL.Const" in module "alg_hmd_r" in file
			"../../lib/alg.olb"
   46F1B        1	"TL.Const" in module "alg_hmc_phi" in file
			"../../lib/alg.olb"
   46F1C        1	"TL.Const" in module "vector" in file
			"../../lib/util.olb"
   46F1D        1	"TL.Const" in module "sigmaproj_tr" in file
			"../../lib/util.olb"
   46F1E        2	"TL.Const" in module "lattice_base" in file
			"../../lib/util.olb"
   46F20        1	"TL.Const" in module "f_stag" in file
			"../../lib/util.olb"
   46F21        2	"TL.Const" in module "eigen_stag" in file
			"../../lib/d_op_stag_types.olb"
   46F23        4	"TL.Const" in module "d_op_clover_supp" in file
			"../../lib/d_op_clover.olb"
   46F27        1	"GLOBALX" in module "clover" in file
			"../../lib/d_op_clover.olb"
   46F28        1	"OWNX" in module "lapack" in file "../../lib/lapack.olb"
   46F29        2	"TL.Const" in module "ritz" in file
			"../../lib/d_op_base.olb"
   46F2B        2	"TL.Const" in module "eigen_wilson" in file
			"../../lib/d_op_wilson_types.olb"
   46F2D        2	"TL.Const" in module "cvtf64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   46F2F        4	"TL.Const" in module "cvtdbl" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   46F33       14	"TL.Const" in module "math64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   46F47       3C	"POWERT" in module "trig" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   46F83       29	"ALOGT" in module "trig" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   46FAC        B	"Init_Table" in file "../../mem/include/link_p2v.lcf"
   46FB7        6	"Finalize_Table" in file
			"../../mem/include/link_p2v.lcf"

  Total Allocation = 2DA0 (hex)

Allocation to Output Section ".text" in module *unnamed* in file
"/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/qcdsp.outtof"

  Offset   Length		Input Section

   46FBD      3A1	"_main" in module "main" in file

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  10
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

	    "/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/main.tof"
   4735E        9	"STATIC_1" in module "main" in file
	    "/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/main.tof"
   47367        9	"STATIC_2" in module "main" in file
	    "/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/main.tof"
   47370        9	"STATIC_3" in module "main" in file
	    "/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/main.tof"
   47379       12	"Ini$main_C" in module "main" in file
	    "/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/main.tof"
   4738B       18	"Fin$main_C" in module "main" in file
	    "/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/main.tof"
   473A3      1AF	"___ct__7AlgHmdRFR7LatticeP9CommonArgP6HmdArg" in
			module "alg_hmd_r" in file "../../lib/alg.olb"
   47552       B9	"___dt__7AlgHmdRFv" in module "alg_hmd_r" in file
			"../../lib/alg.olb"
   4760B      5B4	"_run__7AlgHmdRFv" in module "alg_hmd_r" in file
			"../../lib/alg.olb"
   47BBF       71	"___ct__6AlgHmdFR7LatticeP9CommonArgP6HmdArg" in
			module "alg_hmd" in file "../../lib/alg.olb"
   47C30       3E	"___dt__6AlgHmdFv" in module "alg_hmd" in file
			"../../lib/alg.olb"
   47C6E      3AD	"___ct__9AlgHmcPhiFR7LatticeP9CommonArgP6HmdArg" in
			module "alg_hmc_phi" in file "../../lib/alg.olb"
   4801B      196	"___dt__9AlgHmcPhiFv" in module "alg_hmc_phi" in file
			"../../lib/alg.olb"
   481B1      968	"_run__9AlgHmcPhiFv" in module "alg_hmc_phi" in file
			"../../lib/alg.olb"
   48B19       40	"___ct__3AlgFR7LatticeP9CommonArg" in module "alg_base"
			in file "../../lib/alg.olb"
   48B59       29	"___dt__3AlgFv" in module "alg_base" in file
			"../../lib/alg.olb"
   48B82        5	"_AlgLattice__3AlgFv" in module "alg_base" in file
			"../../lib/alg.olb"
   48B87       2D	"___ct__7VerboseFv" in module "verbose" in file
			"../../lib/util.olb"
   48BB4       11	"___dt__7VerboseFv" in module "verbose" in file
			"../../lib/util.olb"
   48BC5       7B	"_Level__7VerboseFi" in module "verbose" in file
			"../../lib/util.olb"
   48C40       3F	"_Active__7VerboseF16VerboseLevelType" in module
			"verbose" in file "../../lib/util.olb"
   48C7F       2F	"_Func__7VerboseFPcT1" in module "verbose" in file
			"../../lib/util.olb"
   48CAE       2F	"_FuncEnd__7VerboseFPcT1" in module "verbose" in file
			"../../lib/util.olb"
   48CDD       1D	"_Pmalloc__7VerboseFPcN21Pvi" in module "verbose" in
			file "../../lib/util.olb"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  11
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

   48CFA       1D	"_Smalloc__7VerboseFPcN21Pvi" in module "verbose" in
			file "../../lib/util.olb"
   48D17       1B	"_Sfree__7VerboseFPcN21Pv" in module "verbose" in file
			"../../lib/util.olb"
   48D32       3E	"_Flow__7VerboseFPcT1PCce" in module "verbose" in file
			"../../lib/util.olb"
   48D70       29	"_Input__7VerboseFPcT1PCce" in module "verbose" in
			file "../../lib/util.olb"
   48D99       29	"_Result__7VerboseFPcT1PCce" in module "verbose" in
			file "../../lib/util.olb"
   48DC2       5F	"_Warn__7VerboseFPcT1PCce" in module "verbose" in file
			"../../lib/util.olb"
   48E21       29	"_Debug__7VerboseFPcT1PCce" in module "verbose" in
			file "../../lib/util.olb"
   48E4A       1F	"_Debug__7VerboseFPCce" in module "verbose" in file
			"../../lib/util.olb"
   48E69       12	"_LedOn__7VerboseFPcT1" in module "verbose" in file
			"../../lib/util.olb"
   48E7B       12	"_LedOff__7VerboseFPcT1" in module "verbose" in file
			"../../lib/util.olb"
   48E8D       37	"_LedFlash__7VerboseFPcT1i" in module "verbose" in
			file "../../lib/util.olb"
   48EC4       31	"_Clock__7VerboseFPcT1PCce" in module "verbose" in
			file "../../lib/util.olb"
   48EF5      146	".text" in module "vector_util_asm" in file
			"../../lib/util.olb"
   4903B        9	"STATIC_1" in module "vector" in file
			"../../lib/util.olb"
   49044        9	"STATIC_4" in module "vector" in file
			"../../lib/util.olb"
   4904D        9	"STATIC_5" in module "vector" in file
			"../../lib/util.olb"
   49056       20	"___ct__6MatrixFv" in module "vector" in file
			"../../lib/util.olb"
   49076       2A	"___ct__6MatrixFRC6Matrix" in module "vector" in file
			"../../lib/util.olb"
   490A0       14	"___as__6MatrixFf" in module "vector" in file
			"../../lib/util.olb"
   490B4       47	"_Trans__6MatrixFPCf" in module "vector" in file
			"../../lib/util.olb"
   490FB       44	"_RandLink__6MatrixFv" in module "vector" in file
			"../../lib/util.olb"
   4913F        E	"_UnitMatrix__6MatrixFv" in module "vector" in file
			"../../lib/util.olb"
   4914D        A	"_ZeroMatrix__6MatrixFv" in module "vector" in file
			"../../lib/util.olb"
   49157       20	"___ct__6VectorFv" in module "vector" in file

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  12
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

			"../../lib/util.olb"
   49177       19	"_NormSqGlbSum__6VectorFi" in module "vector" in file
			"../../lib/util.olb"
   49190       1A	"_ReDotProductGlbSum__6VectorFPC6Vectori" in module
			"vector" in file "../../lib/util.olb"
   491AA       25	"_CompDotProductGlbSum__6VectorFPC6Vectori" in module
			"vector" in file "../../lib/util.olb"
   491CF       4D	"_RandGaussVector__6VectorF6rfloati" in module "vector"
			in file "../../lib/util.olb"
   4921C       4B	"_normalize__FP6rfloat" in module "unitarize" in file
			"../../lib/util.olb"
   49267      119	"_orthogonalize__FP6rfloatPC6rfloat" in module
			"unitarize" in file "../../lib/util.olb"
   49380      2DB	"_crossProductConj__FP6rfloatPC6rfloatT2" in module
			"unitarize" in file "../../lib/util.olb"
   4965B       1E	"_Unitarize__6MatrixFv" in module "unitarize" in file
			"../../lib/util.olb"
   49679       2B	"_TrLessAntiHermMatrix__6MatrixFRC6Matrix" in module
			"su3_util_asm" in file "../../lib/util.olb"
   496A4       23	"_NegHalfTrSquare__6MatrixCFv" in module "su3_util_asm"
			in file "../../lib/util.olb"
   496C7        8	"_ReTr__6MatrixCFv" in module "su3_util_asm" in file
			"../../lib/util.olb"
   496CF       47	"_Cross2__6MatrixFRC6VectorT1" in module "su3_util_asm"
			in file "../../lib/util.olb"
   49716       36	"_AntiHermMatrix__6MatrixFPCf" in module "su3_util_asm"
			in file "../../lib/util.olb"
   4974C       23	"_Dagger__6MatrixFPCf" in module "su3_util_asm" in
			file "../../lib/util.olb"
   4976F       8F	"RTN8" in module "sproj_tr_asm" in file
			"../../lib/util.olb"
   497FE       8F	"RTN5" in module "sproj_tr_asm" in file
			"../../lib/util.olb"
   4988D       8F	"RTN2" in module "sproj_tr_asm" in file
			"../../lib/util.olb"
   4991C       8F	"RTN6" in module "sproj_tr_asm" in file
			"../../lib/util.olb"
   499AB       8F	"RTN3" in module "sproj_tr_asm" in file
			"../../lib/util.olb"
   49A3A       8F	"RTN7" in module "sproj_tr_asm" in file
			"../../lib/util.olb"
   49AC9       8F	"RTN4" in module "sproj_tr_asm" in file
			"../../lib/util.olb"
   49B58       8F	"RTN1" in module "sproj_tr_asm" in file
			"../../lib/util.olb"
   49BE7       23	"_smalloc__Fi" in module "smalloc" in file
			"../../lib/util.olb"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  13
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

   49C0A       1E	"_sfree__FPv" in module "smalloc" in file
			"../../lib/util.olb"
   49C28       69	"_SigmaprojTrXY__FPfN21iN24" in module "sigmaproj_tr"
			in file "../../lib/util.olb"
   49C91       69	"_SigmaprojTrYX__FPfN21iN24" in module "sigmaproj_tr"
			in file "../../lib/util.olb"
   49CFA       71	"_SigmaprojTrXZ__FPfN21iN24" in module "sigmaproj_tr"
			in file "../../lib/util.olb"
   49D6B       71	"_SigmaprojTrZX__FPfN21iN24" in module "sigmaproj_tr"
			in file "../../lib/util.olb"
   49DDC       71	"_SigmaprojTrXT__FPfN21iN24" in module "sigmaproj_tr"
			in file "../../lib/util.olb"
   49E4D       72	"_SigmaprojTrTX__FPfN21iN24" in module "sigmaproj_tr"
			in file "../../lib/util.olb"
   49EBF       72	"_SigmaprojTrYZ__FPfN21iN24" in module "sigmaproj_tr"
			in file "../../lib/util.olb"
   49F31       71	"_SigmaprojTrZY__FPfN21iN24" in module "sigmaproj_tr"
			in file "../../lib/util.olb"
   49FA2       71	"_SigmaprojTrYT__FPfN21iN24" in module "sigmaproj_tr"
			in file "../../lib/util.olb"
   4A013       71	"_SigmaprojTrTY__FPfN21iN24" in module "sigmaproj_tr"
			in file "../../lib/util.olb"
   4A084       69	"_SigmaprojTrZT__FPfN21iN24" in module "sigmaproj_tr"
			in file "../../lib/util.olb"
   4A0ED       69	"_SigmaprojTrTZ__FPfN21iN24" in module "sigmaproj_tr"
			in file "../../lib/util.olb"
   4A156        D	"___adv__6rfloatFf" in module "rfloat_rnd_op" in file
			"../../lib/util.olb"
   4A163        9	"___amu__6rfloatFf" in module "rfloat_rnd_op" in file
			"../../lib/util.olb"
   4A16C        9	"___ami__6rfloatFf" in module "rfloat_rnd_op" in file
			"../../lib/util.olb"
   4A175        9	"___apl__6rfloatFf" in module "rfloat_rnd_op" in file
			"../../lib/util.olb"
   4A17E        B	"___mi__FRC6rfloat" in module "rfloat_rnd_op" in file
			"../../lib/util.olb"
   4A189       1F	"___pl__FRC6rfloatT1" in module "rfloat" in file
			"../../lib/util.olb"
   4A1A8       1E	"___pl__FdRC6rfloat" in module "rfloat" in file
			"../../lib/util.olb"
   4A1C6       1E	"___pl__FRC6rfloatd" in module "rfloat" in file
			"../../lib/util.olb"
   4A1E4       1F	"___mi__FRC6rfloatT1" in module "rfloat" in file
			"../../lib/util.olb"
   4A203       1E	"___mi__FdRC6rfloat" in module "rfloat" in file
			"../../lib/util.olb"
   4A221       1E	"___mi__FRC6rfloatd" in module "rfloat" in file

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  14
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

			"../../lib/util.olb"
   4A23F       1F	"___ml__FRC6rfloatT1" in module "rfloat" in file
			"../../lib/util.olb"
   4A25E       1E	"___ml__FdRC6rfloat" in module "rfloat" in file
			"../../lib/util.olb"
   4A27C       1E	"___ml__FRC6rfloatd" in module "rfloat" in file
			"../../lib/util.olb"
   4A29A       1F	"___dv__FRC6rfloatT1" in module "rfloat" in file
			"../../lib/util.olb"
   4A2B9       1E	"___dv__FdRC6rfloat" in module "rfloat" in file
			"../../lib/util.olb"
   4A2D7       1E	"___dv__FRC6rfloatd" in module "rfloat" in file
			"../../lib/util.olb"
   4A2F5       16	"___ct__6rfloatFf" in module "rfloat" in file
			"../../lib/util.olb"
   4A30B       17	"___ct__6rfloatFRC6rfloat" in module "rfloat" in file
			"../../lib/util.olb"
   4A322        C	"___dt__6rfloatFv" in module "rfloat" in file
			"../../lib/util.olb"
   4A32E        A	"_norm__8RcomplexCFv" in module "rcomplex_asm" in file
			"../../lib/util.olb"
   4A338        B	"___amu__8RcomplexFf" in module "rcomplex_asm" in file
			"../../lib/util.olb"
   4A343        E	"___mi__FRC8Rcomplex" in module "rcomplex_asm" in file
			"../../lib/util.olb"
   4A351        D	"_conj__FRC8Rcomplex" in module "rcomplex_asm" in file
			"../../lib/util.olb"
   4A35E       1C	"___amu__8RcomplexFRC8Rcomplex" in module "rcomplex_asm"
			 in file "../../lib/util.olb"
   4A37A        C	"___ami__8RcomplexFRC8Rcomplex" in module "rcomplex_asm"
			 in file "../../lib/util.olb"
   4A386       1F	"Ini$rcomplex_C" in module "rcomplex" in file
			"../../lib/util.olb"
   4A3A5       18	"___ct__8RcomplexFfT1" in module "rcomplex" in file
			"../../lib/util.olb"
   4A3BD       19	"___ct__8RcomplexFRC8Rcomplex" in module "rcomplex" in
			file "../../lib/util.olb"
   4A3D6        C	"___dt__8RcomplexFv" in module "rcomplex" in file
			"../../lib/util.olb"
   4A3E2        A	"___as__8RcomplexFRC8Rcomplex" in module "rcomplex" in
			file "../../lib/util.olb"
   4A3EC       42	"_Rand__23GaussianRandomGeneratorFv" in module
			"random_asm" in file "../../lib/util.olb"
   4A42E        F	"_Rand__22UniformRandomGeneratorFv" in module
			"random_asm" in file "../../lib/util.olb"
   4A43D       3A	"_Rand__15RandomGeneratorFv" in module "random_asm" in
			file "../../lib/util.olb"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  15
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

   4A477       64	"_Reset__15RandomGeneratorSFi" in module "random" in
			file "../../lib/util.olb"
   4A4DB       23	"_pmalloc__Fi" in module "pmalloc" in file
			"../../lib/util.olb"
   4A4FE        2	"Freflex__7LatticeFP6VectorT1" in module "lattice_ws"
			in file "../../lib/util.olb"
   4A500       1E	"FsiteOffsetChkb__5FstagCFPCi" in module "lattice_ws"
			in file "../../lib/util.olb"
   4A51E       5D	"___ct__12GwilsonFstagFv" in module "lattice_ws" in
			file "../../lib/util.olb"
   4A57B       52	"___dt__12GwilsonFstagFv" in module "lattice_ws" in
			file "../../lib/util.olb"
   4A5CD        2	"Freflex__7LatticeFP6VectorT1" in module "lattice_wc"
			in file "../../lib/util.olb"
   4A5CF       5D	"___ct__14GwilsonFcloverFv" in module "lattice_wc" in
			file "../../lib/util.olb"
   4A62C       52	"___dt__14GwilsonFcloverFv" in module "lattice_wc" in
			file "../../lib/util.olb"
   4A67E       43	"Ini$lattice_base_C" in module "lattice_base" in file
			"../../lib/util.olb"
   4A6C1       38	"Fin$lattice_base_C" in module "lattice_base" in file
			"../../lib/util.olb"
   4A6F9       37	"_GetLinkOld__7LatticeCFP6MatrixPCiiT3" in module
			"lattice_base" in file "../../lib/util.olb"
   4A730       3C	"_MltFloatImpl__7LatticeF6rfloati" in module
			"lattice_base" in file "../../lib/util.olb"
   4A76C      1E1	"___ct__7LatticeFv" in module "lattice_base" in file
			"../../lib/util.olb"
   4A94D       75	"___dt__7LatticeFv" in module "lattice_base" in file
			"../../lib/util.olb"
   4A9C2        4	"_GaugeField__7LatticeCFv" in module "lattice_base" in
			file "../../lib/util.olb"
   4A9C6       27	"_GaugeField__7LatticeFP6Matrix" in module
			"lattice_base" in file "../../lib/util.olb"
   4A9ED       27	"_CopyGaugeField__7LatticeFP6Matrix" in module
			"lattice_base" in file "../../lib/util.olb"
   4AA14        4	"_StrOrd__7LatticeFv" in module "lattice_base" in file
			"../../lib/util.olb"
   4AA18        4	"_Colors__7LatticeFv" in module "lattice_base" in file
			"../../lib/util.olb"
   4AA1C       14	"_GsiteSize__7LatticeFv" in module "lattice_base" in
			file "../../lib/util.olb"
   4AA30      18D	"_Staple__7LatticeFR6MatrixPii" in module "lattice_base"
			 in file "../../lib/util.olb"
   4ABBD       A6	"_ReTrPlaq__7LatticeCFPiiT2" in module "lattice_base"
			in file "../../lib/util.olb"
   4AC63       92	"_SumReTrPlaqNode__7LatticeCFv" in module "lattice_base"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  16
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

			 in file "../../lib/util.olb"
   4ACF5       2D	"_SumReTrPlaq__7LatticeCFv" in module "lattice_base" in
			file "../../lib/util.olb"
   4AD22       5D	"_Reunitarize__7LatticeFv" in module "lattice_base" in
			file "../../lib/util.olb"
   4AD7F      151	"_Reunitarize__7LatticeFR6rfloatT1" in module
			"lattice_base" in file "../../lib/util.olb"
   4AED0       DE	"_MetropolisAccept__7LatticeF6rfloat" in module
			"lattice_base" in file "../../lib/util.olb"
   4AFAE       9E	"_EvolveGfield__7LatticeFP6Matrix6rfloat" in module
			"lattice_base" in file "../../lib/util.olb"
   4B04C       50	"_MomHamiltonNode__7LatticeFP6Matrix" in module
			"lattice_base" in file "../../lib/util.olb"
   4B09C       81	"_RandGaussAntiHermMatrix__7LatticeFP6Matrix6rfloat" in
			module "lattice_base" in file "../../lib/util.olb"
   4B11D       33	"_SetGfieldOrd__7LatticeFv" in module "lattice_base" in
			file "../../lib/util.olb"
   4B150       33	"_SetGfieldDisOrd__7LatticeFv" in module "lattice_base"
			in file "../../lib/util.olb"
   4B183        6	"_GupdCnt__7LatticeFv" in module "lattice_base" in
			file "../../lib/util.olb"
   4B189       10	"_GupdCntInc__7LatticeFi" in module "lattice_base" in
			file "../../lib/util.olb"
   4B199        F	"_MdTime__7LatticeFv" in module "lattice_base" in file
			"../../lib/util.olb"
   4B1A8       15	"_MdTime__7LatticeF6rfloat" in module "lattice_base" in
			file "../../lib/util.olb"
   4B1BD       19	"_MdTimeInc__7LatticeF6rfloat" in module "lattice_base"
			in file "../../lib/util.olb"
   4B1D6        5	"_Aux0Ptr__7LatticeFv" in module "lattice_base" in
			file "../../lib/util.olb"
   4B1DB        5	"_Aux1Ptr__7LatticeFv" in module "lattice_base" in
			file "../../lib/util.olb"
   4B1E0      11D	"_GsoCheck__7LatticeFv" in module "lattice_base" in
			file "../../lib/util.olb"
   4B2FD       9A	"_SoCheck__7LatticeF6rfloat" in module "lattice_base"
			in file "../../lib/util.olb"
   4B397       16	"_Gamma5__7LatticeFP6VectorT1i" in module "lattice_base"
			 in file "../../lib/util.olb"
   4B3AD       16	"_Ffour2five__7LatticeFP6VectorT1iT3" in module
			"lattice_base" in file "../../lib/util.olb"
   4B3C3       16	"_Ffive2four__7LatticeFP6VectorT1iT3" in module
			"lattice_base" in file "../../lib/util.olb"
   4B3D9       1B	"_FrandGaussVector__7LatticeFP6Vector6rfloati" in
			module "lattice_base" in file "../../lib/util.olb"
   4B3F4        2	"Freflex__7LatticeFP6VectorT1" in module "lattice_base"
			in file "../../lib/util.olb"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  17
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

   4B3F6       5B	"___ct__18GlobalJobParameterFv" in module "gjp" in
			file "../../lib/util.olb"
   4B451       5B	"___dt__18GlobalJobParameterFv" in module "gjp" in
			file "../../lib/util.olb"
   4B4AC      350	"_Initialize__18GlobalJobParameterFRC5DoArg" in module
			"gjp" in file "../../lib/util.olb"
   4B7FC       19	"RTN1" in module "gamma_5_asm" in file
			"../../lib/util.olb"
   4B815        D	"Ini$g_wilson_C" in module "g_wilson" in file
			"../../lib/util.olb"
   4B822        2	"Freflex__7LatticeFP6VectorT1" in module "g_wilson" in
			file "../../lib/util.olb"
   4B824       3C	"___ct__7GwilsonFv" in module "g_wilson" in file
			"../../lib/util.olb"
   4B860       3A	"___dt__7GwilsonFv" in module "g_wilson" in file
			"../../lib/util.olb"
   4B89A        3	"_Gclass__7GwilsonFv" in module "g_wilson" in file
			"../../lib/util.olb"
   4B89D       46	"_GactionGradient__7GwilsonFR6MatrixPii" in module
			"g_wilson" in file "../../lib/util.olb"
   4B8E3       A3	"_GforceSite__7GwilsonFR6MatrixPii" in module "g_wilson"
			 in file "../../lib/util.olb"
   4B986       C5	"_EvolveMomGforce__7GwilsonFP6Matrix6rfloat" in module
			"g_wilson" in file "../../lib/util.olb"
   4BA4B       4E	"_GhamiltonNode__7GwilsonFv" in module "g_wilson" in
			file "../../lib/util.olb"
   4BA99        2	"Freflex__7LatticeFP6VectorT1" in module "f_wilson_t"
			in file "../../lib/util.olb"
   4BA9B       A0	"___ct__12FwilsonTypesFv" in module "f_wilson_t" in
			file "../../lib/util.olb"
   4BB3B       3A	"___dt__12FwilsonTypesFv" in module "f_wilson_t" in
			file "../../lib/util.olb"
   4BB75       4A	"_Gamma5__12FwilsonTypesFP6VectorT1i" in module
			"f_wilson_t" in file "../../lib/util.olb"
   4BBBF        2	"Freflex__7LatticeFP6VectorT1" in module "f_stag_t" in
			file "../../lib/util.olb"
   4BBC1       3C	"___ct__10FstagTypesFv" in module "f_stag_t" in file
			"../../lib/util.olb"
   4BBFD       3A	"___dt__10FstagTypesFv" in module "f_stag_t" in file
			"../../lib/util.olb"
   4BC37       28	"Ini$f_stag_C" in module "f_stag" in file
			"../../lib/util.olb"
   4BC5F       44	"Fin$f_stag_C" in module "f_stag" in file
			"../../lib/util.olb"
   4BCA3        2	"Freflex__7LatticeFP6VectorT1" in module "f_stag" in
			file "../../lib/util.olb"
   4BCA5       D8	"_getUDagX__5FstagCFR6VectorPC6VectorPii" in module

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  18
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

			"f_stag" in file "../../lib/util.olb"
   4BD7D      171	"___ct__5FstagFv" in module "f_stag" in file
			"../../lib/util.olb"
   4BEEE       49	"___dt__5FstagFv" in module "f_stag" in file
			"../../lib/util.olb"
   4BF37        3	"_Fclass__5FstagFv" in module "f_stag" in file
			"../../lib/util.olb"
   4BF3A       1E	"FsiteOffsetChkb__5FstagCFPCi" in module "f_stag" in
			file "../../lib/util.olb"
   4BF58       15	"_FsiteOffset__5FstagCFPCi" in module "f_stag" in file
			"../../lib/util.olb"
   4BF6D        3	"_ExactFlavors__5FstagFv" in module "f_stag" in file
			"../../lib/util.olb"
   4BF70        3	"_SpinComponents__5FstagFv" in module "f_stag" in file
			"../../lib/util.olb"
   4BF73       1D	"_FsiteSize__5FstagFv" in module "f_stag" in file
			"../../lib/util.olb"
   4BF90        3	"_FchkbEvl__5FstagFv" in module "f_stag" in file
			"../../lib/util.olb"
   4BF93       45
		     "_FmatEvlInv__5FstagFP6VectorT1P5CgArgP6rfloat10CnvFrmType"
			 in module "f_stag" in file "../../lib/util.olb"
   4BFD8       1C	"_FmatEvlInv__5FstagFP6VectorT1P5CgArg10CnvFrmType" in
			module "f_stag" in file "../../lib/util.olb"
   4BFF4       3B
	  "_FmatInv__5FstagFP6VectorT1P5CgArgP6rfloat10CnvFrmType12PreserveType"
			 in module "f_stag" in file "../../lib/util.olb"
   4C02F       1E
		  "_FmatInv__5FstagFP6VectorT1P5CgArg10CnvFrmType12PreserveType"
			 in module "f_stag" in file "../../lib/util.olb"
   4C04D       83
	  "_FeigSolv__5FstagFPP6VectorP6rfloatT2PiPP6rfloatP6EigArg10CnvFrmType"
			 in module "f_stag" in file "../../lib/util.olb"
   4C0D0       82	"_SetPhi__5FstagFP6VectorN216rfloat" in module "f_stag"
			in file "../../lib/util.olb"
   4C152       9A	"_FforceSite__5FstagFR6MatrixP6VectorPii" in module
			"f_stag" in file "../../lib/util.olb"
   4C1EC       C0	"_EvolveMomFforce__5FstagFP6MatrixP6Vector6rfloatT3" in
			module "f_stag" in file "../../lib/util.olb"
   4C2AC       22	"_FhamiltonNode__5FstagFP6VectorT1" in module "f_stag"
			in file "../../lib/util.olb"
   4C2CE       7A	"_BhamiltonNode__5FstagFP6Vector6rfloat" in module
			"f_stag" in file "../../lib/util.olb"
   4C348        F	"clover_mat_mlt__FPfPCfT2i" in module "f_clover" in
			file "../../lib/util.olb"
   4C357        C	"__as__6rfloatFf" in module "f_clover" in file
			"../../lib/util.olb"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  19
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

   4C363        A	"__opf__6rfloatCFv" in module "f_clover" in file
			"../../lib/util.olb"
   4C36D        D	"__amu__6rfloatFRC6rfloat" in module "f_clover" in
			file "../../lib/util.olb"
   4C37A        E	"__as__6MatrixFRC6Matrix" in module "f_clover" in file
			"../../lib/util.olb"
   4C388        E	"__apl__6MatrixFRC6Matrix" in module "f_clover" in
			file "../../lib/util.olb"
   4C396        E	"__ami__6MatrixFRC6Matrix" in module "f_clover" in
			file "../../lib/util.olb"
   4C3A4        E	"__amu__6MatrixFf" in module "f_clover" in file
			"../../lib/util.olb"
   4C3B2        D	"DotMEqual__6MatrixFRC6MatrixT1" in module "f_clover"
			in file "../../lib/util.olb"
   4C3BF        B	"Dagger__6MatrixFRC6Matrix" in module "f_clover" in
			file "../../lib/util.olb"
   4C3CA       1E	"__dt__6MatrixFv" in module "f_clover" in file
			"../../lib/util.olb"
   4C3E8        D	"CopyVec__6VectorFPC6Vectori" in module "f_clover" in
			file "../../lib/util.olb"
   4C3F5       11	"NormSqNode__6VectorFi" in module "f_clover" in file
			"../../lib/util.olb"
   4C406       11	"ReDotProductNode__6VectorFPC6Vectori" in module
			"f_clover" in file "../../lib/util.olb"
   4C417       1B	"__ct__5CgArgFv" in module "f_clover" in file
			"../../lib/util.olb"
   4C432       1A	"__dt__5CgArgFv" in module "f_clover" in file
			"../../lib/util.olb"
   4C44C        6	"Freflex__7LatticeFP6VectorT1" in module "f_clover" in
			file "../../lib/util.olb"
   4C452      A5F
		   "_EvolveMomFforceSupp__7FcloverFP6MatrixP6VectorN326rfloatT6"
			 in module "f_clover" in file "../../lib/util.olb"
   4CEB1       BF	"___ct__7FcloverFv" in module "f_clover" in file
			"../../lib/util.olb"
   4CF70       71	"___dt__7FcloverFv" in module "f_clover" in file
			"../../lib/util.olb"
   4CFE1        7	"_Fclass__7FcloverFv" in module "f_clover" in file
			"../../lib/util.olb"
   4CFE8       5A	"_FsiteOffsetChkb__7FcloverCFPCi" in module "f_clover"
			in file "../../lib/util.olb"
   4D042       32	"_FsiteOffset__7FcloverCFPCi" in module "f_clover" in
			file "../../lib/util.olb"
   4D074        7	"_ExactFlavors__7FcloverFv" in module "f_clover" in
			file "../../lib/util.olb"
   4D07B        7	"_SpinComponents__7FcloverFv" in module "f_clover" in
			file "../../lib/util.olb"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  20
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

   4D082       1D	"_FsiteSize__7FcloverFv" in module "f_clover" in file
			"../../lib/util.olb"
   4D09F        7	"_FchkbEvl__7FcloverFv" in module "f_clover" in file
			"../../lib/util.olb"
   4D0A6       34
		   "_FmatEvlInv__7FcloverFP6VectorT1P5CgArgP6rfloat10CnvFrmType"
			 in module "f_clover" in file "../../lib/util.olb"
   4D0DA       1C	"_FmatEvlInv__7FcloverFP6VectorT1P5CgArg10CnvFrmType"
			in module "f_clover" in file "../../lib/util.olb"
   4D0F6       36
	"_FmatInv__7FcloverFP6VectorT1P5CgArgP6rfloat10CnvFrmType12PreserveType"
			 in module "f_clover" in file "../../lib/util.olb"
   4D12C       1E
		"_FmatInv__7FcloverFP6VectorT1P5CgArg10CnvFrmType12PreserveType"
			 in module "f_clover" in file "../../lib/util.olb"
   4D14A      14D
	"_FeigSolv__7FcloverFPP6VectorP6rfloatT2PiPP6rfloatP6EigArg10CnvFrmType"
			 in module "f_clover" in file "../../lib/util.olb"
   4D297       B6	"_SetPhi__7FcloverFP6VectorN216rfloat" in module
			"f_clover" in file "../../lib/util.olb"
   4D34D      5B5	"_EvolveMomFforce__7FcloverFP6MatrixP6Vector6rfloatT3"
			in module "f_clover" in file "../../lib/util.olb"
   4D902       6A	"_FhamiltonNode__7FcloverFP6VectorT1" in module
			"f_clover" in file "../../lib/util.olb"
   4D96C       F5	"_BhamiltonNode__7FcloverFP6Vector6rfloat" in module
			"f_clover" in file "../../lib/util.olb"
   4DA61        F	"MatPc__13DiracOpCloverFP6VectorT1" in module "f_clover"
			 in file "../../lib/util.olb"
   4DA70        F	"MatPcDag__13DiracOpCloverFP6VectorT1" in module
			"f_clover" in file "../../lib/util.olb"
   4DA7F        A	"XnodeSites__18GlobalJobParameterFv" in module
			"f_clover" in file "../../lib/util.olb"
   4DA89        A	"YnodeSites__18GlobalJobParameterFv" in module
			"f_clover" in file "../../lib/util.olb"
   4DA93        A	"ZnodeSites__18GlobalJobParameterFv" in module
			"f_clover" in file "../../lib/util.olb"
   4DA9D        A	"TnodeSites__18GlobalJobParameterFv" in module
			"f_clover" in file "../../lib/util.olb"
   4DAA7        A	"VolNodeSites__18GlobalJobParameterFv" in module
			"f_clover" in file "../../lib/util.olb"
   4DAB1        A	"XnodeBc__18GlobalJobParameterFv" in module "f_clover"
			in file "../../lib/util.olb"
   4DABB        A	"YnodeBc__18GlobalJobParameterFv" in module "f_clover"
			in file "../../lib/util.olb"
   4DAC5        A	"ZnodeBc__18GlobalJobParameterFv" in module "f_clover"
			in file "../../lib/util.olb"
   4DACF        A	"TnodeBc__18GlobalJobParameterFv" in module "f_clover"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  21
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

			in file "../../lib/util.olb"
   4DAD9        C	"CloverCoeff__18GlobalJobParameterFv" in module
			"f_clover" in file "../../lib/util.olb"
   4DAE5       32	"___ct__5ErrorFv" in module "error" in file
			"../../lib/util.olb"
   4DB17       11	"___dt__5ErrorFv" in module "error" in file
			"../../lib/util.olb"
   4DB28       49	"_Pointer__5ErrorFPcN21" in module "error" in file
			"../../lib/util.olb"
   4DB71       49	"_FileA__5ErrorFPcN21" in module "error" in file
			"../../lib/util.olb"
   4DBBA       45	"_NotImplemented__5ErrorFPcT1" in module "error" in
			file "../../lib/util.olb"
   4DBFF       5F	"_NotImplemented__5ErrorFPcT1PCce" in module "error" in
			file "../../lib/util.olb"
   4DC5E       5F	"_General__5ErrorFPcT1PCce" in module "error" in file
			"../../lib/util.olb"
   4DCBD       A6	"_MultStagPhases__FP16ConvertArgStruct" in module
			"convert_func" in file "../../lib/util.olb"
   4DD63      11A	"_RunGConverter__FP16ConvertArgStructPUiT2" in module
			"convert_func" in file "../../lib/util.olb"
   4DE7D      240	"_CanonToAnything__FP16ConvertArgStruct10StrOrdType" in
			module "convert_func" in file "../../lib/util.olb"
   4E0BD      14F	"_FcanonToWilson__FP16ConvertArgStruct" in module
			"convert_func" in file "../../lib/util.olb"
   4E20C      13E	"_FwilsonToCanon__FP16ConvertArgStruct" in module
			"convert_func" in file "../../lib/util.olb"
   4E34A       76	"_Convert__7LatticeF10StrOrdTypeP6VectorT2" in module
			"convert" in file "../../lib/util.olb"
   4E3C0      327	"_Convert__7LatticeF10StrOrdType" in module "convert"
			in file "../../lib/util.olb"
   4E6E7       1E	"_Fconvert__5FstagFP6Vector10StrOrdTypeT2" in module
			"convert" in file "../../lib/util.olb"
   4E705       A2	"_Fconvert__7FcloverFP6Vector10StrOrdTypeT2" in module
			"convert" in file "../../lib/util.olb"
   4E7A7       12	"_negate_link" in module "common" in file
			"../../lib/util.olb"
   4E7B9       12	"_site2cram" in module "common" in file
			"../../lib/util.olb"
   4E7CB       14	"_site2dram" in module "common" in file
			"../../lib/util.olb"
   4E7DF       42	"_StagSCUSetup" in module "stag_scu" in file
			"../../lib/d_op_stag_opt.olb"
   4E821       3A	"_StagSCUCommForward" in module "stag_scu" in file
			"../../lib/d_op_stag_opt.olb"
   4E85B       3A	"_StagSCUCommBackward" in module "stag_scu" in file
			"../../lib/d_op_stag_opt.olb"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  22
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

   4E895        2	"_StagSCUComplete" in module "stag_scu" in file
			"../../lib/d_op_stag_opt.olb"
   4E897       11	"Ini$stag_scu_C" in module "stag_scu" in file
			"../../lib/d_op_stag_opt.olb"
   4E8A8       14	"Fin$stag_scu_C" in module "stag_scu" in file
			"../../lib/d_op_stag_opt.olb"
   4E8BC       27	"_destroy_dirac_buf__Fv" in module "dirac_init" in
			file "../../lib/d_op_stag_opt.olb"
   4E8E3     1338	"_dirac_init__FPCv" in module "dirac_init" in file
			"../../lib/d_op_stag_opt.olb"
   4FC1B      5A9	"comm_init__Fi" in module "dirac_init" in file
			"../../lib/d_op_stag_opt.olb"
   501C4       97
		   "___ct__11DiracOpStagFR7LatticeP6VectorT2P5CgArg10CnvFrmType"
			 in module "d_op_stag" in file
			"../../lib/d_op_stag_opt.olb"
   5025B       5F	"___dt__11DiracOpStagFv" in module "d_op_stag" in file
			"../../lib/d_op_stag_opt.olb"
   502BA       2B	"_DiracArg__11DiracOpStagFP5CgArg" in module "d_op_stag"
			 in file "../../lib/d_op_stag_opt.olb"
   502E5       57	"_MatPcDagMatPc__11DiracOpStagFP6VectorT1P6rfloat" in
			module "d_op_stag" in file "../../lib/d_op_stag_opt.olb"
   5033C       2D	"_Dslash__11DiracOpStagFP6VectorT18ChkbType7DagType" in
			module "d_op_stag" in file "../../lib/d_op_stag_opt.olb"
   50369       C2
		       "_MatInv__11DiracOpStagFP6VectorT1P6rfloat12PreserveType"
			 in module "d_op_stag" in file
			"../../lib/d_op_stag_opt.olb"
   5042B       1A	"_MatInv__11DiracOpStagFP6VectorT112PreserveType" in
			module "d_op_stag" in file "../../lib/d_op_stag_opt.olb"
   50445       1A	"_MatInv__11DiracOpStagFP6rfloat12PreserveType" in
			module "d_op_stag" in file "../../lib/d_op_stag_opt.olb"
   5045F       1A	"_MatInv__11DiracOpStagF12PreserveType" in module
			"d_op_stag" in file "../../lib/d_op_stag_opt.olb"
   50479       14	"_RitzEigMat__11DiracOpStagFP6VectorT1" in module
			"d_op_stag" in file "../../lib/d_op_stag_opt.olb"
   5048D       75	"_RitzMat__11DiracOpStagFP6VectorT1" in module
			"d_op_stag" in file "../../lib/d_op_stag_opt.olb"
   50502      47D
		      "_RitzEig__16DiracOpStagTypesFPP6VectorP6rfloatPiP6EigArg"
			 in module "eigen_stag" in file
			"../../lib/d_op_stag_types.olb"
   5097F       1B	"_RitzLatSize__16DiracOpStagTypesFv" in module
			"eigen_stag" in file "../../lib/d_op_stag_types.olb"
   5099A       39
	      "___ct__16DiracOpStagTypesFR7LatticeP6VectorT2P5CgArg10CnvFrmType"
			 in module "d_op_stag_types" in file

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  23
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

			"../../lib/d_op_stag_types.olb"
   509D3       2D	"___dt__16DiracOpStagTypesFv" in module
			"d_op_stag_types" in file
			"../../lib/d_op_stag_types.olb"
   50A00       23	"iDotM__FR6Matrixi" in module "d_op_clover_supp" in
			file "../../lib/d_op_clover.olb"
   50A23       17	"Ini$d_op_clover_supp_C" in module "d_op_clover_supp"
			in file "../../lib/d_op_clover.olb"
   50A3A       18	"Fin$d_op_clover_supp_C" in module "d_op_clover_supp"
			in file "../../lib/d_op_clover.olb"
   50A52       25	"__dt__6MatrixFv" in module "d_op_clover_supp" in file
			"../../lib/d_op_clover.olb"
   50A77       ED	"_GetLink__13DiracOpCloverCFPCii" in module
			"d_op_clover_supp" in file "../../lib/d_op_clover.olb"
   50B64      183	"_SiteFuv__13DiracOpCloverCFR6MatrixPCiiT3" in module
			"d_op_clover_supp" in file "../../lib/d_op_clover.olb"
   50CE7      1AD	"_SiteCloverMat__13DiracOpCloverCFPCiPf" in module
			"d_op_clover_supp" in file "../../lib/d_op_clover.olb"
   50E94       8E	"_CloverMatChkb__13DiracOpCloverCF8ChkbTypei" in
			module "d_op_clover_supp" in file
			"../../lib/d_op_clover.olb"
   50F22       7B
		 "___ct__13DiracOpCloverFR7LatticeP6VectorT2P5CgArg10CnvFrmType"
			 in module "d_op_clover" in file
			"../../lib/d_op_clover.olb"
   50F9D       5A	"___dt__13DiracOpCloverFv" in module "d_op_clover" in
			file "../../lib/d_op_clover.olb"
   50FF7       8C	"_MatPcDagOrNot__13DiracOpCloverCFP6VectorPC6Vectori"
			in module "d_op_clover" in file
			"../../lib/d_op_clover.olb"
   51083      158	"_DiracArg__13DiracOpCloverFP5CgArg" in module
			"d_op_clover" in file "../../lib/d_op_clover.olb"
   511DB       5A	"_MatPcDagMatPc__13DiracOpCloverFP6VectorT1P6rfloat" in
			module "d_op_clover" in file "../../lib/d_op_clover.olb"
   51235       14	"_Dslash__13DiracOpCloverFP6VectorT18ChkbType7DagType"
			in module "d_op_clover" in file
			"../../lib/d_op_clover.olb"
   51249       9A
		     "_MatInv__13DiracOpCloverFP6VectorT1P6rfloat12PreserveType"
			 in module "d_op_clover" in file
			"../../lib/d_op_clover.olb"
   512E3       1A	"_MatInv__13DiracOpCloverFP6VectorT112PreserveType" in
			module "d_op_clover" in file "../../lib/d_op_clover.olb"
   512FD       1A	"_MatInv__13DiracOpCloverFP6rfloat12PreserveType" in
			module "d_op_clover" in file "../../lib/d_op_clover.olb"
   51317       1A	"_MatInv__13DiracOpCloverF12PreserveType" in module
			"d_op_clover" in file "../../lib/d_op_clover.olb"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  24
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

   51331       5D	"_MatEvlInv__13DiracOpCloverFP6VectorT1P6rfloat" in
			module "d_op_clover" in file "../../lib/d_op_clover.olb"
   5138E        D	"_MatEvlInv__13DiracOpCloverFP6rfloat" in module
			"d_op_clover" in file "../../lib/d_op_clover.olb"
   5139B       6E	"_MatHerm__13DiracOpCloverFP6VectorT1" in module
			"d_op_clover" in file "../../lib/d_op_clover.olb"
   51409       14	"_Mat__13DiracOpCloverFP6VectorT1" in module
			"d_op_clover" in file "../../lib/d_op_clover.olb"
   5141D       14	"_MatDag__13DiracOpCloverFP6VectorT1" in module
			"d_op_clover" in file "../../lib/d_op_clover.olb"
   51431       B3	"_CalcHmdForceVecs__13DiracOpCloverFP6Vector" in
			module "d_op_clover" in file "../../lib/d_op_clover.olb"
   514E4       70	"_clover_init__FP6Clover" in module "clover" in file
			"../../lib/d_op_clover.olb"
   51554       2F	"_clover_end__FP6Clover" in module "clover" in file
			"../../lib/d_op_clover.olb"
   51583       1E	"_mat_hrm_cmpr" in module "lapack" in file
			"../../lib/lapack.olb"
   515A1       59	"_mat_hrm_decm" in module "lapack" in file
			"../../lib/lapack.olb"
   515FA      15E	"_mat_hrm_ldl" in module "lapack" in file
			"../../lib/lapack.olb"
   51758      180	"mat_inv_ldl_cmpr__FP6rfloatPC6rfloati" in module
			"lapack" in file "../../lib/lapack.olb"
   518D8      15A	"_mat_inv" in module "lapack" in file
			"../../lib/lapack.olb"
   51A32       52	".text" in module "wilson_dslash" in file
			"../../lib/d_op_wilson_opt_nos.olb"
   51A84       1C	"ram1" in module "wfm_spproj_segment" in file
			"../../lib/d_op_wilson_opt_nos.olb"
   51AA0       3F	".text" in module "wfm_scu_wait" in file
			"../../lib/d_op_wilson_opt_nos.olb"
   51ADF       35	".text" in module "wfm_scu_init" in file
			"../../lib/d_op_wilson_opt_nos.olb"
   51B20       40	".text" in module "wfm_mat_trick" in file
			"../../lib/d_op_wilson_opt_nos.olb"
   51B60       2B	".text" in module "wfm_dslash" in file
			"../../lib/d_op_wilson_opt_nos.olb"
   51B8B       57	".text" in module "wfm_comm_forward" in file
			"../../lib/d_op_wilson_opt_nos.olb"
   51BE2       58	".text" in module "wfm_comm_backward" in file
			"../../lib/d_op_wilson_opt_nos.olb"
   51C40       40	".text" in module "wfm_cmat_spproj" in file
			"../../lib/d_op_wilson_opt_nos.olb"
   51C80       10	".text" in module "wfm_cb_reg0" in file
			"../../lib/d_op_wilson_opt_nos.olb"
   51C90      108	"_wilson_init__FP6Wilson" in module "wilson_init" in

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  25
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

			file "../../lib/d_op_wilson_opt_nos.olb"
   51D98       70	"_wilson_end__FP6Wilson" in module "wilson_end" in
			file "../../lib/d_op_wilson_opt_nos.olb"
   51E08       1C	"ceiling__Fi" in module "wfm_sublatt_pointers" in file
			"../../lib/d_op_wilson_opt_nos.olb"
   51E24       1C	"floor__Fi" in module "wfm_sublatt_pointers" in file
			"../../lib/d_op_wilson_opt_nos.olb"
   51E40      32D	"_wfm_sublatt_pointers__FiN41P6Wilson" in module
			"wfm_sublatt_pointers" in file
			"../../lib/d_op_wilson_opt_nos.olb"
   5216D       46	"_wfm_copy_forward" in module "wfm_copy_forward" in
			file "../../lib/d_op_wilson_opt_nos.olb"
   521B3       46	"_wfm_copy_backward" in module "wfm_copy_backward" in
			file "../../lib/d_op_wilson_opt_nos.olb"
   521F9       6F	"_GramSchm__FPP6VectoriT1N22" in module "ritz" in file
			"../../lib/d_op_base.olb"
   52268      BF4	"_Ritz__7DiracOpFPP6VectoriR6rfloat6rfloatN34N42T4N22"
			in module "ritz" in file "../../lib/d_op_base.olb"
   52E5C      94A
		       "_Jacobi__7DiracOpFPP6VectoriP6rfloatP8Rcomplex6rfloatT2"
			 in module "jacobi" in file "../../lib/d_op_base.olb"
   537A6      48A	"_InvCg__7DiracOpFP6VectorT16rfloatP6rfloat" in module
			"inv_cg" in file "../../lib/d_op_base.olb"
   53C30       23	"_InvCg__7DiracOpFP6VectorT1P6rfloat" in module "inv_cg"
			 in file "../../lib/d_op_base.olb"
   53C53       23	"_InvCg__7DiracOpFP6rfloat" in module "inv_cg" in file
			"../../lib/d_op_base.olb"
   53C76       A9	"BondCond__FR7LatticeP6Matrix" in module "dirac_op_base"
			 in file "../../lib/d_op_base.olb"
   53D1F        B	"_DiracOpGlbSum__7DiracOpFP6rfloat" in module
			"dirac_op_base" in file "../../lib/d_op_base.olb"
   53D2A      142	"___ct__7DiracOpFR7LatticeP6VectorT2P5CgArg10CnvFrmType"
			 in module "dirac_op_base" in file
			"../../lib/d_op_base.olb"
   53E6C       9C	"___dt__7DiracOpFv" in module "dirac_op_base" in file
			"../../lib/d_op_base.olb"
   53F08      7CB
		    "_RitzEig__18DiracOpWilsonTypesFPP6VectorP6rfloatPiP6EigArg"
			 in module "eigen_wilson" in file
			"../../lib/d_op_wilson_types.olb"
   546D3       3D	"_MultGamma__18DiracOpWilsonTypesFP6VectorPC6VectoriT3"
			in module "d_op_wilson_types" in file
			"../../lib/d_op_wilson_types.olb"
   54710       39
	    "___ct__18DiracOpWilsonTypesFR7LatticeP6VectorT2P5CgArg10CnvFrmType"
			 in module "d_op_wilson_types" in file
			"../../lib/d_op_wilson_types.olb"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  26
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

   54749       2D	"___dt__18DiracOpWilsonTypesFv" in module
			"d_op_wilson_types" in file
			"../../lib/d_op_wilson_types.olb"
   54776       69	"_MatDagMat__18DiracOpWilsonTypesFP6VectorT1" in
			module "d_op_wilson_types" in file
			"../../lib/d_op_wilson_types.olb"
   547DF       6F	"_RitzEigMat__18DiracOpWilsonTypesFP6VectorT1" in
			module "d_op_wilson_types" in file
			"../../lib/d_op_wilson_types.olb"
   5484E       5F	"_RitzMat__18DiracOpWilsonTypesFP6VectorT1" in module
			"d_op_wilson_types" in file
			"../../lib/d_op_wilson_types.olb"
   548AD       4F	"_RitzLatSize__18DiracOpWilsonTypesFv" in module
			"d_op_wilson_types" in file
			"../../lib/d_op_wilson_types.olb"
   548FC        9	".text" in module "double64_rnd" in file
			"../../lib/glb_cpp.olb"
   54905      128	"_glb_sum_five__FP6rfloat" in module "glb_sum_five" in
			file "../../lib/glb_cpp.olb"
   54A2D        C	"Ini$glb_sum_five_C" in module "glb_sum_five" in file
			"../../lib/glb_cpp.olb"
   54A39       A1	"_glb_sum__FP6rfloat" in module "glb_sum" in file
			"../../lib/glb_cpp.olb"
   54ADA        C	"Ini$glb_sum_C" in module "glb_sum" in file
			"../../lib/glb_cpp.olb"
   54AE6       66	"_getPlusData__FPfT1iT3" in module "get_data" in file
			"../../lib/scu.olb"
   54B4C       66	"_getMinusData__FPfT1iT3" in module "get_data" in file
			"../../lib/scu.olb"
   54BB2       3D	"_p2vVector__Fv" in module "p2v" in file
			"../../lib/mem.olb"
   54BEF       66	"_p2vWilsonLib__Fv" in module "p2v" in file
			"../../lib/mem.olb"
   54C55       3D	"_p2vStagDs__Fv" in module "p2v" in file
			"../../lib/mem.olb"
   54C92       5D	"_p2vCloverLib__Fv" in module "p2v" in file
			"../../lib/mem.olb"
   54CEF       1B	"_saveCbufCntrlReg__Fv" in module "cbuf" in file
			"../../lib/cbuf.olb"
   54D0A       10	"_restoreCbufCntrlReg__Fv" in module "cbuf" in file
			"../../lib/cbuf.olb"
   54D1A        A	"_setCbufCntrlReg__FiUi" in module "cbuf" in file
			"../../lib/cbuf.olb"
   54D24        C	"___ct__9SCUDirArgFv" in module "scu_dir_arg" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54D30       23	"___ct__9SCUDirArgFPv6SCUDir5SCUXRiN24" in module
			"scu_dir_arg" in file

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  27
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54D53        C	"___dt__9SCUDirArgFv" in module "scu_dir_arg" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54D5F       31	"_Init__9SCUDirArgFPv6SCUDir5SCUXRiN24" in module
			"scu_dir_arg" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54D90       14	"_Reload__9SCUDirArgFPviN22" in module "scu_dir_arg" in
			file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54DA4       24	".text" in module "sysfunc" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54DC8       25	".text" in module "tclowio" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54DED        8	".text" in module "stub64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54DF5       21	"_vsprintf" in module "stdio" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54E16       19	"PutAChar" in module "mypf" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54E2F       17	"PadOut" in module "mypf" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54E46       1B	"PutChars" in module "mypf" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54E61       29	"WidenOut" in module "mypf" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54E8A       17	"Decimal_WP" in module "mypf" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54EA1       3A	"ConvertInt" in module "mypf" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54EDB       61	"WriteIntField" in module "mypf" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54F3C       27	"PutImage" in module "mypf" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54F63       29	"RoundOff" in module "mypf" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54F8C       1C	"CountZeros" in module "mypf" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   54FA8       A1	"WriteFfloat" in module "mypf" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   55049       A0	"WriteEfloat" in module "mypf" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   550E9      266	"__pf" in module "mypf" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   5534F       E0	"__RealConvF64" in module "cvtf64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   5542F       5E	"__TimesPower" in module "cvtdbl" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  28
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

   5548D       81	"__ConvertDouble" in module "cvtdbl" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcio30bs.olb"
   5550E       4E	".text" in module "qcdsp_tcroot" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   5555C       20	"_strncpy" in module "strncpy" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   5557C        B	"_strlen" in module "strlen" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55587       1B	".text" in module "memmov" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   555A2        2	"__array_pointer_not_from_vec_new" in module "vecnwdl"
			in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   555A4       3E	"___vec_new" in module "vecnwdl" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   555E2       49	"___vec_delete" in module "vecnwdl" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   5562B       26	"dbgrpc" in module "tcrpc" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55651        7	".text" in module "tcrpc" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55658       32	".text" in module "tcinit" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   5568A        2	"___pure_virtual_called" in module "purevirt" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   5568C       21	"___nw__FUi" in module "new" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   556AD        9	"___dl__FPv" in module "new" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   556B6       25	"minsert" in module "memory" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   556DB       21	"mremove" in module "memory" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   556FC       11	"_minit" in module "memory" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   5570D       48	"_malloc" in module "memory" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55755       76	"_free" in module "memory" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   557CB        C	"pfuncs" in module "mathtsup" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   557D7       22	"_exp__FRC8double64" in module "math64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   557F9       22	"_sqrt__FRC8double64" in module "math64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   5581B      2BD	"Ini$shared_math64_cpp" in module "math64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  29
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

   55AD8        6	"_abort" in module "inifin" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55ADE       22	"_exit" in module "inifin" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55B00        2	"__main" in module "inifin" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55B02       11	".text" in module "clock" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55B13       15	".text" in module "trunci32f64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55B28        D	".text" in module "subf32f40" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55B35       56	".text" in module "sqrf64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55B8B        B	".text" in module "roundi32f64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55B96        A	".text" in module "neqf64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55BA0        8	".text" in module "mpyf64f32" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55BA8       70	".text" in module "mpyf64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55C18       A4	".text" in module "mpyf40" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55CBC        D	".text" in module "lssf64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55CC9       3A	".text" in module "fltsup" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55D03        F	".text" in module "floori32f64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55D12        D	".text" in module "floorf40f64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55D1F      16C	".text" in module "f64sqrt" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   55E8B      433	".text" in module "f64exp" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   562BE       5A	".text" in module "divf64" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56318        B	".text" in module "convf64i32" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56323        A	".text" in module "convf64f40" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   5632D        8	".text" in module "convf64f32" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56335        7	".text" in module "convf40f32" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   5633C      13A	".text" in module "addf64" in file

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  30
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56476       1E	".text" in module "addf40" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56494       36	"SQRT32" in module "trig" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   564CA       6E	"ALOG32" in module "trig" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56538        C	"DRU32Z" in module "uns_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56544        D	"MU32U" in module "uns_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56551       2F	"DRU32" in module "uns_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56580        8	"SATI" in module "int_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56588        D	"MS32U" in module "int_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56595        8	"C_SATI" in module "int_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   5659D       2B	"DIV_I30" in module "int_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   565C8       13	"RS32C" in module "int_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   565DB       16	"DS32C" in module "int_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   565F1       16	"RS32Z" in module "int_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56607       14	"DS32W" in module "int_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   5661B       3F	"SDIVF40" in module "fp_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   5665A        8	"DVF32W" in module "fp_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56662       20	"INV_F30" in module "fp_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56682       20	"INV32" in module "fp_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   566A2       50	"EF40U" in module "fp_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   566F2       3A	"SMPYF40" in module "fp_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   5672C        A	"SATF" in module "fp_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56736       4A	"EF40" in module "fp_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56780        9	"INV32Z" in module "fp_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  31
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

   56789        8	"SDF40W" in module "fp_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   56791       38	"DVF32U" in module "fp_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   567C9       2C	"INV40" in module "fp_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"
   567F5        A	"C_SATF" in module "fp_math" in file
			"/usr/local/tartan/v2.1/etc/qcdsp_v5.1.5/tcrt30bs.olb"

  Total Allocation = F842 (hex)


		Output Sections

Section List for Module *unnamed* in file
"/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/qcdsp.outtof"

Number	Physical   Start   Length  Kind		Access	Name
  1	    1000     1000    4321D data			".bss"
  2	   4421D    4421D     2DA0 constant		".data"
  3	   46FBD    46FBD     F842 code			".text"
  4	   567FF   809F95       2C code			"vector_overlay"
  5	   5682B   809800      2DC code			"wfm0_overlay"
  6	   56B07   809C00      1A5 code			"wfm1_overlay"
  7	   56CAC   809E00       C1 code			"cfm_overlay"
  8	   56D6D   809800      282 code			"stag_ds_overlay"
  9	   56FEF        0        0 code			"ghb_overlay"
 10	   56FEF        0        0 code			*unnamed*
 11		        0        0 debug		"debug_directives"
 12		        0        0 debugstring		"debug_strings"
 13		        0        0 line number		"debug_source_location"


		Output Symbols

Global Symbol List for Module *unnamed* in file
"/homeq22/vranas/16k/sfw/phys/exc/phys_v3.11.4/tests/f_hmd/qcdsp.outtof"

Kind  Offset Sect Virtual Physical	Name
 W      E552   3    5550F    5550F "$START$"
        321D   1     421D     421D ".stack"
       1321D   1    1421D    1421D ".sysmem"
         198   6   809D98    56C9F "ab"
         199   6   809D99    56CA0 "ab0"
         19D   6   809D9D    56CA4 "ab0_t"
         19A   6   809D9A    56CA1 "ab1"
         19E   6   809D9E    56CA5 "ab1_t"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  32
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

         19B   6   809D9B    56CA2 "ab2"
         19F   6   809D9F    56CA6 "ab2_t"
         19C   6   809D9C    56CA3 "ab3"
         1A0   6   809DA0    56CA7 "ab3_t"
        E814   3    557D1    557D1 "ADDEXP"
        E813   3    557D0    557D0 "ADDEXP$LAJ"
        F37F   3    5633C    5633C "addf64"
        F3B1   3    5636E    5636E "addf64_cont"
          CF   6   809CCF    56BD6 "af"
          D0   6   809CD0    56BD7 "af0"
          D1   6   809CD1    56BD8 "af1"
          D2   6   809CD2    56BD9 "af2"
          D3   6   809CD3    56BDA "af3"
          19   2    44236    44236 "ARcinit$"
 W      2D9A   2    46FB7    46FB7 "ARFin$"
 W      2D8F   2    46FAC    46FAC "ARIni$"
        F55C   3    56519    56519 "ARTALOG1032"
        F55D   3    5651A    5651A "ARTALOG1032$LAJ"
        F530   3    564ED    564ED "ARTALOG232"
        F531   3    564EE    564EE "ARTALOG232$LAJ"
        F510   3    564CD    564CD "ARTALOG32"
        F511   3    564CE    564CE "ARTALOG32$LAJ"
        F594   3    56551    56551 "ARTBDIVREMU32"
        F595   3    56552    56552 "ARTBDIVREMU32$LAJ"
        F57B   3    56538    56538 "ARTBDIVREMU32Z"
        F57C   3    56539    56539 "ARTBDIVREMU32Z$LAJ"
        F64A   3    56607    56607 "ARTBDIVS32UZ"
        F64B   3    56608    56608 "ARTBDIVS32UZ$LAJ"
        F594   3    56551    56551 "ARTBDIVU32"
        F595   3    56552    56552 "ARTBDIVU32$LAJ"
        F57B   3    56538    56538 "ARTBDIVU32Z"
        F57C   3    56539    56539 "ARTBDIVU32Z$LAJ"
        F634   3    565F1    565F1 "ARTBREMS32Z"
        F635   3    565F2    565F2 "ARTBREMS32Z$LAJ"
           2   2    4421F    4421F "ARTDEFPAGE"
        F7D4   3    56791    56791 "ARTDIVF32U"
        F7D5   3    56792    56792 "ARTDIVF32U$LAJ"
        F69D   3    5665A    5665A "ARTDIVF32UZ"
        F69E   3    5665B    5665B "ARTDIVF32UZ$LAJ"
        F594   3    56551    56551 "ARTDIVREMU32"
        F595   3    56552    56552 "ARTDIVREMU32$LAJ"
        F57B   3    56538    56538 "ARTDIVREMU32Z"
        F57C   3    56539    56539 "ARTDIVREMU32Z$LAJ"
        F64A   3    56607    56607 "ARTDIVS32UZ"
        F64B   3    56608    56608 "ARTDIVS32UZ$LAJ"
        F594   3    56551    56551 "ARTDIVU32"
        F595   3    56552    56552 "ARTDIVU32$LAJ"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  33
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        F57B   3    56538    56538 "ARTDIVU32Z"
        F57C   3    56539    56539 "ARTDIVU32Z$LAJ"
        F77A   3    56737    56737 "ARTEXPF40"
        F779   3    56736    56736 "ARTEXPF40$LAJ"
        F6E6   3    566A3    566A3 "ARTEXPF40U"
        F6E5   3    566A2    566A2 "ARTEXPF40U$LAJ"
        E689   3    55646    55646 "ARTHALT"
        F6C5   3    56682    56682 "ARTINVERSE32"
        F6C6   3    56683    56683 "ARTINVERSE32$LAJ"
        F7C3   3    56780    56780 "ARTINVERSE32Z"
        F7C4   3    56781    56781 "ARTINVERSE32Z$LAJ"
        F80C   3    567C9    567C9 "ARTINVERSE40"
        F80D   3    567CA    567CA "ARTINVERSE40$LAJ"
        E552   3    5550F    5550F "ARTMAIN"
        F5CB   3    56588    56588 "ARTMPYS32U"
        F5CC   3    56589    56589 "ARTMPYS32U$LAJ"
        F587   3    56544    56544 "ARTMPYU32U"
        F588   3    56545    56545 "ARTMPYU32U$LAJ"
        E68D   3    5564A    5564A "ARTREAD"
        F634   3    565F1    565F1 "ARTREMS32Z"
        F635   3    565F2    565F2 "ARTREMS32Z$LAJ"
        F76F   3    5672C    5672C "ARTSATF"
        F770   3    5672D    5672D "ARTSATF$LAJ"
        F5C3   3    56580    56580 "ARTSATI"
        F5C4   3    56581    56581 "ARTSATI$LAJ"
        F65E   3    5661B    5661B "ARTSDIVF40"
        F65F   3    5661C    5661C "ARTSDIVF40$LAJ"
        F7CC   3    56789    56789 "ARTSDIVF40UZ"
        F7CD   3    5678A    5678A "ARTSDIVF40UZ$LAJ"
        F757   3    56714    56714 "ARTSMPYF40"
        F756   3    56713    56713 "ARTSMPYF40$LAJ"
        F735   3    566F2    566F2 "ARTSMPYF40U"
        F736   3    566F3    566F3 "ARTSMPYF40U$LAJ"
        F4E8   3    564A5    564A5 "ARTSQRT32"
        F4E9   3    564A6    564A6 "ARTSQRT32$LAJ"
        EB45   3    55B02    55B02 "ARTSTARTTIME30"
        E683   3    55640    55640 "ARTUTRAP"
        E68B   3    55648    55648 "ARTWRITE"
        E693   3    55650    55650 "ART_RESET_RUNS_HERE"
          E5   5   8098E5    56910 "bank0"
          E6   5   8098E6    56911 "bank1a"
          E7   5   8098E7    56912 "bank1b"
          E8   5   8098E8    56913 "bank2a"
          E9   5   8098E9    56914 "bank2b"
          EA   5   8098EA    56915 "bank3a"
          EB   5   8098EB    56916 "bank3b"
          EC   5   8098EC    56917 "bank4a"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  34
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

          ED   5   8098ED    56918 "bank4b"
          EF   5   8098EF    5691A "cb_cntrl_b"
        F838   3    567F5    567F5 "C_SATF"
        F5D9   3    56596    56596 "C_SATI"
         108   6   809D08    56C0F "c_u0"
         106   6   809D06    56C0D "c_u0_p"
         150   6   809D50    56C57 "c_u1"
         107   6   809D07    56C0E "c_u1_p"
           1   2    4421E    4421E "DEFALT"
          E4   5   8098E4    5690F "direct"
        F319   3    562D6    562D6 "divf64"
        F319   3    562D6    562D6 "divf64_cont"
        F5E0   3    5659D    5659D "DIV_I30"
        E69B   3    55658    55658 "do_cinit"
        E6C4   3    55681    55681 "do_cppfin"
        E6BC   3    55679    55679 "do_cppini"
        ED4F   3    55D0C    55D0C "erangei32f64"
        ED50   3    55D0D    55D0D "erangei32f64_cont"
        F4C8   3    56485    56485 "errorf40"
        F3E3   3    563A0    563A0 "errorf64"
        9A7E   3    50A3B    50A3B "Fin$d_op_clover_supp_C"
        9A7D   3    50A3A    50A3A "Fin$d_op_clover_supp_C$LAJ"
        4CA3   3    4BC60    4BC60 "Fin$f_stag_C"
        4CA2   3    4BC5F    4BC5F "Fin$f_stag_C$LAJ"
        3705   3    4A6C2    4A6C2 "Fin$lattice_base_C"
        3704   3    4A6C1    4A6C1 "Fin$lattice_base_C$LAJ"
         3CF   3    4738C    4738C "Fin$main_C"
         3CE   3    4738B    4738B "Fin$main_C$LAJ"
        78EC   3    4E8A9    4E8A9 "Fin$stag_scu_C"
        78EB   3    4E8A8    4E8A8 "Fin$stag_scu_C$LAJ"
        ED5B   3    55D18    55D18 "floorf40f64_cont"
        E80F   3    557CC    557CC "GETEXP"
        E80E   3    557CB    557CB "GETEXP$LAJ"
        9A67   3    50A24    50A24 "Ini$d_op_clover_supp_C"
        9A66   3    50A23    50A23 "Ini$d_op_clover_supp_C$LAJ"
        4C7B   3    4BC38    4BC38 "Ini$f_stag_C"
        4C7A   3    4BC37    4BC37 "Ini$f_stag_C$LAJ"
        DB1E   3    54ADB    54ADB "Ini$glb_sum_C"
        DB1D   3    54ADA    54ADA "Ini$glb_sum_C$LAJ"
        DA71   3    54A2E    54A2E "Ini$glb_sum_five_C"
        DA70   3    54A2D    54A2D "Ini$glb_sum_five_C$LAJ"
        4858   3    4B815    4B815 "Ini$g_wilson_C"
        4859   3    4B816    4B816 "Ini$g_wilson_C$LAJ"
        36C2   3    4A67F    4A67F "Ini$lattice_base_C"
        36C1   3    4A67E    4A67E "Ini$lattice_base_C$LAJ"
         3BD   3    4737A    4737A "Ini$main_C"
         3BC   3    47379    47379 "Ini$main_C$LAJ"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  35
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        33CA   3    4A387    4A387 "Ini$rcomplex_C"
        33C9   3    4A386    4A386 "Ini$rcomplex_C$LAJ"
        E85F   3    5581C    5581C "Ini$shared_math64_cpp"
        E85E   3    5581B    5581B "Ini$shared_math64_cpp$LAJ"
        78DB   3    4E898    4E898 "Ini$stag_scu_C"
        78DA   3    4E897    4E897 "Ini$stag_scu_C$LAJ"
 W      E678   3    55635    55635 "INT$DINT"
 W      E66E   3    5562B    5562B "INT$INT0"
 W      E66F   3    5562C    5562C "INT$INT1"
 W      E670   3    5562D    5562D "INT$INT2"
 W      E671   3    5562E    5562E "INT$INT3"
 W      E673   3    55630    55630 "INT$RINT0"
 W      E675   3    55632    55632 "INT$RINT1"
 W      E676   3    55633    55633 "INT$TINT0"
 W      E677   3    55634    55634 "INT$TINT1"
 W      E686   3    55643    55643 "INT$TRAP12"
 W      E672   3    5562F    5562F "INT$XINT0"
 W      E674   3    55631    55631 "INT$XINT1"
        F6A5   3    56662    56662 "INV_F30"
         1A3   6   809DA3    56CAA "i_kappa_sq"
         1A2   6   809DA2    56CA9 "kappa_sq"
        ECC4   3    55C81    55C81 "mpyf40"
        EC99   3    55C56    55C56 "mpyf40u"
        EBEB   3    55BA8    55BA8 "mpyf64"
        EBFE   3    55BBB    55BBB "mpyf64_cont"
          CD   6   809CCD    56BD4 "mp_sq_p"
         1A4   6   809DA4    56CAB "m_kappa_sq"
        F495   3    56452    56452 "negf64"
        F4C7   3    56484    56484 "overflowf40"
        F3E2   3    5639F    5639F "overflowf64"
         1A1   6   809DA1    56CA8 "parameters"
         26E   8   809A6E    56FDB "pop_reg"
         25B   8   809A5B    56FC8 "push_reg"
          EE   5   8098EE    56919 "scu_b"
        DE17   3    54DD4    54DD4 "SIOWCHAR"
        EBC4   3    55B81    55B81 "sqrf64"
        F421   3    563DE    563DE "subf64"
        F45B   3    56418    56418 "subf64_cont"
          D6   6   809CD6    56BDD "tas0"
          D5   6   809CD5    56BDC "tas0_p"
         2AC   5   809AAC    56AD7 "tas1"
         2AB   5   809AAB    56AD6 "tas1_p"
          D4   6   809CD4    56BDB "tpsi0_p"
         2AA   5   809AAA    56AD5 "tpsi1_p"
        321A   1     421A     421A "TwoTo23"
          CB   6   809CCB    56BD2 "u0"
          CC   6   809CCC    56BD3 "u1"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  36
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        E686   3    55643    55643 "undeftrap"
        ACCC   3    51C89    51C89 "wfm_r_cb0"
        ACC5   3    51C82    51C82 "wfm_s_cb0"
          CE   6   809CCE    56BD5 "wilson_p"
        EB1C   3    55AD9    55AD9 "_abort"
        EB1B   3    55AD8    55AD8 "_abort$LAJ"
        F4AA   3    56467    56467 "_absf64"
        1C84   3    48C41    48C41 "_Active__7VerboseF16VerboseLevelType"
        1C83   3    48C40    48C40 "_Active__7VerboseF16VerboseLevelType$LAJ"
        F38B   3    56348    56348 "_addf32f64"
        F4B9   3    56476    56476 "_addf40"
        F3A5   3    56362    56362 "_addf64"
        F397   3    56354    56354 "_addf64f32"
        1BC5   3    48B82    48B82 "_AlgLattice__3AlgFv"
        1BC6   3    48B83    48B83 "_AlgLattice__3AlgFv$LAJ"
        2759   3    49716    49716 "_AntiHermMatrix__6MatrixFPCf"
        275A   3    49717    49717 "_AntiHermMatrix__6MatrixFPCf$LAJ"
        E691   3    5564E    5564E "_ART_REMOTE_CALL_HANDLER"
        E690   3    5564D    5564D "_ART_REMOTE_CALL_RETURN"
        E68F   3    5564C    5564C "_ART_REMOTE_CALL_START"
        4219   3    4B1D6    4B1D6 "_Aux0Ptr__7LatticeFv"
        421A   3    4B1D7    4B1D7 "_Aux0Ptr__7LatticeFv$LAJ"
        421E   3    4B1DB    4B1DB "_Aux1Ptr__7LatticeFv"
        421F   3    4B1DC    4B1DC "_Aux1Ptr__7LatticeFv$LAJ"
        5312   3    4C2CF    4C2CF "_BhamiltonNode__5FstagFP6Vector6rfloat"
        5311   3    4C2CE    4C2CE "_BhamiltonNode__5FstagFP6Vector6rfloat$LAJ"
        69B0   3    4D96D    4D96D "_BhamiltonNode__7FcloverFP6Vector6rfloat"
        69AF   3    4D96C    4D96C
				  "_BhamiltonNode__7FcloverFP6Vector6rfloat$LAJ"
        A475   3    51432    51432 "_CalcHmdForceVecs__13DiracOpCloverFP6Vector"
        A474   3    51431    51431
			       "_CalcHmdForceVecs__13DiracOpCloverFP6Vector$LAJ"
        6EC1   3    4DE7E    4DE7E
			    "_CanonToAnything__FP16ConvertArgStruct10StrOrdType"
        6EC0   3    4DE7D    4DE7D
			"_CanonToAnything__FP16ConvertArgStruct10StrOrdType$LAJ"
           0   7   809E00    56CAC "_cfm_dest"
        EB4E   3    55B0B    55B0B "_clock"
        1F08   3    48EC5    48EC5 "_Clock__7VerboseFPcT1PCce"
        1F07   3    48EC4    48EC4 "_Clock__7VerboseFPcT1PCce$LAJ"
        9ED8   3    50E95    50E95 "_CloverMatChkb__13DiracOpCloverCF8ChkbTypei"
        9ED7   3    50E94    50E94
			       "_CloverMatChkb__13DiracOpCloverCF8ChkbTypei$LAJ"
        2D0A   2    46F27    46F27 "_clover_cram_scratch_addr"
        A598   3    51555    51555 "_clover_end__FP6Clover"
        A597   3    51554    51554 "_clover_end__FP6Clover$LAJ"
        A528   3    514E5    514E5 "_clover_init__FP6Clover"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  37
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        A527   3    514E4    514E4 "_clover_init__FP6Clover$LAJ"
           0   7   809E00    56CAC "_clover_mat_mlt_asm"
        1ECA   1     2ECA     2ECA "_cname_none"
        3A5C   3    4AA19    4AA19 "_Colors__7LatticeFv"
        3A5B   3    4AA18    4AA18 "_Colors__7LatticeFv$LAJ"
        202D   3    48FEA    48FEA "_compDotProduct"
        202C   3    48FE9    48FE9 "_compDotProduct$LAJ"
        21EE   3    491AB    491AB "_CompDotProductGlbSum__6VectorFPC6Vectori"
        21ED   3    491AA    491AA
				 "_CompDotProductGlbSum__6VectorFPC6Vectori$LAJ"
        3395   3    4A352    4A352 "_conj__FRC8Rcomplex"
        3394   3    4A351    4A351 "_conj__FRC8Rcomplex$LAJ"
        1FE5   1     2FE5     2FE5 "_converting_str"
        7404   3    4E3C1    4E3C1 "_Convert__7LatticeF10StrOrdType"
        7403   3    4E3C0    4E3C0 "_Convert__7LatticeF10StrOrdType$LAJ"
        738E   3    4E34B    4E34B "_Convert__7LatticeF10StrOrdTypeP6VectorT2"
        738D   3    4E34A    4E34A
				 "_Convert__7LatticeF10StrOrdTypeP6VectorT2$LAJ"
        F378   3    56335    56335 "_convf40f32"
        F370   3    5632D    5632D "_convf64f32"
        F366   3    56323    56323 "_convf64f40"
        F35B   3    56318    56318 "_convf64i32"
        DDEA   3    54DA7    54DA7 "_CoorT"
        DDEB   3    54DA8    54DA8 "_CoorX"
        DDEC   3    54DA9    54DA9 "_CoorY"
        DDED   3    54DAA    54DAA "_CoorZ"
        3A31   3    4A9EE    4A9EE "_CopyGaugeField__7LatticeFP6Matrix"
        3A30   3    4A9ED    4A9ED "_CopyGaugeField__7LatticeFP6Matrix$LAJ"
        2713   3    496D0    496D0 "_Cross2__6MatrixFRC6VectorT1"
        2712   3    496CF    496CF "_Cross2__6MatrixFRC6VectorT1$LAJ"
        23C4   3    49381    49381 "_crossProductConj__FP6rfloatPC6rfloatT2"
        23C3   3    49380    49380 "_crossProductConj__FP6rfloatPC6rfloatT2$LAJ"
        205A   3    49017    49017 "_cTimesV1PlusV2"
        205B   3    49018    49018 "_cTimesV1PlusV2$LAJ"
 W      E552   3    5550F    5550F "_c_int00"
        302A   1     402A     402A "_d64Cube_Root_2"
        3030   1     4030     4030 "_d64Cube_Root_3"
        301C   1     401C     401C "_d64e"
        3004   1     4004     4004 "_d64Eps"
        3006   1     4006     4006 "_d64EpsNeg"
        3024   1     4024     4024 "_d64e_To_Half_Pi"
        3026   1     4026     4026 "_d64e_To_Minus_Half_Pi"
        3022   1     4022     4022 "_d64e_To_Minus_Pi"
        3020   1     4020     4020 "_d64e_To_Pi"
        3036   1     4036     4036 "_d64First"
        300A   1     400A     400A "_d64Integer_First"
        300C   1     400C     400C "_d64Integer_Last"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  38
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        2FFE   1     3FFE     3FFE "_d64Last"
        302C   1     402C     402C "_d64Ln_2"
        3032   1     4032     4032 "_d64Ln_3"
        302E   1     402E     402E "_d64Log_2"
        3034   1     4034     4034 "_d64Log_3"
        301E   1     401E     401E "_d64One_Over_e"
        3012   1     4012     4012 "_d64One_Over_Pi"
        300E   1     400E     400E "_d64Pi"
        301A   1     401A     401A "_d64Pi_Over_Four"
        3018   1     4018     4018 "_d64Pi_Over_Three"
        3016   1     4016     4016 "_d64Pi_Over_Two"
        3002   1     4002     4002 "_d64SmallNeg"
        3000   1     4000     4000 "_d64SmallPos"
        3028   1     4028     4028 "_d64Sqrt_2"
        3014   1     4014     4014 "_d64Two_Over_Pi"
        3010   1     4010     4010 "_d64Two_Pi"
        3008   1     4008     4008 "_d64Xmin"
        278F   3    4974C    4974C "_Dagger__6MatrixFPCf"
        2790   3    4974D    4974D "_Dagger__6MatrixFPCf$LAJ"
        DDE8   3    54DA5    54DA5 "_DbNum"
        1E8E   3    48E4B    48E4B "_Debug__7VerboseFPCce"
        1E8D   3    48E4A    48E4A "_Debug__7VerboseFPCce$LAJ"
        1E65   3    48E22    48E22 "_Debug__7VerboseFPcT1PCce"
        1E64   3    48E21    48E21 "_Debug__7VerboseFPcT1PCce$LAJ"
        ED0F   3    55CCC    55CCC "_defloatf32"
        ED19   3    55CD6    55CD6 "_defloatf40"
        7900   3    4E8BD    4E8BD "_destroy_dirac_buf__Fv"
        78FF   3    4E8BC    4E8BC "_destroy_dirac_buf__Fv$LAJ"
           0   8   809800    56D6D "_dirac"
        92FE   3    502BB    502BB "_DiracArg__11DiracOpStagFP5CgArg"
        92FD   3    502BA    502BA "_DiracArg__11DiracOpStagFP5CgArg$LAJ"
        A0C7   3    51084    51084 "_DiracArg__13DiracOpCloverFP5CgArg"
        A0C6   3    51083    51083 "_DiracArg__13DiracOpCloverFP5CgArg$LAJ"
        CD63   3    53D20    53D20 "_DiracOpGlbSum__7DiracOpFP6rfloat"
        CD62   3    53D1F    53D1F "_DiracOpGlbSum__7DiracOpFP6rfloat$LAJ"
        7927   3    4E8E4    4E8E4 "_dirac_init__FPCv"
        7926   3    4E8E3    4E8E3 "_dirac_init__FPCv$LAJ"
        F301   3    562BE    562BE "_divf32f64"
        F311   3    562CE    562CE "_divf64"
        F309   3    562C6    562C6 "_divf64f32"
           D   4   809FA2    5680C "_dotProduct"
        D93F   3    548FC    548FC "_double64round"
        2324   1     3324     3324 "_do_arg"
        9380   3    5033D    5033D
			    "_Dslash__11DiracOpStagFP6VectorT18ChkbType7DagType"
        937F   3    5033C    5033C
			"_Dslash__11DiracOpStagFP6VectorT18ChkbType7DagType$LAJ"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  39
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        A279   3    51236    51236
			  "_Dslash__13DiracOpCloverFP6VectorT18ChkbType7DagType"
        A278   3    51235    51235
		      "_Dslash__13DiracOpCloverFP6VectorT18ChkbType7DagType$LAJ"
           0   1     1000     1000 "_ERR"
        321C   1     421C     421C "_errno"
        1D45   1     2D45     2D45 "_error_class_name"
        1D4C   1     2D4C     2D4C "_error_file_name"
        3FF2   3    4AFAF    4AFAF "_EvolveGfield__7LatticeFP6Matrix6rfloat"
        3FF1   3    4AFAE    4AFAE "_EvolveGfield__7LatticeFP6Matrix6rfloat$LAJ"
        5496   3    4C453    4C453
		   "_EvolveMomFforceSupp__7FcloverFP6MatrixP6VectorN326rfloatT6"
        5495   3    4C452    4C452
	       "_EvolveMomFforceSupp__7FcloverFP6MatrixP6VectorN326rfloatT6$LAJ"
        5230   3    4C1ED    4C1ED
			    "_EvolveMomFforce__5FstagFP6MatrixP6Vector6rfloatT3"
        522F   3    4C1EC    4C1EC
			"_EvolveMomFforce__5FstagFP6MatrixP6Vector6rfloatT3$LAJ"
        6391   3    4D34E    4D34E
			  "_EvolveMomFforce__7FcloverFP6MatrixP6Vector6rfloatT3"
        6390   3    4D34D    4D34D
		      "_EvolveMomFforce__7FcloverFP6MatrixP6Vector6rfloatT3$LAJ"
        49CA   3    4B987    4B987 "_EvolveMomGforce__7GwilsonFP6Matrix6rfloat"
        49C9   3    4B986    4B986
				"_EvolveMomGforce__7GwilsonFP6Matrix6rfloat$LAJ"
        4FB1   3    4BF6E    4BF6E "_ExactFlavors__5FstagFv"
        4FB0   3    4BF6D    4BF6D "_ExactFlavors__5FstagFv$LAJ"
        60B8   3    4D075    4D075 "_ExactFlavors__7FcloverFv"
        60B7   3    4D074    4D074 "_ExactFlavors__7FcloverFv$LAJ"
        EB22   3    55ADF    55ADF "_exit"
        EB21   3    55ADE    55ADE "_exit$LAJ"
        EECE   3    55E8B    55E8B "_expf64"
        E81B   3    557D8    557D8 "_exp__FRC8double64"
        E81A   3    557D7    557D7 "_exp__FRC8double64$LAJ"
        2266   1     3266     3266 "_e_tab"
        320E   1     420E     420E "_f64Cube_Root_2"
        3212   1     4212     4212 "_f64Cube_Root_3"
        31FE   1     41FE     41FE "_f64e"
        31DA   1     41DA     41DA "_f64Eps"
        31DE   1     41DE     41DE "_f64Epsilon"
        31DC   1     41DC     41DC "_f64EpsNeg"
        31E2   1     41E2     41E2 "_f64Exp_Large"
        31E6   1     41E6     41E6 "_f64Exp_Small"
        3206   1     4206     4206 "_f64e_To_Half_Pi"
        3208   1     4208     4208 "_f64e_To_Minus_Half_Pi"
        3204   1     4204     4204 "_f64e_To_Minus_Pi"
        3202   1     4202     4202 "_f64e_To_Pi"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  40
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        31D4   1     41D4     41D4 "_f64First"
        31EC   1     41EC     41EC "_f64Integer_First"
        31EE   1     41EE     41EE "_f64Integer_Last"
        31D6   1     41D6     41D6 "_f64Last"
        320A   1     420A     420A "_f64Ln_2"
        3214   1     4214     4214 "_f64Ln_3"
        31E0   1     41E0     41E0 "_f64Ln_Xmax"
        31E4   1     41E4     41E4 "_f64Ln_Xmin"
        3210   1     4210     4210 "_f64Log_2"
        3216   1     4216     4216 "_f64Log_3"
        3200   1     4200     4200 "_f64One_Over_e"
        31F4   1     41F4     41F4 "_f64One_Over_Pi"
        31F0   1     41F0     41F0 "_f64Pi"
        31FC   1     41FC     41FC "_f64Pi_Over_Four"
        31FA   1     41FA     41FA "_f64Pi_Over_Three"
        31F8   1     41F8     41F8 "_f64Pi_Over_Two"
        2FA0   1     3FA0     3FA0 "_F64PowersOfTen"
        31D2   1     41D2     41D2 "_f64SmallNeg"
        31D0   1     41D0     41D0 "_f64SmallPos"
        320C   1     420C     420C "_f64Sqrt_2"
        31F6   1     41F6     41F6 "_f64Two_Over_Pi"
        31F2   1     41F2     41F2 "_f64Two_Pi"
        31E8   1     41E8     41E8 "_f64Xbig"
        31D8   1     41D8     41D8 "_f64Xmin"
        31EA   1     41EA     41EA "_f64Ymax"
        7101   3    4E0BE    4E0BE "_FcanonToWilson__FP16ConvertArgStruct"
        7100   3    4E0BD    4E0BD "_FcanonToWilson__FP16ConvertArgStruct$LAJ"
        4FD4   3    4BF91    4BF91 "_FchkbEvl__5FstagFv"
        4FD3   3    4BF90    4BF90 "_FchkbEvl__5FstagFv$LAJ"
        60E3   3    4D0A0    4D0A0 "_FchkbEvl__7FcloverFv"
        60E2   3    4D09F    4D09F "_FchkbEvl__7FcloverFv$LAJ"
        4F7B   3    4BF38    4BF38 "_Fclass__5FstagFv"
        4F7A   3    4BF37    4BF37 "_Fclass__5FstagFv$LAJ"
        6025   3    4CFE2    4CFE2 "_Fclass__7FcloverFv"
        6024   3    4CFE1    4CFE1 "_Fclass__7FcloverFv$LAJ"
        DE09   3    54DC6    54DC6 "_fclose"
        772B   3    4E6E8    4E6E8 "_Fconvert__5FstagFP6Vector10StrOrdTypeT2"
        772A   3    4E6E7    4E6E7
				  "_Fconvert__5FstagFP6Vector10StrOrdTypeT2$LAJ"
        7749   3    4E706    4E706 "_Fconvert__7FcloverFP6Vector10StrOrdTypeT2"
        7748   3    4E705    4E705
				"_Fconvert__7FcloverFP6Vector10StrOrdTypeT2$LAJ"
        5091   3    4C04E    4C04E
	  "_FeigSolv__5FstagFPP6VectorP6rfloatT2PiPP6rfloatP6EigArg10CnvFrmType"
        5090   3    4C04D    4C04D
      "_FeigSolv__5FstagFPP6VectorP6rfloatT2PiPP6rfloatP6EigArg10CnvFrmType$LAJ"
        618E   3    4D14B    4D14B

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  41
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

	"_FeigSolv__7FcloverFPP6VectorP6rfloatT2PiPP6rfloatP6EigArg10CnvFrmType"
        618D   3    4D14A    4D14A
    "_FeigSolv__7FcloverFPP6VectorP6rfloatT2PiPP6rfloatP6EigArg10CnvFrmType$LAJ"
        4407   3    4B3C4    4B3C4 "_Ffive2four__7LatticeFP6VectorT1iT3"
        4406   3    4B3C3    4B3C3 "_Ffive2four__7LatticeFP6VectorT1iT3$LAJ"
        DE0B   3    54DC8    54DC8 "_fflush"
        5196   3    4C153    4C153 "_FforceSite__5FstagFR6MatrixP6VectorPii"
        5195   3    4C152    4C152 "_FforceSite__5FstagFR6MatrixP6VectorPii$LAJ"
        43F1   3    4B3AE    4B3AE "_Ffour2five__7LatticeFP6VectorT1iT3"
        43F0   3    4B3AD    4B3AD "_Ffour2five__7LatticeFP6VectorT1iT3$LAJ"
        52F0   3    4C2AD    4C2AD "_FhamiltonNode__5FstagFP6VectorT1"
        52EF   3    4C2AC    4C2AC "_FhamiltonNode__5FstagFP6VectorT1$LAJ"
        6946   3    4D903    4D903 "_FhamiltonNode__7FcloverFP6VectorT1"
        6945   3    4D902    4D902 "_FhamiltonNode__7FcloverFP6VectorT1$LAJ"
        6BB5   3    4DB72    4DB72 "_FileA__5ErrorFPcN21"
        6BB4   3    4DB71    4DB71 "_FileA__5ErrorFPcN21$LAJ"
        1D48   1     2D48     2D48 "_file_a_str"
        1D46   1     2D46     2D46 "_file_r_str"
        1D47   1     2D47     2D47 "_file_w_str"
        ED0C   3    55CC9    55CC9 "_fixf32"
         DF7   1     1DF7     1DF7 "_fix_gauge_kind__7Lattice"
         DFD   1     1DFD     1DFD "_fix_gauge_ptr__7Lattice"
        ED55   3    55D12    55D12 "_floorf40f64"
        ED46   3    55D03    55D03 "_floori32f64"
        1D76   3    48D33    48D33 "_Flow__7VerboseFPcT1PCce"
        1D75   3    48D32    48D32 "_Flow__7VerboseFPcT1PCce$LAJ"
        501C   3    4BFD9    4BFD9
			     "_FmatEvlInv__5FstagFP6VectorT1P5CgArg10CnvFrmType"
        501B   3    4BFD8    4BFD8
			 "_FmatEvlInv__5FstagFP6VectorT1P5CgArg10CnvFrmType$LAJ"
        4FD7   3    4BF94    4BF94
		     "_FmatEvlInv__5FstagFP6VectorT1P5CgArgP6rfloat10CnvFrmType"
        4FD6   3    4BF93    4BF93
		 "_FmatEvlInv__5FstagFP6VectorT1P5CgArgP6rfloat10CnvFrmType$LAJ"
        611E   3    4D0DB    4D0DB
			   "_FmatEvlInv__7FcloverFP6VectorT1P5CgArg10CnvFrmType"
        611D   3    4D0DA    4D0DA
		       "_FmatEvlInv__7FcloverFP6VectorT1P5CgArg10CnvFrmType$LAJ"
        60EA   3    4D0A7    4D0A7
		   "_FmatEvlInv__7FcloverFP6VectorT1P5CgArgP6rfloat10CnvFrmType"
        60E9   3    4D0A6    4D0A6
	       "_FmatEvlInv__7FcloverFP6VectorT1P5CgArgP6rfloat10CnvFrmType$LAJ"
        5073   3    4C030    4C030
		  "_FmatInv__5FstagFP6VectorT1P5CgArg10CnvFrmType12PreserveType"
        5072   3    4C02F    4C02F
	      "_FmatInv__5FstagFP6VectorT1P5CgArg10CnvFrmType12PreserveType$LAJ"
        5038   3    4BFF5    4BFF5

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  42
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

	  "_FmatInv__5FstagFP6VectorT1P5CgArgP6rfloat10CnvFrmType12PreserveType"
        5037   3    4BFF4    4BFF4
      "_FmatInv__5FstagFP6VectorT1P5CgArgP6rfloat10CnvFrmType12PreserveType$LAJ"
        6170   3    4D12D    4D12D
		"_FmatInv__7FcloverFP6VectorT1P5CgArg10CnvFrmType12PreserveType"
        616F   3    4D12C    4D12C
	    "_FmatInv__7FcloverFP6VectorT1P5CgArg10CnvFrmType12PreserveType$LAJ"
        613A   3    4D0F7    4D0F7
	"_FmatInv__7FcloverFP6VectorT1P5CgArgP6rfloat10CnvFrmType12PreserveType"
        6139   3    4D0F6    4D0F6
    "_FmatInv__7FcloverFP6VectorT1P5CgArgP6rfloat10CnvFrmType12PreserveType$LAJ"
        1FE6   1     2FE6     2FE6 "_fname_fconvert"
        DE07   3    54DC4    54DC4 "_fopen"
        DE08   3    54DC5    54DC5 "_fprintf"
        441D   3    4B3DA    4B3DA
				  "_FrandGaussVector__7LatticeFP6Vector6rfloati"
        441C   3    4B3D9    4B3D9
			      "_FrandGaussVector__7LatticeFP6Vector6rfloati$LAJ"
        E799   3    55756    55756 "_free"
        E798   3    55755    55755 "_free$LAJ"
        602C   3    4CFE9    4CFE9 "_FsiteOffsetChkb__7FcloverCFPCi"
        602B   3    4CFE8    4CFE8 "_FsiteOffsetChkb__7FcloverCFPCi$LAJ"
        4F9C   3    4BF59    4BF59 "_FsiteOffset__5FstagCFPCi"
        4F9B   3    4BF58    4BF58 "_FsiteOffset__5FstagCFPCi$LAJ"
        6086   3    4D043    4D043 "_FsiteOffset__7FcloverCFPCi"
        6085   3    4D042    4D042 "_FsiteOffset__7FcloverCFPCi$LAJ"
        4FB7   3    4BF74    4BF74 "_FsiteSize__5FstagFv"
        4FB6   3    4BF73    4BF73 "_FsiteSize__5FstagFv$LAJ"
        60C6   3    4D083    4D083 "_FsiteSize__7FcloverFv"
        60C5   3    4D082    4D082 "_FsiteSize__7FcloverFv$LAJ"
        1FFE   3    48FBB    48FBB "_fTimesV1MinusV2"
          1B   4   809FB0    5681A "_fTimesV1PlusV2"
        1CF2   3    48CAF    48CAF "_FuncEnd__7VerboseFPcT1"
        1CF1   3    48CAE    48CAE "_FuncEnd__7VerboseFPcT1$LAJ"
        1CC3   3    48C80    48C80 "_Func__7VerboseFPcT1"
        1CC2   3    48C7F    48C7F "_Func__7VerboseFPcT1$LAJ"
        7250   3    4E20D    4E20D "_FwilsonToCanon__FP16ConvertArgStruct"
        724F   3    4E20C    4E20C "_FwilsonToCanon__FP16ConvertArgStruct$LAJ"
        48E1   3    4B89E    4B89E "_GactionGradient__7GwilsonFR6MatrixPii"
        48E0   3    4B89D    4B89D "_GactionGradient__7GwilsonFR6MatrixPii$LAJ"
        4BB9   3    4BB76    4BB76 "_Gamma5__12FwilsonTypesFP6VectorT1i"
        4BB8   3    4BB75    4BB75 "_Gamma5__12FwilsonTypesFP6VectorT1i$LAJ"
        43DB   3    4B398    4B398 "_Gamma5__7LatticeFP6VectorT1i"
        43DA   3    4B397    4B397 "_Gamma5__7LatticeFP6VectorT1i$LAJ"
        483F   3    4B7FC    4B7FC "_gamma_5__FPfT1i"
        4840   3    4B7FD    4B7FD "_gamma_5__FPfT1i$LAJ"
        3A06   3    4A9C3    4A9C3 "_GaugeField__7LatticeCFv"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  43
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        3A05   3    4A9C2    4A9C2 "_GaugeField__7LatticeCFv$LAJ"
        3A0A   3    4A9C7    4A9C7 "_GaugeField__7LatticeFP6Matrix"
        3A09   3    4A9C6    4A9C6 "_GaugeField__7LatticeFP6Matrix$LAJ"
         DF8   1     1DF8     1DF8 "_gauge_field__7Lattice"
        48DE   3    4B89B    4B89B "_Gclass__7GwilsonFv"
        48DD   3    4B89A    4B89A "_Gclass__7GwilsonFv$LAJ"
        1D4B   1     2D4B     2D4B "_general_str"
        6CA2   3    4DC5F    4DC5F "_General__5ErrorFPcT1PCce"
        6CA1   3    4DC5E    4DC5E "_General__5ErrorFPcT1PCce$LAJ"
        373D   3    4A6FA    4A6FA "_GetLinkOld__7LatticeCFP6MatrixPCiiT3"
        373C   3    4A6F9    4A6F9 "_GetLinkOld__7LatticeCFP6MatrixPCiiT3$LAJ"
        9ABB   3    50A78    50A78 "_GetLink__13DiracOpCloverCFPCii"
        9ABA   3    50A77    50A77 "_GetLink__13DiracOpCloverCFPCii$LAJ"
        DB90   3    54B4D    54B4D "_getMinusData__FPfT1iT3"
        DB8F   3    54B4C    54B4C "_getMinusData__FPfT1iT3$LAJ"
        DB2A   3    54AE7    54AE7 "_getPlusData__FPfT1iT3"
        DB29   3    54AE6    54AE6 "_getPlusData__FPfT1iT3$LAJ"
        4CE9   3    4BCA6    4BCA6 "_getUDagX__5FstagCFR6VectorPC6VectorPii"
        4CE8   3    4BCA5    4BCA5 "_getUDagX__5FstagCFR6VectorPC6VectorPii$LAJ"
        4927   3    4B8E4    4B8E4 "_GforceSite__7GwilsonFR6MatrixPii"
        4926   3    4B8E3    4B8E3 "_GforceSite__7GwilsonFR6MatrixPii$LAJ"
        4A8F   3    4BA4C    4BA4C "_GhamiltonNode__7GwilsonFv"
        4A8E   3    4BA4B    4BA4B "_GhamiltonNode__7GwilsonFv$LAJ"
      809800   9   809800   8607EF "_ghb_dest"
         1FD   1     11FD     11FD "_GJP"
        13C4   1     23C4     23C4 "_gjp_local_axis"
        13BA   1     23BA     23BA "_gjp_scu_dir"
        D949   3    54906    54906 "_glb_sum_five__FP6rfloat"
        D948   3    54905    54905 "_glb_sum_five__FP6rfloat$LAJ"
        DA7D   3    54A3A    54A3A "_glb_sum__FP6rfloat"
        DA7C   3    54A39    54A39 "_glb_sum__FP6rfloat$LAJ"
        B23D   3    521FA    521FA "_GramSchm__FPP6VectoriT1N22"
        B23C   3    521F9    521F9 "_GramSchm__FPP6VectoriT1N22$LAJ"
        3A60   3    4AA1D    4AA1D "_GsiteSize__7LatticeFv"
        3A5F   3    4AA1C    4AA1C "_GsiteSize__7LatticeFv$LAJ"
        4224   3    4B1E1    4B1E1 "_GsoCheck__7LatticeFv"
        4223   3    4B1E0    4B1E0 "_GsoCheck__7LatticeFv$LAJ"
        41CD   3    4B18A    4B18A "_GupdCntInc__7LatticeFi"
        41CC   3    4B189    4B189 "_GupdCntInc__7LatticeFi$LAJ"
        41C7   3    4B184    4B184 "_GupdCnt__7LatticeFv"
        41C6   3    4B183    4B183 "_GupdCnt__7LatticeFv$LAJ"
         E5B   1     1E5B     1E5B "_g_dir_offset__7Lattice"
         DFC   1     1DFC     1DFC "_g_upd_cnt__7Lattice"
        1D4A   1     2D4A     2D4A "_hardware_str"
         BE5   1     1BE5     1BE5 "_inextp__15RandomGenerator"
         BE6   1     1BE6     1BE6 "_inext__15RandomGenerator"
        44F0   3    4B4AD    4B4AD "_Initialize__18GlobalJobParameterFRC5DoArg"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  44
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        44EF   3    4B4AC    4B4AC
				"_Initialize__18GlobalJobParameterFRC5DoArg$LAJ"
        DDA3   3    54D60    54D60 "_Init__9SCUDirArgFPv6SCUDir5SCUXRiN24"
        DDA2   3    54D5F    54D5F "_Init__9SCUDirArgFPv6SCUDir5SCUXRiN24$LAJ"
        1DB4   3    48D71    48D71 "_Input__7VerboseFPcT1PCce"
        1DB3   3    48D70    48D70 "_Input__7VerboseFPcT1PCce$LAJ"
        CC97   3    53C54    53C54 "_InvCg__7DiracOpFP6rfloat"
        CC96   3    53C53    53C53 "_InvCg__7DiracOpFP6rfloat$LAJ"
        C7EA   3    537A7    537A7 "_InvCg__7DiracOpFP6VectorT16rfloatP6rfloat"
        C7E9   3    537A6    537A6
				"_InvCg__7DiracOpFP6VectorT16rfloatP6rfloat$LAJ"
        CC74   3    53C31    53C31 "_InvCg__7DiracOpFP6VectorT1P6rfloat"
        CC73   3    53C30    53C30 "_InvCg__7DiracOpFP6VectorT1P6rfloat$LAJ"
         DF9   1     1DF9     1DF9 "_is_allocated__7Lattice"
         DFA   1     1DFA     1DFA "_is_initialized__7Lattice"
        BEA0   3    52E5D    52E5D
		       "_Jacobi__7DiracOpFPP6VectoriP6rfloatP8Rcomplex6rfloatT2"
        BE9F   3    52E5C    52E5C
		   "_Jacobi__7DiracOpFPP6VectoriP6rfloatP8Rcomplex6rfloatT2$LAJ"
        1ED1   3    48E8E    48E8E "_LedFlash__7VerboseFPcT1i"
        1ED0   3    48E8D    48E8D "_LedFlash__7VerboseFPcT1i$LAJ"
        1EBF   3    48E7C    48E7C "_LedOff__7VerboseFPcT1"
        1EBE   3    48E7B    48E7B "_LedOff__7VerboseFPcT1$LAJ"
        1EAD   3    48E6A    48E6A "_LedOn__7VerboseFPcT1"
        1EAC   3    48E69    48E69 "_LedOn__7VerboseFPcT1$LAJ"
        1C09   3    48BC6    48BC6 "_Level__7VerboseFi"
        1C08   3    48BC5    48BC5 "_Level__7VerboseFi$LAJ"
        DDF3   3    54DB0    54DB0 "_Ln2NumNodes"
        F50D   3    564CA    564CA "_log"
        F559   3    56516    56516 "_log10"
        F52D   3    564EA    564EA "_log2"
        ECFF   3    55CBC    55CBC "_lssf64"
           1   3    46FBE    46FBE "_main"
           0   3    46FBD    46FBD "_main$LAJ"
        E751   3    5570E    5570E "_malloc"
        E750   3    5570D    5570D "_malloc$LAJ"
        D7BA   3    54777    54777 "_MatDagMat__18DiracOpWilsonTypesFP6VectorT1"
        D7B9   3    54776    54776
			       "_MatDagMat__18DiracOpWilsonTypesFP6VectorT1$LAJ"
        A461   3    5141E    5141E "_MatDag__13DiracOpCloverFP6VectorT1"
        A460   3    5141D    5141D "_MatDag__13DiracOpCloverFP6VectorT1$LAJ"
        A3D2   3    5138F    5138F "_MatEvlInv__13DiracOpCloverFP6rfloat"
        A3D1   3    5138E    5138E "_MatEvlInv__13DiracOpCloverFP6rfloat$LAJ"
        A375   3    51332    51332
				"_MatEvlInv__13DiracOpCloverFP6VectorT1P6rfloat"
        A374   3    51331    51331
			    "_MatEvlInv__13DiracOpCloverFP6VectorT1P6rfloat$LAJ"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  45
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        A3DF   3    5139C    5139C "_MatHerm__13DiracOpCloverFP6VectorT1"
        A3DE   3    5139B    5139B "_MatHerm__13DiracOpCloverFP6VectorT1$LAJ"
        94A3   3    50460    50460 "_MatInv__11DiracOpStagF12PreserveType"
        94A2   3    5045F    5045F "_MatInv__11DiracOpStagF12PreserveType$LAJ"
        9489   3    50446    50446
				 "_MatInv__11DiracOpStagFP6rfloat12PreserveType"
        9488   3    50445    50445
			     "_MatInv__11DiracOpStagFP6rfloat12PreserveType$LAJ"
        946F   3    5042C    5042C
			       "_MatInv__11DiracOpStagFP6VectorT112PreserveType"
        946E   3    5042B    5042B
			   "_MatInv__11DiracOpStagFP6VectorT112PreserveType$LAJ"
        93AD   3    5036A    5036A
		       "_MatInv__11DiracOpStagFP6VectorT1P6rfloat12PreserveType"
        93AC   3    50369    50369
		   "_MatInv__11DiracOpStagFP6VectorT1P6rfloat12PreserveType$LAJ"
        A35B   3    51318    51318 "_MatInv__13DiracOpCloverF12PreserveType"
        A35A   3    51317    51317 "_MatInv__13DiracOpCloverF12PreserveType$LAJ"
        A341   3    512FE    512FE
			       "_MatInv__13DiracOpCloverFP6rfloat12PreserveType"
        A340   3    512FD    512FD
			   "_MatInv__13DiracOpCloverFP6rfloat12PreserveType$LAJ"
        A327   3    512E4    512E4
			     "_MatInv__13DiracOpCloverFP6VectorT112PreserveType"
        A326   3    512E3    512E3
			 "_MatInv__13DiracOpCloverFP6VectorT112PreserveType$LAJ"
        A28D   3    5124A    5124A
		     "_MatInv__13DiracOpCloverFP6VectorT1P6rfloat12PreserveType"
        A28C   3    51249    51249
		 "_MatInv__13DiracOpCloverFP6VectorT1P6rfloat12PreserveType$LAJ"
        9329   3    502E6    502E6
			      "_MatPcDagMatPc__11DiracOpStagFP6VectorT1P6rfloat"
        9328   3    502E5    502E5
			  "_MatPcDagMatPc__11DiracOpStagFP6VectorT1P6rfloat$LAJ"
        A21F   3    511DC    511DC
			    "_MatPcDagMatPc__13DiracOpCloverFP6VectorT1P6rfloat"
        A21E   3    511DB    511DB
			"_MatPcDagMatPc__13DiracOpCloverFP6VectorT1P6rfloat$LAJ"
        A03B   3    50FF8    50FF8
			   "_MatPcDagOrNot__13DiracOpCloverCFP6VectorPC6Vectori"
        A03A   3    50FF7    50FF7
		       "_MatPcDagOrNot__13DiracOpCloverCFP6VectorPC6Vectori$LAJ"
        A5C6   3    51583    51583 "_mat_hrm_cmpr"
        A5C7   3    51584    51584 "_mat_hrm_cmpr$LAJ"
        A5E5   3    515A2    515A2 "_mat_hrm_decm"
        A5E4   3    515A1    515A1 "_mat_hrm_decm$LAJ"
        A63E   3    515FB    515FB "_mat_hrm_ldl"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  46
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        A63D   3    515FA    515FA "_mat_hrm_ldl$LAJ"
        A91C   3    518D9    518D9 "_mat_inv"
        A91B   3    518D8    518D8 "_mat_inv$LAJ"
        A44D   3    5140A    5140A "_Mat__13DiracOpCloverFP6VectorT1"
        A44C   3    51409    51409 "_Mat__13DiracOpCloverFP6VectorT1$LAJ"
         BE7   1     1BE7     1BE7 "_ma__15RandomGenerator"
        DDE7   3    54DA4    54DA4 "_MbNum"
        1F72   3    48F2F    48F2F "_mDotMEqual"
        1F9C   3    48F59    48F59 "_mDotMPlus"
        4201   3    4B1BE    4B1BE "_MdTimeInc__7LatticeF6rfloat"
        4200   3    4B1BD    4B1BD "_MdTimeInc__7LatticeF6rfloat$LAJ"
        41EC   3    4B1A9    4B1A9 "_MdTime__7LatticeF6rfloat"
        41EB   3    4B1A8    4B1A8 "_MdTime__7LatticeF6rfloat$LAJ"
        41DD   3    4B19A    4B19A "_MdTime__7LatticeFv"
        41DC   3    4B199    4B199 "_MdTime__7LatticeFv$LAJ"
         DFE   1     1DFE     1DFE "_md_time__7Lattice"
        E5CA   3    55587    55587 "_memmove"
        3F14   3    4AED1    4AED1 "_MetropolisAccept__7LatticeF6rfloat"
        3F13   3    4AED0    4AED0 "_MetropolisAccept__7LatticeF6rfloat$LAJ"
        E73F   3    556FC    556FC "_minit"
        E740   3    556FD    556FD "_minit$LAJ"
        3774   3    4A731    4A731 "_MltFloatImpl__7LatticeF6rfloati"
        3773   3    4A730    4A730 "_MltFloatImpl__7LatticeF6rfloati$LAJ"
        4090   3    4B04D    4B04D "_MomHamiltonNode__7LatticeFP6Matrix"
        408F   3    4B04C    4B04C "_MomHamiltonNode__7LatticeFP6Matrix$LAJ"
           0   4   809F95    567FF "_moveMem"
        EC5B   3    55C18    55C18 "_mpyf32f40"
        ECAB   3    55C68    55C68 "_mpyf40"
        EC6F   3    55C2C    55C2C "_mpyf40f32"
        EBF6   3    55BB3    55BB3 "_mpyf64"
        EBE3   3    55BA0    55BA0 "_mpyf64f32"
        D717   3    546D4    546D4
			 "_MultGamma__18DiracOpWilsonTypesFP6VectorPC6VectoriT3"
        D716   3    546D3    546D3
		     "_MultGamma__18DiracOpWilsonTypesFP6VectorPC6VectoriT3$LAJ"
        6D01   3    4DCBE    4DCBE "_MultStagPhases__FP16ConvertArgStruct"
        6D00   3    4DCBD    4DCBD "_MultStagPhases__FP16ConvertArgStruct$LAJ"
        2263   1     3263     3263 "_m_rec_buf"
        77EA   3    4E7A7    4E7A7 "_negate_link"
        77EB   3    4E7A8    4E7A8 "_negate_link$LAJ"
        F49C   3    56459    56459 "_negf64"
        26E7   3    496A4    496A4 "_NegHalfTrSquare__6MatrixCFv"
        26E8   3    496A5    496A5 "_NegHalfTrSquare__6MatrixCFv$LAJ"
        EBD9   3    55B96    55B96 "_neqf64"
         E56   1     1E56     1E56 "_node_sites__7Lattice"
        2260   3    4921D    4921D "_normalize__FP6rfloat"
        225F   3    4921C    4921C "_normalize__FP6rfloat$LAJ"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  47
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        21BB   3    49178    49178 "_NormSqGlbSum__6VectorFi"
        21BA   3    49177    49177 "_NormSqGlbSum__6VectorFi$LAJ"
        3371   3    4A32E    4A32E "_norm__8RcomplexCFv"
        3372   3    4A32F    4A32F "_norm__8RcomplexCFv$LAJ"
        6BFE   3    4DBBB    4DBBB "_NotImplemented__5ErrorFPcT1"
        6BFD   3    4DBBA    4DBBA "_NotImplemented__5ErrorFPcT1$LAJ"
        6C43   3    4DC00    4DC00 "_NotImplemented__5ErrorFPcT1PCce"
        6C42   3    4DBFF    4DBFF "_NotImplemented__5ErrorFPcT1PCce$LAJ"
        1D49   1     2D49     2D49 "_not_implemented_str"
        231C   1     331C     331C "_nsurf_sites"
        DDF2   3    54DAF    54DAF "_NumNodes"
        2264   1     3264     3264 "_n_even_sites"
        2265   1     3265     3265 "_n_odd_sites"
        200F   3    48FCC    48FCC "_oneMinusfTimesMatrix"
        22AB   3    49268    49268 "_orthogonalize__FP6rfloatPC6rfloat"
        22AA   3    49267    49267 "_orthogonalize__FP6rfloatPC6rfloat$LAJ"
        22C1   1     32C1     32C1 "_o_tab"
        DCD6   3    54C93    54C93 "_p2vCloverLib__Fv"
        DCD5   3    54C92    54C92 "_p2vCloverLib__Fv$LAJ"
        DC99   3    54C56    54C56 "_p2vStagDs__Fv"
        DC98   3    54C55    54C55 "_p2vStagDs__Fv$LAJ"
        DBF6   3    54BB3    54BB3 "_p2vVector__Fv"
        DBF5   3    54BB2    54BB2 "_p2vVector__Fv$LAJ"
        DC33   3    54BF0    54BF0 "_p2vWilsonLib__Fv"
        DC32   3    54BEF    54BEF "_p2vWilsonLib__Fv$LAJ"
        1D21   3    48CDE    48CDE "_Pmalloc__7VerboseFPcN21Pvi"
        1D20   3    48CDD    48CDD "_Pmalloc__7VerboseFPcN21Pvi$LAJ"
        351F   3    4A4DC    4A4DC "_pmalloc__Fi"
        351E   3    4A4DB    4A4DB "_pmalloc__Fi$LAJ"
        1D4D   1     2D4D     2D4D "_pointer_str"
        6B6C   3    4DB29    4DB29 "_Pointer__5ErrorFPcN21"
        6B6B   3    4DB28    4DB28 "_Pointer__5ErrorFPcN21$LAJ"
        F021   3    55FDE    55FDE "_powf64"
        DDF8   3    54DB5    54DB5 "_printf"
        DE15   3    54DD2    54DD2 "_putchar"
        2261   1     3261     3261 "_p_rec_buf"
        2262   1     3262     3262 "_p_send_buf"
        40E0   3    4B09D    4B09D
			    "_RandGaussAntiHermMatrix__7LatticeFP6Matrix6rfloat"
        40DF   3    4B09C    4B09C
			"_RandGaussAntiHermMatrix__7LatticeFP6Matrix6rfloat$LAJ"
        2213   3    491D0    491D0 "_RandGaussVector__6VectorF6rfloati"
        2212   3    491CF    491CF "_RandGaussVector__6VectorF6rfloati$LAJ"
        213F   3    490FC    490FC "_RandLink__6MatrixFv"
        213E   3    490FB    490FB "_RandLink__6MatrixFv$LAJ"
        3481   3    4A43E    4A43E "_Rand__15RandomGeneratorFv"
        3480   3    4A43D    4A43D "_Rand__15RandomGeneratorFv$LAJ"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  48
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        3472   3    4A42F    4A42F "_Rand__22UniformRandomGeneratorFv"
        3471   3    4A42E    4A42E "_Rand__22UniformRandomGeneratorFv$LAJ"
        3430   3    4A3ED    4A3ED "_Rand__23GaussianRandomGeneratorFv"
        342F   3    4A3EC    4A3EC "_Rand__23GaussianRandomGeneratorFv$LAJ"
        21D4   3    49191    49191 "_ReDotProductGlbSum__6VectorFPC6Vectori"
        21D3   3    49190    49190 "_ReDotProductGlbSum__6VectorFPC6Vectori$LAJ"
        ED28   3    55CE5    55CE5 "_refloatf32"
        ED34   3    55CF1    55CF1 "_refloatf40"
        DDD3   3    54D90    54D90 "_Reload__9SCUDirArgFPviN22"
        DDD4   3    54D91    54D91 "_Reload__9SCUDirArgFPviN22$LAJ"
        34BB   3    4A478    4A478 "_Reset__15RandomGeneratorSFi"
        34BA   3    4A477    4A477 "_Reset__15RandomGeneratorSFi$LAJ"
        DD4D   3    54D0A    54D0A "_restoreCbufCntrlReg__Fv"
        DD4E   3    54D0B    54D0B "_restoreCbufCntrlReg__Fv$LAJ"
        1DDD   3    48D9A    48D9A "_Result__7VerboseFPcT1PCce"
        1DDC   3    48D99    48D99 "_Result__7VerboseFPcT1PCce$LAJ"
        3C01   3    4ABBE    4ABBE "_ReTrPlaq__7LatticeCFPiiT2"
        3C00   3    4ABBD    4ABBD "_ReTrPlaq__7LatticeCFPiiT2$LAJ"
        270A   3    496C7    496C7 "_ReTr__6MatrixCFv"
        270B   3    496C8    496C8 "_ReTr__6MatrixCFv$LAJ"
        3DC3   3    4AD80    4AD80 "_Reunitarize__7LatticeFR6rfloatT1"
        3DC2   3    4AD7F    4AD7F "_Reunitarize__7LatticeFR6rfloatT1$LAJ"
        3D66   3    4AD23    4AD23 "_Reunitarize__7LatticeFv"
        3D65   3    4AD22    4AD22 "_Reunitarize__7LatticeFv$LAJ"
        94BD   3    5047A    5047A "_RitzEigMat__11DiracOpStagFP6VectorT1"
        94BC   3    50479    50479 "_RitzEigMat__11DiracOpStagFP6VectorT1$LAJ"
        D823   3    547E0    547E0
				  "_RitzEigMat__18DiracOpWilsonTypesFP6VectorT1"
        D822   3    547DF    547DF
			      "_RitzEigMat__18DiracOpWilsonTypesFP6VectorT1$LAJ"
        9546   3    50503    50503
		      "_RitzEig__16DiracOpStagTypesFPP6VectorP6rfloatPiP6EigArg"
        9545   3    50502    50502
		  "_RitzEig__16DiracOpStagTypesFPP6VectorP6rfloatPiP6EigArg$LAJ"
        CF4C   3    53F09    53F09
		    "_RitzEig__18DiracOpWilsonTypesFPP6VectorP6rfloatPiP6EigArg"
        CF4B   3    53F08    53F08
		"_RitzEig__18DiracOpWilsonTypesFPP6VectorP6rfloatPiP6EigArg$LAJ"
        99C3   3    50980    50980 "_RitzLatSize__16DiracOpStagTypesFv"
        99C2   3    5097F    5097F "_RitzLatSize__16DiracOpStagTypesFv$LAJ"
        D8F1   3    548AE    548AE "_RitzLatSize__18DiracOpWilsonTypesFv"
        D8F0   3    548AD    548AD "_RitzLatSize__18DiracOpWilsonTypesFv$LAJ"
        94D1   3    5048E    5048E "_RitzMat__11DiracOpStagFP6VectorT1"
        94D0   3    5048D    5048D "_RitzMat__11DiracOpStagFP6VectorT1$LAJ"
        D892   3    5484F    5484F "_RitzMat__18DiracOpWilsonTypesFP6VectorT1"
        D891   3    5484E    5484E
				 "_RitzMat__18DiracOpWilsonTypesFP6VectorT1$LAJ"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  49
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        B2AC   3    52269    52269
			  "_Ritz__7DiracOpFPP6VectoriR6rfloat6rfloatN34N42T4N22"
        B2AB   3    52268    52268
		      "_Ritz__7DiracOpFPP6VectoriR6rfloat6rfloatN34N42T4N22$LAJ"
        EBCE   3    55B8B    55B8B "_roundi32f64"
        ED62   3    55D1F    55D1F "_rsqrtf64"
        6DA7   3    4DD64    4DD64 "_RunGConverter__FP16ConvertArgStructPUiT2"
        6DA6   3    4DD63    4DD63
				 "_RunGConverter__FP16ConvertArgStructPUiT2$LAJ"
         64F   3    4760C    4760C "_run__7AlgHmdRFv"
         64E   3    4760B    4760B "_run__7AlgHmdRFv$LAJ"
        11F5   3    481B2    481B2 "_run__9AlgHmcPhiFv"
        11F4   3    481B1    481B1 "_run__9AlgHmcPhiFv$LAJ"
        21AD   1     31AD     31AD "_sarg"
        21A5   1     31A5     31A5 "_sargpB"
        219D   1     319D     319D "_sargpF"
        DD33   3    54CF0    54CF0 "_saveCbufCntrlReg__Fv"
        DD32   3    54CEF    54CEF "_saveCbufCntrlReg__Fv$LAJ"
        2A33   1     3A33     3A33 "_scope_lock__7DiracOp"
        2260   1     3260     3260 "_scratch_addr"
        DE03   3    54DC0    54DC0 "_SCUPoll__FP9SCUDirArg"
        DDFB   3    54DB8    54DB8 "_SCUReboot"
        DDFC   3    54DB9    54DB9 "_SCURemap__F6SCUDir"
        DDFA   3    54DB7    54DB7 "_SCUReset"
        DDFD   3    54DBA    54DBA "_SCUSetDMA__FP9SCUDirArg"
        DDFE   3    54DBB    54DBB "_SCUSetDMA__FPP9SCUDirArgi"
        DE01   3    54DBE    54DBE "_SCUTransAddr__FP9SCUDirArg"
        DE02   3    54DBF    54DBF "_SCUTransAddr__FPP9SCUDirArgi"
        DE05   3    54DC2    54DC2 "_SCUTransComplete__FP9SCUDirArg"
        DE06   3    54DC3    54DC3 "_SCUTransComplete__FPP9SCUDirArgi"
        DE04   3    54DC1    54DC1 "_SCUTransComplete__Fv"
        DDFF   3    54DBC    54DBC "_SCUTrans__FP9SCUDirArg"
        DE00   3    54DBD    54DBD "_SCUTrans__FPP9SCUDirArgi"
        DDF4   3    54DB1    54DB1 "_Seed"
        DDF5   3    54DB2    54DB2 "_SeedS"
        DDF7   3    54DB4    54DB4 "_SeedST"
        DDF6   3    54DB3    54DB3 "_SeedT"
        DD5D   3    54D1A    54D1A "_setCbufCntrlReg__FiUi"
        DD5E   3    54D1B    54D1B "_setCbufCntrlReg__FiUi$LAJ"
        4194   3    4B151    4B151 "_SetGfieldDisOrd__7LatticeFv"
        4193   3    4B150    4B150 "_SetGfieldDisOrd__7LatticeFv$LAJ"
        4161   3    4B11E    4B11E "_SetGfieldOrd__7LatticeFv"
        4160   3    4B11D    4B11D "_SetGfieldOrd__7LatticeFv$LAJ"
        5114   3    4C0D1    4C0D1 "_SetPhi__5FstagFP6VectorN216rfloat"
        5113   3    4C0D0    4C0D0 "_SetPhi__5FstagFP6VectorN216rfloat$LAJ"
        62DB   3    4D298    4D298 "_SetPhi__7FcloverFP6VectorN216rfloat"
        62DA   3    4D297    4D297 "_SetPhi__7FcloverFP6VectorN216rfloat$LAJ"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  50
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        1D5B   3    48D18    48D18 "_Sfree__7VerboseFPcN21Pv"
        1D5A   3    48D17    48D17 "_Sfree__7VerboseFPcN21Pv$LAJ"
        2C4E   3    49C0B    49C0B "_sfree__FPv"
        2C4D   3    49C0A    49C0A "_sfree__FPv$LAJ"
        2E91   3    49E4E    49E4E "_SigmaprojTrTX__FPfN21iN24"
        2E90   3    49E4D    49E4D "_SigmaprojTrTX__FPfN21iN24$LAJ"
        3057   3    4A014    4A014 "_SigmaprojTrTY__FPfN21iN24"
        3056   3    4A013    4A013 "_SigmaprojTrTY__FPfN21iN24$LAJ"
        3131   3    4A0EE    4A0EE "_SigmaprojTrTZ__FPfN21iN24"
        3130   3    4A0ED    4A0ED "_SigmaprojTrTZ__FPfN21iN24$LAJ"
        2E20   3    49DDD    49DDD "_SigmaprojTrXT__FPfN21iN24"
        2E1F   3    49DDC    49DDC "_SigmaprojTrXT__FPfN21iN24$LAJ"
        2C6C   3    49C29    49C29 "_SigmaprojTrXY__FPfN21iN24"
        2C6B   3    49C28    49C28 "_SigmaprojTrXY__FPfN21iN24$LAJ"
        2D3E   3    49CFB    49CFB "_SigmaprojTrXZ__FPfN21iN24"
        2D3D   3    49CFA    49CFA "_SigmaprojTrXZ__FPfN21iN24$LAJ"
        2FE6   3    49FA3    49FA3 "_SigmaprojTrYT__FPfN21iN24"
        2FE5   3    49FA2    49FA2 "_SigmaprojTrYT__FPfN21iN24$LAJ"
        2CD5   3    49C92    49C92 "_SigmaprojTrYX__FPfN21iN24"
        2CD4   3    49C91    49C91 "_SigmaprojTrYX__FPfN21iN24$LAJ"
        2F03   3    49EC0    49EC0 "_SigmaprojTrYZ__FPfN21iN24"
        2F02   3    49EBF    49EBF "_SigmaprojTrYZ__FPfN21iN24$LAJ"
        30C8   3    4A085    4A085 "_SigmaprojTrZT__FPfN21iN24"
        30C7   3    4A084    4A084 "_SigmaprojTrZT__FPfN21iN24$LAJ"
        2DAF   3    49D6C    49D6C "_SigmaprojTrZX__FPfN21iN24"
        2DAE   3    49D6B    49D6B "_SigmaprojTrZX__FPfN21iN24$LAJ"
        2F75   3    49F32    49F32 "_SigmaprojTrZY__FPfN21iN24"
        2F74   3    49F31    49F31 "_SigmaprojTrZY__FPfN21iN24$LAJ"
        77FC   3    4E7B9    4E7B9 "_site2cram"
        77FD   3    4E7BA    4E7BA "_site2cram$LAJ"
        780E   3    4E7CB    4E7CB "_site2dram"
        780F   3    4E7CC    4E7CC "_site2dram$LAJ"
        9D2B   3    50CE8    50CE8 "_SiteCloverMat__13DiracOpCloverCFPCiPf"
        9D2A   3    50CE7    50CE7 "_SiteCloverMat__13DiracOpCloverCFPCiPf$LAJ"
        9BA8   3    50B65    50B65 "_SiteFuv__13DiracOpCloverCFR6MatrixPCiiT3"
        9BA7   3    50B64    50B64
				 "_SiteFuv__13DiracOpCloverCFR6MatrixPCiiT3$LAJ"
        DDEE   3    54DAB    54DAB "_SizeT"
        DDEF   3    54DAC    54DAC "_SizeX"
        DDF0   3    54DAD    54DAD "_SizeY"
        DDF1   3    54DAE    54DAE "_SizeZ"
        1D3E   3    48CFB    48CFB "_Smalloc__7VerboseFPcN21Pvi"
        1D3D   3    48CFA    48CFA "_Smalloc__7VerboseFPcN21Pvi$LAJ"
        2C2B   3    49BE8    49BE8 "_smalloc__Fi"
        2C2A   3    49BE7    49BE7 "_smalloc__Fi$LAJ"
        4341   3    4B2FE    4B2FE "_SoCheck__7LatticeF6rfloat"
        4340   3    4B2FD    4B2FD "_SoCheck__7LatticeF6rfloat$LAJ"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  51
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        4FB4   3    4BF71    4BF71 "_SpinComponents__5FstagFv"
        4FB3   3    4BF70    4BF70 "_SpinComponents__5FstagFv$LAJ"
        60BF   3    4D07C    4D07C "_SpinComponents__7FcloverFv"
        60BE   3    4D07B    4D07B "_SpinComponents__7FcloverFv$LAJ"
        27B3   3    49770    49770 "_sprojTrTm__FPfN21iN24"
        27B2   3    4976F    4976F "_sprojTrTm__FPfN21iN24$LAJ"
        2A7E   3    49A3B    49A3B "_sprojTrTp__FPfN21iN24"
        2A7D   3    49A3A    49A3A "_sprojTrTp__FPfN21iN24$LAJ"
        28D1   3    4988E    4988E "_sprojTrXm__FPfN21iN24"
        28D0   3    4988D    4988D "_sprojTrXm__FPfN21iN24$LAJ"
        2B9C   3    49B59    49B59 "_sprojTrXp__FPfN21iN24"
        2B9B   3    49B58    49B58 "_sprojTrXp__FPfN21iN24$LAJ"
        2B0D   3    49ACA    49ACA "_sprojTrYm__FPfN21iN24"
        2B0C   3    49AC9    49AC9 "_sprojTrYm__FPfN21iN24$LAJ"
        29EF   3    499AC    499AC "_sprojTrYp__FPfN21iN24"
        29EE   3    499AB    499AB "_sprojTrYp__FPfN21iN24$LAJ"
        2960   3    4991D    4991D "_sprojTrZm__FPfN21iN24"
        295F   3    4991C    4991C "_sprojTrZm__FPfN21iN24$LAJ"
        2842   3    497FF    497FF "_sprojTrZp__FPfN21iN24"
        2841   3    497FE    497FE "_sprojTrZp__FPfN21iN24$LAJ"
        EC83   3    55C40    55C40 "_sqrf40"
        EB78   3    55B35    55B35 "_sqrf64"
        F4D7   3    56494    56494 "_sqrt"
        EE13   3    55DD0    55DD0 "_sqrtf64"
        E83D   3    557FA    557FA "_sqrt__FRC8double64"
        E83C   3    557F9    557F9 "_sqrt__FRC8double64$LAJ"
        2320   1     3320     3320 "_srbuf_offset"
        DE0A   3    54DC7    54DC7 "_sscanf"
        789F   3    4E85C    4E85C "_StagSCUCommBackward"
        789E   3    4E85B    4E85B "_StagSCUCommBackward$LAJ"
        7865   3    4E822    4E822 "_StagSCUCommForward"
        7864   3    4E821    4E821 "_StagSCUCommForward$LAJ"
        78D9   3    4E896    4E896 "_StagSCUComplete"
        78D8   3    4E895    4E895 "_StagSCUComplete$LAJ"
        7823   3    4E7E0    4E7E0 "_StagSCUSetup"
        7822   3    4E7DF    4E7DF "_StagSCUSetup$LAJ"
           0   8   809800    56D6D "_stag_ds_dest"
        3A74   3    4AA31    4AA31 "_Staple__7LatticeFR6MatrixPii"
        3A73   3    4AA30    4AA30 "_Staple__7LatticeFR6MatrixPii$LAJ"
        E5BF   3    5557C    5557C "_strlen"
        E5C0   3    5557D    5557D "_strlen$LAJ"
        E59F   3    5555C    5555C "_strncpy"
        E5A0   3    5555D    5555D "_strncpy$LAJ"
        3A58   3    4AA15    4AA15 "_StrOrd__7LatticeFv"
        3A57   3    4AA14    4AA14 "_StrOrd__7LatticeFv$LAJ"
         DFB   1     1DFB     1DFB "_str_ord__7Lattice"
        EB6B   3    55B28    55B28 "_subf32f40"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  52
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        F435   3    563F2    563F2 "_subf32f64"
        F44F   3    5640C    5640C "_subf64"
        F441   3    563FE    563FE "_subf64f32"
        3CA7   3    4AC64    4AC64 "_SumReTrPlaqNode__7LatticeCFv"
        3CA6   3    4AC63    4AC63 "_SumReTrPlaqNode__7LatticeCFv$LAJ"
        3D39   3    4ACF6    4ACF6 "_SumReTrPlaq__7LatticeCFv"
        3D38   3    4ACF5    4ACF5 "_SumReTrPlaq__7LatticeCFv$LAJ"
        DDF9   3    54DB6    54DB6 "_sync"
        F008   3    55FC5    55FC5 "_tentof64"
        20F8   3    490B5    490B5 "_Trans__6MatrixFPCf"
        20F7   3    490B4    490B4 "_Trans__6MatrixFPCf$LAJ"
        26BD   3    4967A    4967A "_TrLessAntiHermMatrix__6MatrixFRC6Matrix"
        26BC   3    49679    49679
				  "_TrLessAntiHermMatrix__6MatrixFRC6Matrix$LAJ"
        EB56   3    55B13    55B13 "_trunci32f64"
        1F38   3    48EF5    48EF5 "_uDotXEqual"
        DDE9   3    54DA6    54DA6 "_UniqueID"
        269F   3    4965C    4965C "_Unitarize__6MatrixFv"
        269E   3    4965B    4965B "_Unitarize__6MatrixFv$LAJ"
        2182   3    4913F    4913F "_UnitMatrix__6MatrixFv"
        2183   3    49140    49140 "_UnitMatrix__6MatrixFv$LAJ"
        1FD4   3    48F91    48F91 "_vecAddEquVec"
        1FE2   3    48F9F    48F9F "_vecMinusEquVec"
        1FF0   3    48FAD    48FAD "_vecNegative"
        1FC5   3    48F82    48F82 "_vecTimesEquFloat"
           0   4   809F95    567FF "_vector_dest"
         22B   1     122B     122B "_VRB"
        DE39   3    54DF6    54DF6 "_vsprintf"
        DE38   3    54DF5    54DF5 "_vsprintf$LAJ"
        1E06   3    48DC3    48DC3 "_Warn__7VerboseFPcT1PCce"
        1E05   3    48DC2    48DC2 "_Warn__7VerboseFPcT1PCce$LAJ"
           0   5   809800    5682B "_wfm0_dest"
           0   6   809C00    56B07 "_wfm1_dest"
         20C   5   809A0C    56A37 "_wfm_cmat_spproj"
         1FE   5   8099FE    56A29 "_wfm_cmat_spproj_count"
        AC27   3    51BE4    51BE4 "_wfm_comm_backward"
        ABD0   3    51B8D    51B8D "_wfm_comm_forward"
        B1F7   3    521B4    521B4 "_wfm_copy_backward"
        B1F6   3    521B3    521B3 "_wfm_copy_backward$LAJ"
        B1B1   3    5216E    5216E "_wfm_copy_forward"
        B1B0   3    5216D    5216D "_wfm_copy_forward$LAJ"
        ABA3   3    51B60    51B60 "_wfm_dslash"
         108   5   809908    56933 "_wfm_mat_trick"
          F3   5   8098F3    5691E "_wfm_mat_trick_count"
        AB24   3    51AE1    51AE1 "_wfm_scu_init"
        AAE5   3    51AA2    51AA2 "_wfm_scu_wait"
          60   6   809C60    56B67 "_wfm_spproj"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  53
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        AAC7   3    51A84    51A84 "_wfm_spproj_segment"
        AE84   3    51E41    51E41 "_wfm_sublatt_pointers__FiN41P6Wilson"
        AE83   3    51E40    51E40 "_wfm_sublatt_pointers__FiN41P6Wilson$LAJ"
          13   6   809C13    56B1A "_wfm_trick_segment_1"
          38   6   809C38    56B3F "_wfm_trick_segment_2"
          4B   6   809C4B    56B52 "_wfm_trick_segment_3"
           2   5   809802    5682D "_wfm_trick_spproj"
        2702   1     3702     3702 "_wfm_wire_map"
        AA75   3    51A32    51A32 "_wilson_dslash"
        ADDC   3    51D99    51D99 "_wilson_end__FP6Wilson"
        ADDB   3    51D98    51D98 "_wilson_end__FP6Wilson$LAJ"
        ACD4   3    51C91    51C91 "_wilson_init__FP6Wilson"
        ACD3   3    51C90    51C90 "_wilson_init__FP6Wilson$LAJ"
        312C   1     412C     412C "_XbigDbl"
        2190   3    4914D    4914D "_ZeroMatrix__6MatrixFv"
        2191   3    4914E    4914E "_ZeroMatrix__6MatrixFv$LAJ"
        E5E6   3    555A3    555A3 "__array_pointer_not_from_vec_new"
        E5E5   3    555A2    555A2 "__array_pointer_not_from_vec_new$LAJ"
        E694   3    55651    55651 "__ARTDHALT"
        E4D1   3    5548E    5548E "__ConvertDouble"
        E4D0   3    5548D    5548D "__ConvertDouble$LAJ"
        DE30   3    54DED    54DED "__ConvertF64"
        E586   3    55543    55543 "__Exit"
        F50D   3    564CA    564CA "__log"
        EB43   3    55B00    55B00 "__main"
        EB44   3    55B01    55B01 "__main$LAJ"
        E12D   3    550EA    550EA "__pf"
        E12C   3    550E9    550E9 "__pf$LAJ"
        E393   3    55350    55350 "__RealConvF64"
        E392   3    5534F    5534F "__RealConvF64$LAJ"
        E473   3    55430    55430 "__TimesPower"
        E472   3    5542F    5542F "__TimesPower$LAJ"
        319A   3    4A157    4A157 "___adv__6rfloatFf"
        3199   3    4A156    4A156 "___adv__6rfloatFf$LAJ"
        31AF   3    4A16C    4A16C "___ami__6rfloatFf"
        31B0   3    4A16D    4A16D "___ami__6rfloatFf$LAJ"
        33BD   3    4A37A    4A37A "___ami__8RcomplexFRC8Rcomplex"
        33BE   3    4A37B    4A37B "___ami__8RcomplexFRC8Rcomplex$LAJ"
        31A6   3    4A163    4A163 "___amu__6rfloatFf"
        31A7   3    4A164    4A164 "___amu__6rfloatFf$LAJ"
        337B   3    4A338    4A338 "___amu__8RcomplexFf"
        337C   3    4A339    4A339 "___amu__8RcomplexFf$LAJ"
        33A1   3    4A35E    4A35E "___amu__8RcomplexFRC8Rcomplex"
        33A2   3    4A35F    4A35F "___amu__8RcomplexFRC8Rcomplex$LAJ"
        31B8   3    4A175    4A175 "___apl__6rfloatFf"
        31B9   3    4A176    4A176 "___apl__6rfloatFf$LAJ"
        20E4   3    490A1    490A1 "___as__6MatrixFf"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  54
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        20E3   3    490A0    490A0 "___as__6MatrixFf$LAJ"
        3425   3    4A3E2    4A3E2 "___as__8RcomplexFRC8Rcomplex"
        3426   3    4A3E3    4A3E3 "___as__8RcomplexFRC8Rcomplex$LAJ"
        4C05   3    4BBC2    4BBC2 "___ct__10FstagTypesFv"
        4C04   3    4BBC1    4BBC1 "___ct__10FstagTypesFv$LAJ"
        9208   3    501C5    501C5
		   "___ct__11DiracOpStagFR7LatticeP6VectorT2P5CgArg10CnvFrmType"
        9207   3    501C4    501C4
	       "___ct__11DiracOpStagFR7LatticeP6VectorT2P5CgArg10CnvFrmType$LAJ"
        4ADF   3    4BA9C    4BA9C "___ct__12FwilsonTypesFv"
        4ADE   3    4BA9B    4BA9B "___ct__12FwilsonTypesFv$LAJ"
        3562   3    4A51F    4A51F "___ct__12GwilsonFstagFv"
        3561   3    4A51E    4A51E "___ct__12GwilsonFstagFv$LAJ"
        9F66   3    50F23    50F23
		 "___ct__13DiracOpCloverFR7LatticeP6VectorT2P5CgArg10CnvFrmType"
        9F65   3    50F22    50F22
	     "___ct__13DiracOpCloverFR7LatticeP6VectorT2P5CgArg10CnvFrmType$LAJ"
        3613   3    4A5D0    4A5D0 "___ct__14GwilsonFcloverFv"
        3612   3    4A5CF    4A5CF "___ct__14GwilsonFcloverFv$LAJ"
        99DE   3    5099B    5099B
	      "___ct__16DiracOpStagTypesFR7LatticeP6VectorT2P5CgArg10CnvFrmType"
        99DD   3    5099A    5099A
	  "___ct__16DiracOpStagTypesFR7LatticeP6VectorT2P5CgArg10CnvFrmType$LAJ"
        D754   3    54711    54711
	    "___ct__18DiracOpWilsonTypesFR7LatticeP6VectorT2P5CgArg10CnvFrmType"
        D753   3    54710    54710
	"___ct__18DiracOpWilsonTypesFR7LatticeP6VectorT2P5CgArg10CnvFrmType$LAJ"
        443A   3    4B3F7    4B3F7 "___ct__18GlobalJobParameterFv"
        4439   3    4B3F6    4B3F6 "___ct__18GlobalJobParameterFv$LAJ"
        1B5D   3    48B1A    48B1A "___ct__3AlgFR7LatticeP9CommonArg"
        1B5C   3    48B19    48B19 "___ct__3AlgFR7LatticeP9CommonArg$LAJ"
        6B29   3    4DAE6    4DAE6 "___ct__5ErrorFv"
        6B28   3    4DAE5    4DAE5 "___ct__5ErrorFv$LAJ"
        4DC1   3    4BD7E    4BD7E "___ct__5FstagFv"
        4DC0   3    4BD7D    4BD7D "___ct__5FstagFv$LAJ"
         C03   3    47BC0    47BC0 "___ct__6AlgHmdFR7LatticeP9CommonArgP6HmdArg"
         C02   3    47BBF    47BBF
			       "___ct__6AlgHmdFR7LatticeP9CommonArgP6HmdArg$LAJ"
        20BA   3    49077    49077 "___ct__6MatrixFRC6Matrix"
        20B9   3    49076    49076 "___ct__6MatrixFRC6Matrix$LAJ"
        209A   3    49057    49057 "___ct__6MatrixFv"
        2099   3    49056    49056 "___ct__6MatrixFv$LAJ"
        3339   3    4A2F6    4A2F6 "___ct__6rfloatFf"
        3338   3    4A2F5    4A2F5 "___ct__6rfloatFf$LAJ"
        334F   3    4A30C    4A30C "___ct__6rfloatFRC6rfloat"
        334E   3    4A30B    4A30B "___ct__6rfloatFRC6rfloat$LAJ"
        219B   3    49158    49158 "___ct__6VectorFv"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  55
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        219A   3    49157    49157 "___ct__6VectorFv$LAJ"
         3E7   3    473A4    473A4
				  "___ct__7AlgHmdRFR7LatticeP9CommonArgP6HmdArg"
         3E6   3    473A3    473A3
			      "___ct__7AlgHmdRFR7LatticeP9CommonArgP6HmdArg$LAJ"
        CD6E   3    53D2B    53D2B
			"___ct__7DiracOpFR7LatticeP6VectorT2P5CgArg10CnvFrmType"
        CD6D   3    53D2A    53D2A
		    "___ct__7DiracOpFR7LatticeP6VectorT2P5CgArg10CnvFrmType$LAJ"
        5EF5   3    4CEB2    4CEB2 "___ct__7FcloverFv"
        5EF4   3    4CEB1    4CEB1 "___ct__7FcloverFv$LAJ"
        4868   3    4B825    4B825 "___ct__7GwilsonFv"
        4867   3    4B824    4B824 "___ct__7GwilsonFv$LAJ"
        37B0   3    4A76D    4A76D "___ct__7LatticeFv"
        37AF   3    4A76C    4A76C "___ct__7LatticeFv$LAJ"
        1BCB   3    48B88    48B88 "___ct__7VerboseFv"
        1BCA   3    48B87    48B87 "___ct__7VerboseFv$LAJ"
        33E9   3    4A3A6    4A3A6 "___ct__8RcomplexFfT1"
        33E8   3    4A3A5    4A3A5 "___ct__8RcomplexFfT1$LAJ"
        3401   3    4A3BE    4A3BE "___ct__8RcomplexFRC8Rcomplex"
        3400   3    4A3BD    4A3BD "___ct__8RcomplexFRC8Rcomplex$LAJ"
         CB2   3    47C6F    47C6F
				"___ct__9AlgHmcPhiFR7LatticeP9CommonArgP6HmdArg"
         CB1   3    47C6E    47C6E
			    "___ct__9AlgHmcPhiFR7LatticeP9CommonArgP6HmdArg$LAJ"
        DD74   3    54D31    54D31 "___ct__9SCUDirArgFPv6SCUDir5SCUXRiN24"
        DD73   3    54D30    54D30 "___ct__9SCUDirArgFPv6SCUDir5SCUXRiN24$LAJ"
        DD68   3    54D25    54D25 "___ct__9SCUDirArgFv"
        DD67   3    54D24    54D24 "___ct__9SCUDirArgFv$LAJ"
        E6F1   3    556AE    556AE "___dl__FPv"
        E6F0   3    556AD    556AD "___dl__FPv$LAJ"
        4C41   3    4BBFE    4BBFE "___dt__10FstagTypesFv"
        4C40   3    4BBFD    4BBFD "___dt__10FstagTypesFv$LAJ"
        929F   3    5025C    5025C "___dt__11DiracOpStagFv"
        929E   3    5025B    5025B "___dt__11DiracOpStagFv$LAJ"
        4B7F   3    4BB3C    4BB3C "___dt__12FwilsonTypesFv"
        4B7E   3    4BB3B    4BB3B "___dt__12FwilsonTypesFv$LAJ"
        35BF   3    4A57C    4A57C "___dt__12GwilsonFstagFv"
        35BE   3    4A57B    4A57B "___dt__12GwilsonFstagFv$LAJ"
        9FE1   3    50F9E    50F9E "___dt__13DiracOpCloverFv"
        9FE0   3    50F9D    50F9D "___dt__13DiracOpCloverFv$LAJ"
        3670   3    4A62D    4A62D "___dt__14GwilsonFcloverFv"
        366F   3    4A62C    4A62C "___dt__14GwilsonFcloverFv$LAJ"
        9A17   3    509D4    509D4 "___dt__16DiracOpStagTypesFv"
        9A16   3    509D3    509D3 "___dt__16DiracOpStagTypesFv$LAJ"
        D78D   3    5474A    5474A "___dt__18DiracOpWilsonTypesFv"
        D78C   3    54749    54749 "___dt__18DiracOpWilsonTypesFv$LAJ"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  56
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        4495   3    4B452    4B452 "___dt__18GlobalJobParameterFv"
        4494   3    4B451    4B451 "___dt__18GlobalJobParameterFv$LAJ"
        1B9D   3    48B5A    48B5A "___dt__3AlgFv"
        1B9C   3    48B59    48B59 "___dt__3AlgFv$LAJ"
        6B5B   3    4DB18    4DB18 "___dt__5ErrorFv"
        6B5A   3    4DB17    4DB17 "___dt__5ErrorFv$LAJ"
        4F32   3    4BEEF    4BEEF "___dt__5FstagFv"
        4F31   3    4BEEE    4BEEE "___dt__5FstagFv$LAJ"
         C74   3    47C31    47C31 "___dt__6AlgHmdFv"
         C73   3    47C30    47C30 "___dt__6AlgHmdFv$LAJ"
        3366   3    4A323    4A323 "___dt__6rfloatFv"
        3365   3    4A322    4A322 "___dt__6rfloatFv$LAJ"
         596   3    47553    47553 "___dt__7AlgHmdRFv"
         595   3    47552    47552 "___dt__7AlgHmdRFv$LAJ"
        CEB0   3    53E6D    53E6D "___dt__7DiracOpFv"
        CEAF   3    53E6C    53E6C "___dt__7DiracOpFv$LAJ"
        5FB4   3    4CF71    4CF71 "___dt__7FcloverFv"
        5FB3   3    4CF70    4CF70 "___dt__7FcloverFv$LAJ"
        48A4   3    4B861    4B861 "___dt__7GwilsonFv"
        48A3   3    4B860    4B860 "___dt__7GwilsonFv$LAJ"
        3991   3    4A94E    4A94E "___dt__7LatticeFv"
        3990   3    4A94D    4A94D "___dt__7LatticeFv$LAJ"
        1BF8   3    48BB5    48BB5 "___dt__7VerboseFv"
        1BF7   3    48BB4    48BB4 "___dt__7VerboseFv$LAJ"
        341A   3    4A3D7    4A3D7 "___dt__8RcomplexFv"
        3419   3    4A3D6    4A3D6 "___dt__8RcomplexFv$LAJ"
        105F   3    4801C    4801C "___dt__9AlgHmcPhiFv"
        105E   3    4801B    4801B "___dt__9AlgHmcPhiFv$LAJ"
        DD97   3    54D54    54D54 "___dt__9SCUDirArgFv"
        DD96   3    54D53    54D53 "___dt__9SCUDirArgFv$LAJ"
        32FD   3    4A2BA    4A2BA "___dv__FdRC6rfloat"
        32FC   3    4A2B9    4A2B9 "___dv__FdRC6rfloat$LAJ"
        331B   3    4A2D8    4A2D8 "___dv__FRC6rfloatd"
        331A   3    4A2D7    4A2D7 "___dv__FRC6rfloatd$LAJ"
        32DE   3    4A29B    4A29B "___dv__FRC6rfloatT1"
        32DD   3    4A29A    4A29A "___dv__FRC6rfloatT1$LAJ"
        3247   3    4A204    4A204 "___mi__FdRC6rfloat"
        3246   3    4A203    4A203 "___mi__FdRC6rfloat$LAJ"
        31C2   3    4A17F    4A17F "___mi__FRC6rfloat"
        31C1   3    4A17E    4A17E "___mi__FRC6rfloat$LAJ"
        3265   3    4A222    4A222 "___mi__FRC6rfloatd"
        3264   3    4A221    4A221 "___mi__FRC6rfloatd$LAJ"
        3228   3    4A1E5    4A1E5 "___mi__FRC6rfloatT1"
        3227   3    4A1E4    4A1E4 "___mi__FRC6rfloatT1$LAJ"
        3387   3    4A344    4A344 "___mi__FRC8Rcomplex"
        3386   3    4A343    4A343 "___mi__FRC8Rcomplex$LAJ"
        32A2   3    4A25F    4A25F "___ml__FdRC6rfloat"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  57
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

        32A1   3    4A25E    4A25E "___ml__FdRC6rfloat$LAJ"
        32C0   3    4A27D    4A27D "___ml__FRC6rfloatd"
        32BF   3    4A27C    4A27C "___ml__FRC6rfloatd$LAJ"
        3283   3    4A240    4A240 "___ml__FRC6rfloatT1"
        3282   3    4A23F    4A23F "___ml__FRC6rfloatT1$LAJ"
        E6D0   3    5568D    5568D "___nw__FUi"
        E6CF   3    5568C    5568C "___nw__FUi$LAJ"
        31EC   3    4A1A9    4A1A9 "___pl__FdRC6rfloat"
        31EB   3    4A1A8    4A1A8 "___pl__FdRC6rfloat$LAJ"
        320A   3    4A1C7    4A1C7 "___pl__FRC6rfloatd"
        3209   3    4A1C6    4A1C6 "___pl__FRC6rfloatd$LAJ"
        31CD   3    4A18A    4A18A "___pl__FRC6rfloatT1"
        31CC   3    4A189    4A189 "___pl__FRC6rfloatT1$LAJ"
        E6CE   3    5568B    5568B "___pure_virtual_called"
        E6CD   3    5568A    5568A "___pure_virtual_called$LAJ"
        E626   3    555E3    555E3 "___vec_delete"
        E625   3    555E2    555E2 "___vec_delete$LAJ"
        E5E8   3    555A5    555A5 "___vec_new"
        E5E7   3    555A4    555A4 "___vec_new$LAJ"
        18B7   1     28B7     28B7 "___vtbl__10FstagTypes"
         CC7   1     1CC7     1CC7 "___vtbl__10FstagTypes__12GwilsonFstag"
        193D   1     293D     293D "___vtbl__10FstagTypes__5Fstag"
        232E   1     332E     332E "___vtbl__11DiracOpStag"
        17EB   1     27EB     27EB "___vtbl__12FwilsonTypes"
         DA4   1     1DA4     1DA4 "___vtbl__12FwilsonTypes__14GwilsonFclover"
        1B5A   1     2B5A     2B5A "___vtbl__12FwilsonTypes__7Fclover"
         CD0   1     1CD0     1CD0 "___vtbl__12GwilsonFstag"
        258E   1     358E     358E "___vtbl__13DiracOpClover"
         DB0   1     1DB0     1DB0 "___vtbl__14GwilsonFclover"
         BD9   1     1BD9     1BD9 "___vtbl__15RandomGenerator"
        24BA   1     34BA     34BA "___vtbl__16DiracOpStagTypes"
        2C47   1     3C47     3C47 "___vtbl__18DiracOpWilsonTypes"
         BC7   1     1BC7     1BC7 "___vtbl__22UniformRandomGenerator"
         BD0   1     1BD0     1BD0 "___vtbl__23GaussianRandomGenerator"
         956   1     1956     1956 "___vtbl__3Alg"
        1D4E   1     2D4E     2D4E "___vtbl__5Error"
        1946   1     2946     2946 "___vtbl__5Fstag"
         C34   1     1C34     1C34 "___vtbl__5Fstag__12GwilsonFstag"
         613   1     1613     1613 "___vtbl__6AlgHmd"
         46A   1     146A     146A "___vtbl__7AlgHmdR"
        2A34   1     3A34     3A34 "___vtbl__7DiracOp"
        1B66   1     2B66     2B66 "___vtbl__7Fclover"
         D11   1     1D11     1D11 "___vtbl__7Fclover__14GwilsonFclover"
        16B1   1     26B1     26B1 "___vtbl__7Gwilson"
         DFF   1     1DFF     1DFF "___vtbl__7Lattice"
        1860   1     2860     2860 "___vtbl__7Lattice__10FstagTypes"
        1794   1     2794     2794 "___vtbl__7Lattice__12FwilsonTypes"

Tartan Linker, SPARC/C40, v5.1.0	5/31/99 17:50:02		Page  58
Copyright (c) 1986-1992 by Tartan, Inc., All Rights Reserved

         C70   1     1C70     1C70 "___vtbl__7Lattice__12GwilsonFstag"
         D4D   1     1D4D     1D4D "___vtbl__7Lattice__14GwilsonFclover"
        18E6   1     28E6     28E6 "___vtbl__7Lattice__5Fstag"
        1B03   1     2B03     2B03 "___vtbl__7Lattice__7Fclover"
        16C6   1     26C6     26C6 "___vtbl__7Lattice__7Gwilson"
         988   1     1988     1988 "___vtbl__7Verbose"
         653   1     1653     1653 "___vtbl__9AlgHmcPhi"

Error summary:
	Information: 3, Warnings: 6, Errors: 3, Fatal: 0

Tartan Linker elapsed time 3.0 seconds, CPU time 0.31 seconds
CPS_END_NAMESPACE
