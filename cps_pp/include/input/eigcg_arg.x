/*------------------------------------------------------------------*/
/*!\file
  \brief  Definitions of the EigCGArg structure.
*/
/*------------------------------------------------------------------*/

/*! A structure holding the parameters relevant to the eigenvalue measurements in eigCG.*/

class EigCGArg {

    int nev; /*number of eigen vectos to solve in each step. it generate 2*nev low modes*/
	int m; /*restart steps, should bigger than 2nev */
	int max_def_len; /*maximum deflaiton length, the code will create such a space to store it */
	Float max_eig_cut; /*store only the eigen value that is smaller than this value*/
	bool always_restart; /*if true, then it will do CG restart even when accumulating low modes */
	int restart_len; /*number of restarts*/
	Float restart[10]; /* restart the init-CG stop condition */

    memfun EigCGArg();
    memfun ~EigCGArg();
};


