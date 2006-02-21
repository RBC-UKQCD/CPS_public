
/*! A structure to hold the solver parameters.*/
class CgArg {

    Float mass;			/*!<  The mass parameter. */

    int max_num_iter;		/*!<  The maximum number of solver
				 iterations to do. */

    Float stop_rsd;		/*!<  The target residual. */
	Float true_rsd;

    enum RitzMatType RitzMatOper; /*!< Which operator to determine eigenvalues
				    of, if any. */
    enum InverterType Inverter;     /*!< Which solver to use.*/

    int bicgstab_n;             /* BiCGstab(n) parameter. */

    memfun CgArg();
};

