struct CommonArg {

    void *results;  /*!< This can point to a string containing the name
		      of a file to which to write results, or an object to
    		    which results might be passed.*/

    char *filename;   /*!< The name of a file to which to write results. */
    char *label;       /*!< A label which uniquely identifies this task. */ 
};

