struct EqStateArg {
  int dir;      /* the special direction which determines what hyperplane(s)
		 * the sum of plaq is applied on. There are two sum
		 * operations, one on the hyperplane perpendicular to this
		 * direction, the other on all the hyperplanes parallel to
		 * this direction.
		 */
};
