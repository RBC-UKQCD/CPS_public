set mode='co'
echo pushd src/comms/bgl/scu/
pushd src/comms/bgl/scu/
echo cp bgl_net.C.$mode bgl_net.C
cp bgl_net.C.$mode bgl_net.C
echo cp bgl_sys.C.$mode bgl_sys.C
cp bgl_sys.C.$mode bgl_sys.C
echo popd
popd
echo pushd src/util/dirac_op/d_op_wilson/bgl/
pushd src/util/dirac_op/d_op_wilson/bgl/
echo cp wfm_comm.C.$mode wfm_comm.C
cp wfm_comm.C.$mode wfm_comm.C
echo popd
popd
#echo cp Makefile.rules.$mode Makefiles.rules
#cp Makefile.rules.$mode Makefile.rules
#cp Makefile.rules.$mode ../ANL/Makefile.rules

