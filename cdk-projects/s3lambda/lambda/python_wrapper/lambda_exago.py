import lambda_exago
exago.initialize("opflow")
opf = exago.OPFLOW()
opf.read_mat_power_data('datafiles/case9/case9mod.m')
opf.solve()
opf.print_solution()
del opf
exago.finalize()
