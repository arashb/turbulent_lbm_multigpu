all:
	scons

doc:
	doxygen ./docs/Doxyfile

cleandoc:
	rm -rf ./docs/html ./docs/latex

cleanvtk:
	rm  vtkOutput/*.vtk

cleanbenchmark:
	rm  benchmarkOutput/*.ini

cleanprofile:
	rm  profileOutput/*.ini

clean:
	rm -rf build/* vtkOutput/*.vtk benchmarkOutput/*.ini profileOutput/*.ini
