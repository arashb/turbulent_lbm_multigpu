all:
	scons

doc:
	doxygen ./docs/Doxyfile

cleandoc:
	rm -rf ./docs/html ./docs/latex

cleanvtk:
	rm -rf vtkOutput/*.vtk

clean:
	rm -rf build/* vtkOutput/*.vtk
