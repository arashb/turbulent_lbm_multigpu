all:
	scons

doc:
	doxygen ./docs/Doxyfile

cleandoc:
	rm -rf ./docs/html ./docs/latex

cleanvtk:
	rm  ./output/vtk/*.vtk

cleanbenchmark:
	rm  ./output/benchmark/*.ini

cleanprofile:
	rm  ./output/profile/*.ini

cleanlog: 
	rm  ./output/log/*.log

cleanbuild:
	rm -rf build/*

clean: cleanbuild cleanlog cleanprofile cleanbenchmark cleanvtk  
