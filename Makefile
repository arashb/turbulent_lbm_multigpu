all:
	scons

doc:
	doxygen ./docs/Doxyfile

cleandoc:
	rm -rf ./docs/html ./docs/latex

cleanvtk:
	rm -f ./output/vtk/*

cleanbenchmark:
	rm -f ./output/benchmark/*

cleanprofile:
	rm -f ./output/profile/*

cleanlog: 
	rm -f ./output/log/*

cleanbuild:
	rm -rf build/*

cleanresultsoutput:
	rm -rf output

clean: cleanbuild cleanlog cleanprofile cleanbenchmark cleanvtk 
