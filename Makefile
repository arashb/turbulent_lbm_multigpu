all:
	scons

doc:
	doxygen ./docs/Doxyfile

cleandoc:
	rm -rf ./docs/html ./docs/latex

clean:
	rm -rf build/* vtkOutput/*
