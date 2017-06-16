all:compgeom_.so

compgeom_.so:
	g++ -shared -fPIC -o compgeom_.so \
		-I /home/vlad/PyEnvs/py2dp/include/python2.7 \
		-std=c++0x compgeom.cxx -lboost_python -lCGAL \
		-lCGAL_Core -lCGAL_Kernel_cpp -lgmp -lmpfr

clean:
	rm -Rf compgeom_.so

.PHONY: clean all

