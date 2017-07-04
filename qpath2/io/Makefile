all: io_.so

io_.so: io_.cxx
	g++ -shared -fPIC -o io_.so \
		-I /home/vlad/PyEnvs/py2dp/include/python2.7 \
		`pkg-config --cflags openslide` \
		-std=c++0x io_.cxx -lboost_python \
		`pkg-config --libs openslide`


clean:
	rm -Rf io_.so

.PHONY: clean all

