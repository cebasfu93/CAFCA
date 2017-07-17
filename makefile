LB3D.x : LB3D.c
	gcc -lm -lfftw3 LB3D.c -o LB3D.x
	./LB3D.x

coords.txt charges.txt types.txt names.txt : init.py
	python init.py -i HCl.mol2

clean:
	rm -f coords.txt charges.txt types.txt names.txt LB3D.x
