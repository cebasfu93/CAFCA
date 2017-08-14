CAFCA3d.gif CAFCA.gif temp/ temp3d/ : test_plot.py potx.outc poty.outc potz.outc rspace.outc accx.outc accy.outc accz.outc densx.outc densy.outc densz.outc constants.outc
	#python plot.py
	rm -rf temp/ temp3d/
	python test_plot.py

constants.outc rspace.outc densx.outc densy.outc densz.outc accx.outc accy.outc accz.outc potx.outc poty.outc potz.outc : LB3D.x
	./LB3D.x

LB3D.x : LB3D.c funciones.h funciones.c constantes.h atomos.outpy constants.outpy
	gcc -lm -lfftw3 LB3D.c -o LB3D.x

atomos.outpy constants.outpy : init.py atom.py molecule.py
	python init.py -i HCl.mol2
	#python init.py -i MET.mol2
	#python init.py -i ACT.mol2
	python molecule.py

clean:
	rm -f atomos.outpy constants.outpy
	rm -f atom.pyc init.pyc
	rm -f atomos.outc potx.outc poty.outc potz.outc accx.outc accy.outc accz.outc densx.outc densy.outc densz.outc rspace.outc constants.outc
	rm -f LB3D.x
	rm -f CAFCA.gif CAFCA3d.gif
	rm -rf temp/ temp3d/
