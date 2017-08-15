#Celular Automata for Fluid Charge Adjustments
from optparse import OptionParser

#Opciones del script
parser=OptionParser()
parser.add_option("-i","--input", action="store", type="string", dest="InputFile", default='HCl.mol2', help="Path of the .mol2 input file.")
parser.add_option("-L","--realdistance", action="store", type="float", dest="Distance", default=0.5, help="Real-space distance from atoms to simulation box (nm).")
parser.add_option("-V","--veldistance", action="store", type="float", dest="Distancev", default=0.5, help="Velocity distance from atoms to simulation box (nm/ps).")
parser.add_option("-l","--dx", action="store", type="float", dest="Dx", default=0.01, help="Resolution of the real space (nm)")
parser.add_option("-v","--dv", action="store", type="float", dest="Dv", default=0.01, help="Resolution of the velocity space (nm/ps)")


(options, args) = parser.parse_args()

inp_name=options.InputFile
distx=options.Distance
distv=options.Distancev
dx=options.Dx
dv=options.Dv
