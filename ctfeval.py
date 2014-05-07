#!/usr/bin/env python

#pythonlib
import os
#appion
from appionlib import basicScript
from appionlib import apDisplay
from appionlib.apCtf import ctfdisplay

class CtfEval(basicScript.BasicScript):
	"""
	appion Loop function that
	runs Craig's ace2 program
	to estimate the CTF in images
	"""

	#======================
	def start(self):
		#def runCTFdisplayTools(imgdata, ctfvalues, opimagedir, fftpath=None, fftfreq=None):
		### RUN CTF DISPLAY TOOLS
		ctfdisplaydict = ctfdisplay.makeCtfImages(imgdata, ctfvalues, fftpath, fftfreq)
		if ctfdisplaydict is None:
			raise
		### save the classic images as well
		### new powerspec file
		psfile = os.path.join(opimagedir, ctfdisplaydict['powerspecfile'])
		if not os.path.isfile(ctfdisplaydict['powerspecfile']):
			apDisplay.printWarning("Powerspec file not created")
		else:
			print ctfdisplaydict['powerspecfile']
			shutil.move(ctfdisplaydict['powerspecfile'], psfile)
			ctfvalues['graph1'] = os.path.basename(psfile)
		### new 1d plot file
		plotfile = os.path.join(opimagedir, ctfdisplaydict['plotsfile'])
		shutil.move(ctfdisplaydict['plotsfile'], plotfile)
		ctfvalues['graph2'] = os.path.basename(plotfile)
		ctfvalues['confidence_30_10'] = ctfdisplaydict['conf3010']
		ctfvalues['confidence_5_peak'] = ctfdisplaydict['conf5peak']
		ctfvalues['overfocus_conf_30_10'] = ctfdisplaydict['overconf3010']
		ctfvalues['overfocus_conf_5_peak'] = ctfdisplaydict['overconf5peak']
		ctfvalues['resolution_80_percent'] = ctfdisplaydict['res80']
		ctfvalues['resolution_50_percent'] = ctfdisplaydict['res50']
		if not 'confidence_d' in ctfvalues or ctfvalues['confidence_d'] is None:
			ctfvalues['confidence_d'] = ctfdisplaydict['conf5peak']
		if not 'confidence' in ctfvalues or ctfvalues['confidence'] is None:
			ctfvalues['confidence'] = ctfdisplaydict['conf3010']

		### override the confidence
		ctfvalues['confidence'] = max(ctfvalues['confidence'], ctfvalues['confidence_d'], ctfdisplaydict['conf5peak'], ctfdisplaydict['conf3010'])

		return ctfvalues


	#======================
	def setupParserOptions(self):
		### values
		self.parser.add_option("-i", "--image", dest="imagefile",
			help="Image file for processing (MRC preferred)", metavar="FILE")
		self.parser.add_option("-p", "--apix", dest="apix", type="float",
			help="Pixel size of image in Angstroms", metavar="##")
		self.parser.add_option("-s", "--cs", dest="cs", type="float",
			help="Spherical aberration in millimeters (mm)", metavar="##")
		self.parser.add_option("-v", "--kv", "--voltage", dest="kv", type="int",
			help="Voltage potential of microscope in kiloVolts (kV)", metavar="##")
		self.parser.add_option("-r", "--reslimit", dest="reslimit", default=5.0, type="float",
			help="Maximum resolution to check ctf in Angstroms (default 5 A)", metavar="##")

		self.parser.add_option("-1", "--defocus1", "--mindef", dest="defocus1", type="float",
			help="Minimum defocus 1 (major axis) in microns", metavar="##")
		self.parser.add_option("-2", "--defocus2", "--maxdef", dest="defocus2", type="float",
			help="Maximum defocus 2 (minor axis) in microns", metavar="##")

		self.parser.add_option("-c", "--ampcontrast", dest="ampcontrast", type="float",
			help="Amplitude contrast value values in between 0.0 and 0.5", metavar="##")
		self.ampContrastTypes = ('spider', 'ctffind', 'xmipp')
		self.parser.add_option("-t", "--ampcontrasttype", dest="ampcontrasttype", default="ctffind",
			type="choice", choices=self.ampContrastTypes,
			help="Type of amplitude contrast value", metavar="TYPE")

		self.parser.add_option("-a", "--astigangle", dest="astigangle", type="float",
			help="Astigmatism angle in degrees", metavar="##")
		self.astigAngleTypes = ('????', '???1', '???2')
		self.parser.add_option("-m", "--astigangletype", dest="astigtype", default="????",
			type="choice", choices=self.astigAngleTypes,
			help="Method of measuring astig angle", metavar="METHOD")
	
		return

	#======================
	def checkConflicts(self):
		if not os.path.isfile(self.params['imagefile']):
			apDisplay.printError("File does not exist")
		keys = self.params.keys()
		keys.sort()
		for key in keys:
			if self.params[key] is None:
				apDisplay.printError("%s is not defined"%(key))


		return


if __name__ == '__main__':
	ctfEval = CtfEval()
	ctfEval.start()
	ctfEval.close()



