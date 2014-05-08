#!/usr/bin/env python

#pythonlib
import os
#appion
from pyami import mrc
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
		imgdata = {
			'filename': os.path.abspath(self.params['imagefile']),
			'image': mrc.read(self.params['imagefile']),
		}
		ctfdata = {
			'volts': self.params['kv']*1e3,
			'cs': self.params['cs'],
			'apix': self.params['apix'],
			'defocus1': self.params['defocus1']*1e-6,
			'defocus2': self.params['defocus2']*1e-6,
			'angle_astigmatism': self.params['astigangle'],
			'amplitude_contrast': self.params['ampcontrast'],
		}
		a = ctfdisplay.CtfDisplay()
		a.debug = self.params['debug']
		ctfdisplay.ctftools.debug = self.params['debug']
		ctfdisplaydict = a.CTFpowerspec(imgdata, ctfdata, None, None, True)
		if ctfdisplaydict is None:
			raise

		ctfdata['confidence_30_10'] = ctfdisplaydict['conf3010']
		ctfdata['confidence_5_peak'] = ctfdisplaydict['conf5peak']
		ctfdata['overfocus_conf_30_10'] = ctfdisplaydict['overconf3010']
		ctfdata['overfocus_conf_5_peak'] = ctfdisplaydict['overconf5peak']
		ctfdata['resolution_80_percent'] = ctfdisplaydict['res80']
		ctfdata['resolution_50_percent'] = ctfdisplaydict['res50']
		### override the confidence
		ctfdata['confidence'] = max(ctfdisplaydict['conf5peak'], ctfdisplaydict['conf3010'])

		return ctfdata

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
		#self.ampContrastTypes = ('spider', 'ctffind', 'xmipp')
		#self.parser.add_option("-t", "--ampcontrasttype", dest="ampcontrasttype", default="ctffind",
		#	type="choice", choices=self.ampContrastTypes,
		#	help="Type of amplitude contrast value", metavar="TYPE")

		self.parser.add_option("-a", "--astigangle", dest="astigangle", type="float",
			help="Astigmatism angle in degrees", metavar="##")
		#self.astigAngleTypes = ('????', '???1', '???2')
		#self.parser.add_option("-m", "--astigangletype", dest="astigtype", default="????",
		#	type="choice", choices=self.astigAngleTypes,
		#	help="Method of measuring astig angle", metavar="METHOD")
	
		self.parser.add_option("-d", "--debug", dest="debug", default=False,
			action="store_true", help="Turn on debugging features")

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

		if not os.path.isfile(self.params['imagefile']):
			apDisplay.printError("Could not read micrograph for processing")

		if self.params['kv'] > 400.0 or self.params['kv'] < 60:
			apDisplay.printError("atypical high tension value %.1f kiloVolts"%(volts))	

		if self.params['cs'] > 7.0 or self.params['cs'] < 0.4:
			apDisplay.printError("atypical C_s value %.1f mm"%(cs))

		if self.params['apix'] > 20.0 or self.params['apix'] < 0.1:
			apDisplay.printError("atypical pixel size value %.1f Angstroms"%(pixelsize))

		if self.params['defocus1'] > 15.0 or self.params['defocus1'] < 0.1:
			apDisplay.printError("atypical defocus #1 value %.1f microns (underfocus is positve)"%(focus1))

		if self.params['defocus2'] > 15.0 or self.params['defocus2'] < 0.1:
			apDisplay.printError("atypical defocus #2 value %.1f microns (underfocus is positve)"%(focus1))

		if abs(self.params['astigangle']) > 0.01 and abs(self.params['astigangle']) < 3.0:
			apDisplay.printWarning("atypical angle astigmatism value %.3f"%(angle))

		if self.params['ampcontrast'] < 0.0 or self.params['ampcontrast'] > 0.5:
			apDisplay.printError("atypical amplitude contrast value %.3f"%(ampcon))

		return

if __name__ == '__main__':
	ctfEval = CtfEval()
	ctfEval.start()
	ctfEval.close()


