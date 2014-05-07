
## python
import os
import re
import time
## appion
from appionlib import apDisplay
from appionlib import apFile
from appionlib import spyder

"""
A large collection of SPIDER functions for that are common to other SPIDER functions

I try to keep the trend
image file:
	*****img.spi
image stack file:
	*****stack.spi
doc/keep/reject file:
	*****doc.spi
file with some data:
	*****data.spi

that way its easy to tell what type of file it is
"""

#===============================
def stackToSpiderStack(stack,spiderstack,apix,boxsize,lp=0,hp=0,bin=1,numpart=0):
	"""
	convert an input stack (i.e. imagic) and write it out to a spider-formatted stack
	"""
	emancmd  = "proc2d "
	if not os.path.isfile(stack):
		apDisplay.printError("stackfile does not exist: "+stack)
	emancmd += stack+" "


	emancmd += "apix="+str(apix)+" "
	if lp > 0:
		emancmd += "lp="+str(lp)+" "
	if hp > 0:
		emancmd += "hp="+str(hp)+" "
	if bin > 1:
		clipboxsize = boxsize*bin
		emancmd += "shrink="+str(bin)+" "
		emancmd += "clip="+str(clipboxsize)+","+str(clipboxsize)+" "
	if numpart > 0:
		emancmd += "last="+str(numpart-1)+" "
	emancmd += "spiderswap edgenorm"
	starttime = time.time()
	apDisplay.printColor("Running spider stack conversion this can take a while", "cyan")
	from appionlib import apEMAN
	apEMAN.executeEmanCmd(emancmd, verbose=True)
	apDisplay.printColor("finished eman in "+apDisplay.timeString(time.time()-starttime), "cyan")
	return

#===============================
def spiderOutputLine(int1, int2, float1, float2, float3, float4, float5, float6=1.0):
	line = "%04d" % int1
	line += " %1d" % int2
	line += " "+apDisplay.leftPadString("%3.6f" % float1, n=11)
	line += " "+apDisplay.leftPadString("%3.6f" % float2, n=11)
	line += " "+apDisplay.leftPadString("%3.6f" % float3, n=11)
	line += " "+apDisplay.leftPadString("%3.6f" % float4, n=11)
	line += " "+apDisplay.leftPadString("%3.6f" % float5, n=11)
	line += " "+apDisplay.leftPadString("%3.6f" % float6, n=11)
	line += "\n"
	return line

#===============================
def spiderOutLine(num, floatlist):
	line = "%04d" % num
	line += " %1d" % len(floatlist)
	for fnum in floatlist:
		line += " "+apDisplay.leftPadString("%3.6f" % fnum, n=11)
	line += "\n"
	return line

#===============================
def spiderInLine(line):
	sline = line.strip()
	if sline[0] == ";":
		return None
	bits = sline.split()
	rownum = int(bits[0])
	numfloats = int(bits[1])
	floatlist = []
	for i in range(numfloats):
		floatlist.append(float(bits[i+2]))
	spidict = {
		'row': rownum,
		'count': numfloats,
		'floatlist': floatlist,
	}
	return spidict

#===============================
def addParticleToStack(partnum, partfile, stackfile, dataext=".spi"):
	mySpider = spyder.SpiderSession(dataext=dataext, logo=False)
	mySpider.toSpiderQuiet("CP", 
		spyder.fileFilter(partfile), #particle file
		spyder.fileFilter(stackfile)+"@%06d"%(partnum), #stack file
	)
	mySpider.close()
	return

#===============================
def averageStack(stackfile, numpart, avgfile, varfile, dataext=".spi"):
	mySpider = spyder.SpiderSession(dataext=dataext, logo=True, log=False)
	mySpider.toSpider("AS R", 
		spyder.fileFilter(stackfile)+"@******", #stack file
		"1-%6d"%(numpart), #num of particles
		"A", #use all particles
		spyder.fileFilter(avgfile), #average file
		spyder.fileFilter(varfile), #variance file
	)
	mySpider.close()
	return

#===============================
def createMask(maskfile, maskrad, boxsize, dataext=".spi"):
	"""
	We should use imagefun.filled_circle() instead
	"""
	apDisplay.printMsg("Creating mask with diameter %.1f and boxsize %d"%(maskrad*2.0,boxsize))
	mySpider = spyder.SpiderSession(dataext=dataext, logo=False)
	mySpider.toSpiderQuiet("MO", 
		spyder.fileFilter(maskfile), 
		"%d,%d" % (boxsize, boxsize), 
		"C", 
		str(maskrad),
	)
	mySpider.close()
	if not os.path.isfile(maskfile) or apFile.fileSize(maskfile) < 2:
		apDisplay.printError("Failed to create mask file")
	return

#===============================
def createBoxMask(maskfile,boxsize,maskx,masky,falloff,imask=None,dataext=".spi"):
	"""
	creates a rectangular mask with soft falloff
	"""
	apDisplay.printMsg("Creating %i x %i box mask"%(boxsize,boxsize))
	mySpi = spyder.SpiderSession(dataext=dataext, logo=False, nproc=1, log=False)
	# create blank image for mask
	mySpi.toSpiderQuiet("BL","_1","%i,%i"%(boxsize,boxsize),"N","1")
	# mask it in X
	mySpi.toSpiderQuiet("MA X", "_1", "_2",
		"%i"%maskx, "C", "E", "0",
		"%i,%i"%(boxsize/2,boxsize/2),
		"%.2f"%falloff)
	# inner mask in X
	if imask is not None:
		mySpi.toSpiderQuiet("MA X","_2","_3",
			"%i,%i"%(boxsize/2,imask),
			"C","E","0",
			"%i,%i"%(boxsize/2,boxsize/2),
			"%.2f"%(falloff/4))
		mySpi.toSpiderQuiet("CP","_3","_2")
	# mask in Y
	mySpi.toSpiderQuiet("MA Y","_2",
		spyder.fileFilter(maskfile),
		"%i"%masky, "C", "E", "0",
		"%i,%i"%(boxsize/2,boxsize/2),
		"%.2f"%falloff)
	
	mySpi.close()

	if not os.path.isfile(maskfile) or apFile.fileSize(maskfile) < 2:
		apDisplay.printError("Failed to create mask file")
	return

#===============================
def intListToString(strintlist):
	### convert to ints and sort
	intlist = []
	for item in strintlist:
		intitem = int(item)
		intlist.append(intitem)
	intlist.sort()

	### make list of factors
	intstr = ""
	lastitem = None
	for item in intlist:
		if lastitem == item-1:
			lastend = "-"+str(lastitem) 
			end = intstr[-len(lastend):]
			if end == lastend:
				intstr = re.sub(lastend, "-"+str(item), intstr)
			else:
				intstr += "-"+str(item)
		else:
			intstr += ","+str(item)
		lastitem = item
	intstr = intstr[1:]
	intkey = re.sub(",", "_", intstr)
	return intstr, intkey

