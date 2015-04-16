import mysamp,numpy
myclient=mysamp.demo()
x = numpy.arange(0, 100)
y=x**2
#Stop and start TOPCAT connection to SAMP hub
#Broadcast table to TOPCAT

myclient['t0'] = {'x':x, 'y':y }
