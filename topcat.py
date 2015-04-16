	
if True:
    # The atexit module is just to try and close off zombie clients and hubs.
    # This could be done with a class deconstructor (if that's what they're called in python), but I've been told they're unreliable.
     
    import atpy, atexit, numpy, os, sampy
     
    class sampCommunications:
     
        def __init__(self):
            # Before we start, let's kill off any zombies
            try: samp.kill()
            except: pass
           
            # sampy seems to fall over sometimes if 'localhost' isn't specified, even though it shouldn't
            self.hub = sampy.SAMPHubServer(addr='localhost')
           
            self.hub.start()
           
            self.metadata = {
                                'samp.name' : 'ATPy module',
                                'samp.icon.url' : 'http://dribbble.com/system/users/2302/avatars/thumb/Panda2.png?1281017958', # Gratuitous cute panda
                                'samp.description.text' : 'ATPy Module',
                                'client.version' : '0.9.4'
                            }
           
            self.client = sampy.SAMPIntegratedClient(metadata=self.metadata, addr='localhost')
            self.client.connect()
           
            # Bind interaction functions - we will register that we want to listen for table.highlight.row (the SAMP highlight protocol), and all the typical SAMP-stuff
            # We could include this to listen for other sorts of subscriptions like pointAt(), etc.
            self.client.bindReceiveNotification('table.highlight.row', self.highlightRow)
            self.client.bindReceiveNotification('samp.app.*', self.receiveNotification)
            self.client.bindReceiveCall('samp.app.*', self.receiveCall)
            self.client.bindReceiveResponse('*', self.receiveResponse)
           
            # Try to avoid zombie hubs and clients
            atexit.register(self.kill)
           
           
        # Try to die quietly
        def kill(self):
            if self.client.isConnected(): self.client.disconnect()
            if self.hub._is_running: self.hub.stop()
           
     
        # Generic response protocols
        def receiveNotification(self, privateKey, senderID, mType, params, extra):
            print '[SAMP] Notification ', privateKey, senderID, mType, params, extra
               
        def receiveResponse(self, privateKey, senderID, msgID, response):
            print '[SAMP] Response ', privateKey, senderID, msgID, response
           
        def receiveCall(self, privateKey, senderID, msgID, mType, params, extra):
            print '[SAMP] Call ', privateKey, senderID, msgID, mType, params, extra
            self.client.ereply(msgID, sampy.SAMP_STATUS_OK, result = {'txt' : 'printed' })
           
        # Let's see if topcat is running - this will return the first version of topcat running, if there are multiple/zombies
        def isTopcatRunning(self):
       
            neighbours = self.client.getRegisteredClients()
           
            for neighbour in neighbours:
                metadata = self.client.getMetadata(neighbour)
               
                try:
                    if (metadata['samp.name'] == 'topcat'):
                        self.topcat = neighbour
                        return True
                   
                except KeyError:
                    continue
                       
            self.topcat = None
            return False
           
        # Broadcast a table file to TOPCAT
        def broadcastTable(self, table):
       
            metadata = {
                        'samp.mtype' : 'table.load.votable',
                        'samp.params' : {
                                         'name' : table,
                                         'table-id' : table,
                                         'url' : 'file://' + os.getcwd() + '/' + table
                                        }
                       }
                       
            if self.isTopcatRunning(): return self.client.notify(self.topcat, metadata)
            else: return False
       
        def highlightRow(self, privateKey, senderID, mType, params, extra):
            print '[SAMP] Highlighted row', privateKey, senderID, mType, params, extra
           
            if self.isTopcatRunning() and (senderID == self.topcat):
               
                try:
                    filename, row = [params['url'], int(params['row'])]
                except KeyError:
                    print '[SAMP] Highlighted row was missing vital information'
                else:
                    print 'TOPCAT tells us that row %s of file %s was highlighted!' % (row, filename)
           
            else:
                # Interactions with Aladin?
                pass
           
     
    # In the code that I use, my spectra is analysed by a function, and before it starts it initiates a sampCommunications() class
    # This gives TOPCAT a second to connect to the hub
    samp = sampCommunications()
     
    # Some data is generated from my program. Let's use some mock data.
    x = numpy.arange(0, 100)
    y = x**2
     
    # Put this data into a VO table and save it
    table = atpy.Table()
     
    table.add_column('x', x)
    table.add_column('y(x) = x^2', y)
     
    table.write('test.xml', verbose=False)
     
     
    # Now broadcast the table to TOPCAT
    samp.broadcastTable('test.xml')

