'''
Created on 9 Oct 2018

@author: thomasgumbricht
'''
import urllib.request
from html.parser import HTMLParser
import os
from sys import exit
from shutil import move
import geoimagine.support.karttur_dt as mj_dt
from geoimagine.kartturmain import Composition, LayerCommon, RasterLayer
#import geoimagine.smap.hdf5_2_geotiff as hdf5_2_geotiff
from geoimagine.smap import hdf5_2_geotiff as hdf5_2_geotiff

class SmapComposition:
    '''
    class for sentinel compositions
    '''
    def __init__(self, compD):  
        for key in compD:
            if '_' in compD[key]:
                exitstr = 'the "%s" parameter can not contain underscore (_): %s ' %(key, compD[key])
                exit(exitstr) 
            setattr(self, key, compD[key])
        if not hasattr(self, 'folder'):
            exitstr = 'All SMAP compositions must contain a folder'
            exit(exitstr)
            
class SmapTile(LayerCommon):
    '''Class for sentinel tiles'''
    def __init__(self, smapid, composition, locusD, datumD, filepath, FN): 
        """The constructor expects an instance of the composition class."""
        LayerCommon.__init__(self)
        self.smapid = smapid
        self.comp = composition
        
        self.locus = locusD['locus']
        self.locuspath = locusD['path']

        self.path = filepath
        self.FN = FN

        self.datum = lambda: None
        for key, value in datumD.items():
            setattr(self.datum, key, value)
        if self.datum.acqdate:
            self._SetDOY()
            self._SetAcqdateDOY()
        self._SetPath()
        self._SetQuery()
        
    def _SetPath(self):
        """Sets the complete path to sentinel tiles"""
        
        self.FP = os.path.join('/Volumes',self.path.volume, self.comp.system, self.comp.source, self.comp.division, self.comp.folder, self.locuspath, self.datum.acqdatestr)
        self.FPN = os.path.join(self.FP,self.FN)
        if ' ' in self.FPN:
            exitstr = 'EXITING smap FPN contains space %s' %(self.FPN)
            exit(exitstr)
            
    def _SetQuery(self):
        self.query = {'smapid':self.smapid, 'tilefilename':self.FN,'source':self.comp.source,'product':self.comp.product,
                 'version':self.comp.version,'acqdate':self.datum.acqdate, 'doy':self.datum.doy, 'folder':self.comp.folder}

class ProcessSmap:
    '''class for SMAP specific processing
    '''   
    def __init__(self, process, session, verbose):
        self.verbose = verbose
        self.process = process   
        self.session = session  

        #Direct to SMAP sub-processes
        if self.process.proc.processid.lower() == 'searchsmapproducts':
            self._SearchSmapProducts()
        elif self.process.proc.processid.lower() == 'smapsearchtodb':
            self._SearchToPostgres()
        elif self.process.proc.processid.lower() == 'downloadsmapdaac':
            self._DownLoadSmapDaac()
        elif self.process.proc.processid.lower() == 'extractsmaphdf':
            self._ExtractSmapHdf()
        elif self.process.proc.processid.lower() == 'checksmap':
            self._CheckSmap()
        else:
            exitstr = 'Exiting, processid %(p)s missing in ProcessSmap' %{'p':self.process.proc.processid}
            exit(exitstr)
 
    def _SearchSmapProducts(self):
        '''IMPORTANT the user credentials must be in a hidden file in users home directory called ".netrc"
        '''
        #Set todays date
        today = mj_dt.Today()
        #Set the serverurl, pand the SMAP roduct and version to search for 
        self.serverurl = self.process.params.serverurl
        self.product = self.process.params.product
        self.version = self.process.params.version
        #check that the version is correctly stated
        if not len(self.version) == 3:
            exit('The smap version must be 3 digits, e.g. "005" or "006"')
        if not self.version.isdigit():
            exit('The smap version must be 3 digits, e.g. "005" or "006"')
        #Set the sensorpath on the server   
        sensorurl = 'SMAP'
        #put the remote search path for the requested dates together
        prodPath ='%s.%s' %(self.product,self.version)
        #create the localpath where the search data (html) will be saved
        localPath = os.path.join('/volumes',self.process.dstpath.volume,'DAAC-SMAP',prodPath)
        if not os.path.exists(localPath):
            os.makedirs(localPath)
        #change to the local directory
        cmd ='cd %s;' %(localPath)
        os.system(cmd)
        #Loop over the dates defined in process
        for datum in self.process.srcperiod.datumD:
            print ('searching',datum)
            #search for the data
            if self.process.srcperiod.datumD[datum]['acqdate'] > today:
                #skip all dates later than today (there can be no images from the future)
                continue
            #convert date to pointed string used on the server
            dateStr = mj_dt.DateToStrPointDate(self.process.srcperiod.datumD[datum]['acqdate'])
            #define the complete url to the SMAP data
            url = os.path.join(self.serverurl,sensorurl,prodPath,dateStr)

            #
            localFPN = os.path.join(localPath,dateStr)
            if os.path.exists(localFPN) and not self.process.overwrite:
                continue
            #Run the wget command including definition of the cookie needed for accessing the server
            cmd ='cd %s;' %(localPath)
            cmd ='%(cmd)s /usr/local/bin/wget -L --load-cookies --spider --no-parent ~/.cookies --save-cookies ~/.cookies %(url)s' %{'cmd':cmd, 'url':url}

            os.system(cmd)
          
    def _SearchToPostgres(self):
        '''Load search holdings to local db
            Does not utilize the layer class but take parameters directly from xml
        '''
        #Set todays date
        today = mj_dt.Today()
        #set the paths
        prodPath ='%s.%s' %(self.process.params.product, self.process.params.version)
        localPath = os.path.join('/Volumes',self.process.srcpath.volume,'DAAC-SMAP',prodPath)
        #Loop over the dates
        for datum in self.process.srcperiod.datumD:
            if self.process.srcperiod.datumD[datum]['acqdate'] > today:
                #skip all dates later than today (there can be no images from the future)
                continue
            #convert date to pointed string used on the server
            dateStr = mj_dt.DateToStrPointDate(self.process.srcperiod.datumD[datum]['acqdate'])
            localFPN = os.path.join(localPath,dateStr)
            #Create a sub-folder called done, when the search results are transferred to the db the html will be moved into the done folder
            tarFPN = os.path.join(localPath,'done',dateStr)
            if not os.path.exists(os.path.split(tarFPN)[0]):
                os.makedirs(os.path.split(tarFPN)[0])
            if os.path.exists(localFPN):    
                self._ReadSMAPhtml(self.session,localFPN,tarFPN,self.process.srcperiod.datumD[datum]['acqdate'])
            else:
                print ('SMAP file missing', localFPN)
    
    def _IniTileDownload(self,statusD):
        '''
        '''
        self.dlL = []
        #create a temp folder to which the download will be directed, only when the download is complete will the data be moved in place
        if not os.path.exists(self.tempFP):
            os.makedirs(self.tempFP)
        #if asscript, the whole downloading will be written as a shell script
        if self.process.params.asscript:
            shFN = 'download_%(prod)s.sh' %{'prod':self.process.params.product}
            shFP = os.path.join(self.tempFP, 'script')
            if not os.path.exists(shFP):
                os.makedirs(shFP)
            self.downloadShFPN = os.path.join(shFP,shFN)
            self.downloadScriptF = open(self.downloadShFPN,'w')
            #cmd = 'mkdir -p %(fp)s;\n' %{'fp':shFP}
            #self.dowloadScriptF.write(cmd)
        #Get the tiles
        tiles = self.session._SelectSmapData(self.process.srcperiod, self.process.params, statusD)
        return tiles
    
    def _CheckSmap(self):
        #Set the expected layers and parameters for filling the db
        queryD = {}
        queryD['product'] = {'val':self.process.params.product, 'op':'=' }
        queryD['retrieve'] = {'val':'Y', 'op':'=' }
        self.paramL = ['source', 'product', 'folder', 'band', 'prefix', 'suffix', 'celltype', 'dataunit', 'scalefac', 'offsetadd', 'cellnull', 'measure', 'retrieve', 'hdffolder', 'hdfgrid']
        self.compL = ['source', 'product', 'folder', 'band', 'prefix', 'suffix', 'celltype', 'dataunit', 'scalefac', 'offsetadd', 'cellnull', 'measure']
        self.extractL = self.session._SelectSMAPTemplate( queryD, self.paramL )

        #from geoimagine.support.modis import DisentangleModisTileName as convFN 
        #First loop over the src folder structure to find all tiles at this position
        #Construct a dummy tile, to get the FP
        smapid = 'smapid' 
        hdfFN = '*.%(hdr)s' % {'hdr':self.process.srcpath.hdrfiletype}

        product = self.process.params.product
        version = self.process.params.version
        source = '%(p)s.%(v)s' %{'p':product,'v':version}

        acqdate = mj_dt.Today()

        tile = (hdfFN, smapid, source, product, version, 'original', acqdate)

        smapTile = self._ConstructDaacTile(tile,self.process.srcpath)
        datepath = os.path.split(smapTile.FPN)[0]
        locuspath = os.path.split(datepath)[0]

        for root, directories, filenames in os.walk(locuspath):
            for filename in filenames:
                
                if filename.endswith(self.process.srcpath.hdrfiletype):

                    queryD = {'smapfilename':filename}
                    paramL = ['smapid', 'smapfilename', 'source', 'product', 'version', 'acqdate']
                    tile = self.session._SelectSingleSMAPDaacTile(queryD,paramL)
                    smapid, smapfilename, source, product, version, acqdate = tile
                    tile = (smapfilename, smapid, source, product, version, 'original', acqdate)
                    smapTile = self._ConstructDaacTile(tile,self.process.srcpath)

                    #Replace is needed for adjusting between SMAP and Karttur default naming conventions
                    source = source.replace('-E','_E')
                    source = source.replace('-S','_S')
                    if os.path.exists(smapTile.FPN): 
                        self.session._InsertSmapData(smapTile.query)
                        statusD = {'smapid': smapid,'column':'downloaded', 'status': 'Y'}
                        self.session._UpdateSmapStatus(statusD)
 
                        #Only tiles found on file are checked, should it be updated
                        self._SearchExtractLayers(acqdate)

                        if self.nrExploded == len(self.extractL):
                            statusD = {'smapid': smapid,'column':'organized', 'status': 'Y'}
                            self.session._UpdateSmapStatus(statusD)
                            statusD = {'smapid': smapid,'column':'exploded', 'status': 'Y'}
                            self.session._UpdateSmapStatus(statusD)
                    else:
                        pass
                        #This should not happen  
                                                 
    def _SearchExtractLayers(self,acqdate):
        '''Search for extracted layers for specific SMAP tile
        '''
        self.nrExploded = 0
        # self.explodeD is not used
        self.explodeD = {}
        for extcomp in self.extractL:
            paramD = dict(zip(self.paramL,extcomp))
            compD = dict(zip(self.compL,extcomp))
                        
            comp = Composition(compD, self.process.system.dstsystem, self.process.system.dstdivision)
            #Set the datum
            
            acqdatestr = mj_dt.DateToStrDate(acqdate)

            datumD = {'acqdatestr': acqdatestr, 'acqdate':acqdate}

            #Construct the locus dictionary
            locusD = {'locus':'global','path':'global'}
            filepath = lambda: None
            filepath.volume = self.process.dstpath.volume; filepath.hdrfiletype = self.process.dstpath.hdrfiletype
            
            #Create a standard raster layer
            layer = RasterLayer(comp, locusD, datumD, filepath)

            if not layer._Exists() or self.process.overwrite:
                self.explodeD[paramD['band']] = {'layer':layer,'params':paramD}
            elif layer._Exists():
                self.session._InsertLayer(layer,self.process.overwrite,self.process.delete)
                self.nrExploded += 1                       
                                                
    def _DownLoadSmapDaac(self):
        '''
        '''
        #create a temp folder to which the download will be directed, only when the download is complete will the data be moved in place
        self.tempFP = os.path.join('/Volumes',self.process.dstpath.volume, 'smap', 'temp')
        statusD = {}
        # TGTODO downloaded must be in xml, defaulted to N and not obligatory
        statusD['downloaded'] = self.process.params.downloaded
        #tiles = self.session._SelectSmapData(self.process.srcperiod, self.process.params, statusD)
        tiles = self._IniTileDownload(statusD)
        
        for tile in tiles:
            self._AddDownload(tile,self.process.dstpath)           
        self._AccessSMAP()
        if self.process.params.asscript:
            self.dowloadScriptF.close()
            
    def _AddDownload(self,tile,sdpath):

        smapTile = self._ConstructDaacTile(tile,sdpath)
        smapfilename, smapid, source, product, version, folder, acqdate = tile
        source = source.replace('-E','_E')
        if os.path.exists(smapTile.FPN): 
            self.session._InsertSMAPtile(smapTile.query)
            statusD = {'smapid': smapid,'column':'downloaded', 'status': 'Y'}
            self.session._UpdateSmapStatus(statusD)
        else:

            if self.process.params.asscript:
                cmd = 'mkdir -p %(FP)s;\n' %{'FP':smapTile.FP}
                self.downloadScriptF.write(cmd)
            datedir = mj_dt.DateToStrPointDate(acqdate)
            localTempFPN = os.path.join(self.tempFP,smapTile.FN)
            self.dlL.append({'query':smapTile.query,'productversion':source,'datedir':datedir,'fn':smapfilename,'dstFPN':smapTile.FPN,'tempFPN':localTempFPN,'smapid':smapid})
  
    def _ConstructDaacTile(self,tile,sdpath):
        '''
        '''
        smapfilename, smapid, source, product, version, folder, acqdate = tile
        #construct the composition
        compD = {'source':source, 'product':product, 'version':version, 'folder':folder, 'system':'smap', 'division':'region'}
        #Invoke the composition
        comp = SmapComposition(compD)
        #Set the datum
        datumD = {'acqdatestr': mj_dt.DateToStrDate(acqdate), 'acqdate':acqdate}
        #Set the filename
        FN = smapfilename
        #Set the locus         
        loc = 'global'
        #Set the locuspath
        locusPath = 'global'
        #Construct the locus dictionary
        locusD = {'locus':loc, 'path':locusPath}
        #Invoke and return a SentinelTile             
        return SmapTile(smapid, comp, locusD, datumD, sdpath, FN)
    
    def _ConstructSmapLayer(self,compD,acqdate,compFormatD):
        '''
        '''
        comp = Composition(compD, self.process.system.dstsystem, self.process.system.dstdivision)
        comp._Update(compFormatD)
        datumD = {'acqdatestr': mj_dt.DateToStrDate(acqdate), 'acqdate':acqdate}
            
        #Set the locus         
        loc = 'global'
            
        #Set the locuspath
        locusPath = 'global'
            
        #Construct the locus dictionary
        locusD = {'locus':loc, 'path':locusPath}
            
        filepath = lambda: None
        filepath.volume = self.process.dstpath.volume; filepath.hdrfiletype = self.process.dstpath.hdr
            
        #Create a standard reaster layer
        bandR = RasterLayer(comp, locusD, datumD, filepath)

        return bandR
                    
    def _ReadSMAPhtml(self,session,FPN,tarFPN,acqdate):
        queryD = self._ParseSmapWgetHTML(FPN)
        session._InsertSmapData(queryD)
        move(FPN,tarFPN)
        
    def _ParseSmapWgetHTML(self, FPN):
        tmpFP = os.path.split(FPN)[0]
        tmpFP = os.path.split(tmpFP)[0]
        tmpFP = os.path.join(tmpFP,'tmpcsv')
        if not os.path.exists(tmpFP):
            os.makedirs(tmpFP)

        FPN = 'file://%(fpn)s' %{'fpn':FPN}
        req = urllib.request.Request(FPN)
        with urllib.request.urlopen(req) as response:
            html = response.read()
        parser = MjHTMLParser()

        parser.queryD = {}
        parser.feed(str(html)) 
        return (parser.queryD)

    def _AccessSMAP(self):   
        '''This is similar to _AccessMODIS
        '''
        serverurl = self.process.params.serverurl
        for tile in self.dlL:
            remotepath = os.path.join(serverurl,'SMAP',tile['productversion'],tile['datedir'])
            url = os.path.join(remotepath,tile['fn']) 

            home = os.path.expanduser("~")
            cookieFPN = os.path.join(home,'.smap_cookies')
            cmd = "curl -n -L -c %(c)s -b %(c)s  %(r)s --output %(l)s;" %{'u':self.process.params.remoteuser, 'c':cookieFPN, 'r':url, 'l':tile['tempFPN']}
            cmd = "%(cmd)s mv %(output)s %(dstFPN)s;" %{'cmd':cmd,'output':tile['tempFPN'], 'dstFPN':tile['dstFPN']}
            if self.process.params.asscript:
                cmdL = cmd.split(';')
                for c in cmdL:
                    if len(c) > 1:
                        writeln = '%(c)s;\n' %{'c':c}
                        self.downloadScriptF.write(writeln)
            else:
                os.system(cmd)
                statusD = {'smapid': tile['smapid'],'column':'downloaded', 'status': 'Y'}
                self.session._UpdateSmapStatus(statusD)            

    def _ExtractSmapHdf(self):
        '''Extract the SMAP hdf file
        '''
        #Set asscript to True, this will create a shell file for downloading all missing tiles, if any
        #self.process.params.asscript = True
        self.tempFP = os.path.join('/Volumes',self.process.srcpath.volume, 'smap', 'temp')
        if self.process.params.asscript:
            
            shFP = os.path.join(self.tempFP, 'script')
            if not os.path.exists(shFP):
                os.makedirs(shFP)
            shFN = 'explode_%(prod)s.sh' %{'prod':self.process.params.product}
            explodeShFPN = os.path.join(shFP,shFN)
            shFN = 'download_%(prod)s.sh' %{'prod':self.process.params.product}
            downloadShFPN = os.path.join(shFP,shFN)
            self.explodeScriptF = open(explodeShFPN,'w')
            self.downloadScriptF = open(downloadShFPN,'w')
            
        #Get the tiles
        statusD = {}
        statusD['downloaded'] = 'Y'
        if not self.process.overwrite and self.process.params.exploded:
            statusD['exploded'] = 'Y'
            
        tiles = self._IniTileDownload(statusD)

        #Search template for layers to extract
        #Get the layers to extract for this product + version
        self.paramL = ['source', 'product', 'folder', 'band', 'prefix', 'suffix', 'celltype', 'dataunit', 'cellnull', 'scalefac', 'measure', 'offsetadd', 'region', 'fileext', 'hdffolder', 'hdfgrid']
        queryD = {'source': '%(p)s.%(v)s' %{'p':self.process.params.product, 'v':self.process.params.version},'retrieve':'Y'}
        self.extractLayerL = self.session._SelectTemplateLayersOnSource(queryD, self.paramL)

        if len(self.extractLayerL) == 0:
            exitstr = 'No layers to exract for smap', queryD
            exit(exitstr)
        missingFlag = False

        for tile in tiles:
            #Construct the smap tile
            smapTile = self._ConstructDaacTile(tile,self.process.srcpath)
            smapfilename, smapid, source, product, version, folder, acqdate = tile
            if not smapTile._Exists():
                warnstr = ('warning the smaptile missing: %s' %(smapTile.FPN))
                print (warnstr)
                self._AddDownload(tile,self.process.srcpath)
                missingFlag = True
                continue 
            nrExploded = self._ExplodeH5(smapTile, acqdate, product)
            print ('    smap.h5, nrexploded',   smapfilename,nrExploded)

            if nrExploded == len(self.extractLayerL):    
                statusD = {'smapid': smapid,'column':'organized', 'status': 'Y'}
                self.session._UpdateSmapStatus(statusD)
                statusD = {'smapid': smapid,'column':'exploded', 'status': 'Y'}
                self.session._UpdateSmapStatus(statusD)
        #Write the missing tiles to the access shell script
        self._AccessSMAP()
        if self.process.params.asscript:
            self.explodeScriptF.close()
            self.downloadScriptF.close()
            printstr = 'To explode tiles you can run the scriptfile %(fpn)s' %{'fpn':explodeShFPN}
            print (printstr)
            if missingFlag:
                printstr = 'To download missing tiles you can run the scriptfile %(fpn)s' %{'fpn':self.downloadShFPN}
                print (printstr)
            
            
            #Loop over the layers to extract in this tile
            '''
            print (extractLayerL)
            for extraclayer in extractLayerL:
                extractD = dict(zip(paramL,extraclayer))

                dstLayer = self._ConstructLayer(extractD,acqdate)

                if dstLayer._Exists():
                    self.session._InsertLayer(dstLayer, self.process.overwrite, self.process.delete)
                    continue
                
                cmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_translate '
                cmd = '%(cmd)s -a_ullr -17367530.45 7314540.11 17367530.45 -7314540.11 ' %{'cmd':cmd}
                cmd = '%(cmd)s -a_srs "+proj=cea +lon_0=0 +lat_ts=30 +ellps=WGS84 +units=m" ' %{'cmd':cmd}
                cmd = '%(cmd)s  -a_nodata -9999  ' %{'cmd':cmd}
                cmd = '%(cmd)s HDF5:"%(hdf)s"://%(folder)s/%(grid)s %(dst)s' %{'cmd':cmd,
                        'hdf':smapTile.FPN,'folder':extractD['hdffolder'],'grid':extractD['hdfgrid'], 
                        'dst':dstLayer.FPN}

                os.system(cmd)
                BREAK HERE
                cmd = '/Library/Frameworks/GDAL.framework/Programs/gdal_translate -a_srs "EPSG:3410" '
                cmd = '%(cmd)s HDF5:"%(hdf)s"://%(folder)s/%(grid)s %(dst)s' %{'cmd':cmd,
                        'hdf':smapTile.FPN,'folder':extractD['hdffolder'],'grid':extractD['hdfgrid'], 
                        'dst':dstLayer.FPN}
                print (cmd)

                #Some products have more than one hdffolder and thus more than one lon/lat pair 
                if product == 'SPL3SMP' and extractD['hdffolder'] == 'Soil_Moisture_Retrieval_Data_PM':
                    queryD = {'hdfgrid':'longitude_pm', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 
                elif product == 'SPL3SMP-E' and extractD['hdffolder'] == 'Soil_Moisture_Retrieval_Data_PM':
                    queryD = {'hdfgrid':'longitude_pm', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 

                else:
                    queryD = {'hdfgrid':'longitude', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 
                lonLayer = self.session._SelectTemplateLayersOnGrid(queryD, paramL)
                print ('lonLayer',lonLayer)
                if lonLayer == None:
                    exitstr = 'No lon/lat data found for SMAP extraction',queryD
                    exit(exitstr)

                lonD = dict(zip(paramL,lonLayer))
                extractD['longrid'] = lonD['hdfgrid']
                if product == 'SPL3SMP' and extractD['hdffolder'] == 'Soil_Moisture_Retrieval_Data_PM':
                    queryD = {'hdfgrid':'latitude_pm', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 
                elif product == 'SPL3SMP-E' and extractD['hdffolder'] == 'Soil_Moisture_Retrieval_Data_PM':
                    queryD = {'hdfgrid':'latitude_pm', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 

                else:
                    queryD = {'hdfgrid':'latitude', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 
                latLayer = self.session._SelectTemplateLayersOnGrid(queryD, paramL)
                latD = dict(zip(paramL,latLayer))

                extractD['latgrid'] = latD['hdfgrid']
                extractD['lonlatfolder'] = latD['hdffolder']
                hdf5_2_geotiff.Retrieve(smapTile.FPN, extractD, dstLayer.FPN)
                
                '''
                   
    def _ExplodeH5(self, smapTile, acqdate, product):
        #  
        nrExploded = 0 
        for extraclayer in self.extractLayerL:
            extractD = dict(zip(self.paramL,extraclayer))

            dstLayer = self._ConstructLayer(extractD,acqdate)

            if dstLayer._Exists():
                if self.process.overwrite:

                    os.remove(dstLayer.FPN)
                else:
                    nrExploded += 1
                    self.session._InsertLayer(dstLayer, self.process.overwrite, self.process.delete)
                    continue
            self._Hdf5_2_geotiff(extractD,product,smapTile,dstLayer)

            if os.path.isfile(dstLayer.FPN):
                nrExploded += 1
                self.session._InsertLayer(dstLayer, self.process.overwrite, self.process.delete)
            
            '''
            #The giving of fixed coordinates is not good
            -17349514.3530680164694786,-7296524.6913595553487539 : 17349514.3530680164694786,7296524.6913595534861088
            cmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_translate '

            cmd = '%(cmd)s -a_ullr -17349514.353 7296524.691 17349514.353 -7296524.691 ' %{'cmd':cmd}


            #SET proj to EASE GRID 2 (epsg:6033)
            cmd = '%(cmd)s -a_srs "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" ' %{'cmd':cmd}
            
            cmd = '%(cmd)s  -a_nodata -9999  ' %{'cmd':cmd}
            cmd = '%(cmd)s HDF5:"%(hdf)s"://%(folder)s/%(grid)s %(dst)s' %{'cmd':cmd,
                    'hdf':smapTile.FPN,'folder':extractD['hdffolder'],'grid':extractD['hdfgrid'], 
                    'dst':dstLayer.FPN}
            print (cmd)

            if self.process.params.asscript:
                cmd = '%(cmd)s;\n' %{'cmd':cmd}
                self.explodeScriptF.write(cmd)
            else:   
                os.system(cmd)
                #register band
                if os.path.isfile(dstLayer.FPN):
                    nrExploded += 1
                    self.session._InsertLayer(dstLayer, self.process.overwrite, self.process.delete)
            '''
        return nrExploded
    
    def _Hdf5_2_geotiff(self,extractD,product,smapTile,dstLayer):
        #Some products have more than one hdffolder and thus more than one lon/lat pair 
        if product == 'SPL3SMP' and extractD['hdffolder'] == 'Soil_Moisture_Retrieval_Data_PM':
            queryD = {'hdfgrid':'longitude_pm', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 
        elif product == 'SPL3SMP-E' and extractD['hdffolder'] == 'Soil_Moisture_Retrieval_Data_PM':
            queryD = {'hdfgrid':'longitude_pm', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 

        else:
            queryD = {'hdfgrid':'longitude', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 
        lonLayer = self.session._SelectTemplateLayersOnGrid(queryD, self.paramL)

        if lonLayer == None:
            exitstr = 'No lon/lat data found for SMAP extraction',queryD
            exit(exitstr)

        lonD = dict(zip(self.paramL,lonLayer))
        extractD['longrid'] = lonD['hdfgrid']
        if product == 'SPL3SMP' and extractD['hdffolder'] == 'Soil_Moisture_Retrieval_Data_PM':
            queryD = {'hdfgrid':'latitude_pm', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 
        elif product == 'SPL3SMP-E' and extractD['hdffolder'] == 'Soil_Moisture_Retrieval_Data_PM':
            queryD = {'hdfgrid':'latitude_pm', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 

        else:
            queryD = {'hdfgrid':'latitude', 'product':extractD['product'], 'source':extractD['source'], 'hdffolder':extractD['hdffolder']} 
        latLayer = self.session._SelectTemplateLayersOnGrid(queryD, self.paramL)
        latD = dict(zip(self.paramL,latLayer))

        extractD['latgrid'] = latD['hdfgrid']
        extractD['lonlatfolder'] = latD['hdffolder']
        print ('retrieving, ',smapTile.FPN,dstLayer.FPN)
        hdf5_2_geotiff.Retrieve(smapTile.FPN, extractD, dstLayer.FPN)
  
    def _ExplodeHDFMODIS(self, hdfFPN, explodeD):
        #  
        nrExploded = 0 
        for band in explodeD:
            tarFPN = explodeD[band]['layer'].FPN
            hdffolder = explodeD[band]['params']['hdffolder']
            hdfgrid = explodeD[band]['params']['hdfgrid']
            #copy the file to memory and extract the hdf straight from memory? 
            cmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_translate '
            cmd = '%(cmd)s HDF4_EOS:EOS_GRID:"%(hdf)s":%(folder)s:%(band)s %(tar)s' %{'cmd':cmd,'hdf':hdfFPN,'folder':hdffolder,'band':hdfgrid, 'tar':tarFPN}

            if self.process.params.asscript:
                cmd = '%(cmd)s;\n' %{'cmd':cmd}
                self.explodeScriptF.write(cmd)
                if self.process.proc.processid.lower() == 'explodemodisregion':
                    self.regionscriptF.write(cmd)
            else:   
                os.system(cmd)
                #register band
                if os.path.isfile(tarFPN):
                    nrExploded += 1
                    self.session._InsertLayer(explodeD[band]['layer'],self.process.overwrite,self.process.delete)
                    #explodeD[band]['layer'].RegisterLayer(self.process.proj.system)
                    #_InsertLayer(self,layer,overwrite,delete)
        return nrExploded
    
    
    def _ConstructLayer(self,extractD,acqdate):
        '''
        '''
        compD = extractD
        comp = Composition(compD, 'smap', 'region')
        datumD = {'acqdatestr': mj_dt.DateToStrDate(acqdate), 'acqdate':acqdate}
            
        #Set the locus         
        loc = extractD['region']
        
        #Set the locuspath
        locusPath = extractD['region']
        
        #Construct the locus dictionary
        locusD = {'locus':loc,  'path':locusPath}
        
        filepath = lambda: None
        filepath.volume = self.process.dstpath.volume; filepath.hdrfiletype = extractD['fileext']
        
        #Create a standard raster layer
        return RasterLayer(comp, locusD, datumD, filepath)
  
        '''
        
        
        srcFPN  = self._GetBandFPN(senTilePath, searchstr,'.jp2')
        
        cmd = ['/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_translate']
        cmd.extend([ '-tr', '%(tr)d' %{'tr':resol} ,' %(tr)d' %{'tr':resol} ])
        cmd.extend(['-ot', celltype, '-a_nodata', '%(cn)d' %{'cn':cellnull} ])
        cmd.extend(['%(src)s' %{'src':srcFPN}, '%(dst)s' %{'dst':bandR.FPN} ]) 
        
        ThisProc = subprocess.check_call(cmd)
        print ('subprocess result', ThisProc)
        '''
    
class MjHTMLParser(HTMLParser):
            
    def handle_starttag(self, tag, attrs):
        # Only parse the 'anchor' tag.
        if tag == "a":
            # Check the list of defined attributes.
            for name, value in attrs:
                # If href is defined, print it.
                if name == "href" and 'SMAP' in value:
                    if value[0:6] == '/SMAP/':
                        source = value.split('/')[2]
                        product,version = source.split('.')
                        self.queryD['source'] = source.replace('_','-')
                        self.queryD['product'] = product.replace('_','-')
                        self.queryD['version'] = version
                    elif value[0:4] == 'SMAP' and os.path.splitext(value)[1] == '.h5':
                        smapfilename = value
                        fnParts = value.split('_')
                        if len(fnParts) == 8 and '_E_' in smapfilename:
                            sensor, level, type, code, enhanced,  acqdatestr, Rcode, vext = value.split('_')
                        elif len(fnParts) == 7:
                            sensor, level, type, code,acqdatestr, Rcode, vext = value.split('_')
                        else:
                            errorstringnotthere
                        acqdate = mj_dt.yyyymmddDate(acqdatestr)
                        self.queryD['smapid'] = os.path.splitext(smapfilename)[0]
                        self.queryD['smapfilename'] = smapfilename
                        self.queryD['acqdate'] = acqdate
                        self.queryD['doy'] = mj_dt.DateToDOY(acqdate)
                    elif value[0:4] == 'SMAP' and os.path.splitext(value)[1] == '.xml':
                        metafilename = value
                    elif value[0:4] == 'SMAP' and os.path.splitext(value)[1] == '.qa':
                        qafilename = value
