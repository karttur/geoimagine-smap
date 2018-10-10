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

class ProcessSmap:
    'class for modis specific processing'   
    def __init__(self, process, session, verbose):
        #self.session = session
        self.verbose = verbose
        self.process = process   
        self.session = session         
        #direct to subprocess
        print ('processid',self.process.proc.processid )

        if self.process.proc.processid.lower() == 'searchsmapproducts':
            self._SearchSmapProducts()
        
        if self.process.proc.processid == 'loadDataPool':
            self._LoadDataPool()
       
    def _SearchSmapProducts(self):
        '''IMPORTANT the user credentials must be in a hidden file in users home directory called ".netrc"
        '''
        today = mj_dt.Today()
        self.serverurl = self.process.params.serverurl
        self.product = self.process.params.product
        self.version = self.process.params.version
        
        if not len(self.version) == 3:
            exit('The modis version must be 3 digits, e.g. "005" or "006"')
        if not self.version.isdigit():
            exit('The modis version must be 3 digits, e.g. "005" or "006"')
            
        sensorurl = 'SMAP'
        
        #put the remote search path for the requested dates together
        
        prodPath ='%s.%s' %(self.product,self.version)
        localPath = os.path.join('/volumes',self.process.dstpath.volume,'DAAC-SMAP',prodPath)

        if not os.path.exists(localPath):
            os.makedirs(localPath)
        cmd ='cd %s;' %(localPath)
        os.system(cmd)
        #dates involved in the search
        #for datum in self.process.srcperiod.datumD:
        #    print (datum, self.process.srcperiod.datumD[datum]['acqdate'])

        for datum in self.process.srcperiod.datumD:
            print ('searching',datum)
            #search the datapool
            if self.process.srcperiod.datumD[datum]['acqdate'] > today:
                continue
            dateStr = mj_dt.DateToStrPointDate(self.process.srcperiod.datumD[datum]['acqdate'])
        
            url = os.path.join(self.serverurl,sensorurl,prodPath,dateStr)
            localFPN = os.path.join(localPath,dateStr)
            #urlStr = '"%s"' %(urlStr)
            if os.path.exists(localFPN) and not self.process.overwrite:
                continue
            print ('url',url)
            print ('localFPN',localFPN)
            
            #cmd ='/usr/local/bin/wget --spider --no-parent  %(url)s' %{'url':url}
            cmd ='cd %s;' %(localPath)
            cmd ='%(cmd)s /usr/local/bin/wget -L --load-cookies --spider --no-parent ~/.cookies --save-cookies ~/.cookies %(url)s' %{'cmd':cmd, 'url':url}

            #print (cmd)

            os.system(cmd)
            
    def _LoadDataPool(self,session):
        '''Load dotapool holdings to local db
            Does not utilize the layer class but take parameters directly from xml
        '''
        prodPath ='%s.%s' %(self.process.params.product, self.process.params.version)
        localPath = path.join('/Volumes',self.process.srcpath.volume,prodPath) 
        for date in self.process.srcperiod.datumL:
            print ('    Loading',self.process.params.product, self.process.params.version, date)
            dateStr = '%(y)s.%(m)s.%(d)s' %{'y':date[0:4],'m':date[4:6],'d':date[6:8]}
            #dateStr = mj_dt.DateToStrPointDate(date)         
            localFPN = path.join(localPath,dateStr)
            print (localFPN)

            tarFPN = path.join(localPath,'done',dateStr)
            if not path.exists(path.split(tarFPN)[0]):
                makedirs(path.split(tarFPN)[0])
            if path.exists(localFPN):    
                #print (self.process.srcperiod.datumD[date]['acqdate'])
                self._ReadMODIShtml(session,localFPN,tarFPN,self.process.srcperiod.datumD[date]['acqdate'])
            else:
                print ('MODIS bulk file missing', localFPN)
                
    def _ReadMODIShtml(self,session,FPN,tarFPN,acqdate):
        tmpFPN,headL = self._ParseModisWgetHTML(FPN)
        session._LoadBulkTiles(self.process.params,acqdate,tmpFPN,headL)
        #move the done file to a subdir called done
        move(FPN,tarFPN)
        
    def _ParseModisWgetHTML(self, FPN):
        headL = ['tileid','tilefilename','source','product','version','acqdate','h','v','hvtile']
        tmpFP = path.split(FPN)[0]
        tmpFP = path.split(tmpFP)[0]
        tmpFP = path.join(tmpFP,'tmpcsv')
        if not path.exists(tmpFP):
            makedirs(tmpFP)
        tmpFPN = path.join(tmpFP,'tilelist.csv')
        FPN = 'file://%(fpn)s' %{'fpn':FPN}
        req = urllib.request.Request(FPN)
        with urllib.request.urlopen(req) as response:
            html = response.read()
        parser = MjHTMLParser()
        parser.SetLists(headL)
        parser.feed(str(html)) 
        WriteCSV(parser.hdfL,tmpFPN)
        return tmpFPN, headL
    

class MjHTMLParser(HTMLParser):
    def SetLists(self,headL):
        self.hdfL = []
        self.xmlL = []   
        self.hdfL.append(headL)
        self.xmlL.append(headL)
    
    def SplitModisFileName(self,value):
        FNparts = value.split('.')
        product = FNparts[0]
        #doy = int(FNparts[1][5:8])
        acqYYYYdoy = FNparts[1][1:8]
        acqdate = mj_dt.yyyydoyDate(acqYYYYdoy)
        #htile = int(FNparts[2][1:3])
        #vtile = int(FNparts[2][4:6])
        version = FNparts[3]
        #prodid = FNparts[4]
        #filetype = FNparts[len(FNparts)-1]
        source = '%(prod)sv%(v)s' %{'prod':product,'v':FNparts[3]}
        tileid = '%(prod)s-%(v)s-%(yyyydoy)s-%(hv)s' %{'prod':product,'v':FNparts[3], 'yyyydoy':FNparts[1][1:8],'hv':FNparts[2] }
        hdfL = [tileid, value, source, product, version, acqdate]
        #D = {'tileid':tileid,'version':version,'tilefilename':value,'source':'MODIS','product':product,'acqdate':acqdate,'doy':doy,'folder':'orignal','htile':htile,'vtile':vtile}
        self.hdfL.append(hdfL)
        
    def handle_starttag(self, tag, attrs):
        # Only parse the 'anchor' tag.
        if tag == "a":
            # Check the list of defined attributes.
            for name, value in attrs:
                # If href is defined, print it.
                if name == "href":
                    ext = path.splitext(value)[1]
                    if ext.lower() == '.hdf':
                        self.SplitModisFileName(value)

def WriteCSV(csvL,tmpFPN):
    import csv
    with open(tmpFPN, 'w') as csvfile:
        wr = csv.writer(csvfile)
        for x,row in enumerate(csvL):
            if x > 0:
                hvtile = row[0].split('-')[3]
                h = int(hvtile[1:3])
                v = int(hvtile[4:6])
                row.extend([h,v,hvtile])
            wr.writerow(row)
        
    
    