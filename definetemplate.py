'''
Created on 26 Oct 2018

@author: thomasgumbricht

Stand alone module for converting output from gdalinfo to the smap.template db table
'''

def ReadSubdatasets(srcL):
    '''ReadSubdatasets
    '''
    #Celltype translator from smap hdf metadata to kartturs GDAL based coding
    srcFPN,EPSG = srcL
    celltypeD = {'16-bit unsigned integer': 'UINT16', '32-bit floating-point': 'FLOAT32', '8-bit unsigned character': 'BYTE'}
    celltypeD['64-bit floating-point'] = 'FLOAT64'
    celltypeD['32-bit unsigned integer'] = 'UINT32'


    #items that are inferred to the db table
    orderL = ['<columns>', 'source', 'product', 'folder', 'prefix', 'suffix', 'fileext', 'region', 'compid', 'hdfgrid', 'hdffolder', 'scalefac', 'offsetadd','title','label','band', 'retrieve', 'measure', 'cellnull', 'celltype', 'dataunit','epsg','</columns>',]
    #dict to hold all layers found in the hdf metadata
    formatD = {}
    #open the textfile that is just a copy and paste text of the output from gdalinfo
    with open(srcFPN) as f:
        content = f.readlines()
    #clean out whitespace and end of line
    content = [x.strip() for x in content]
    for line in content:
        #Identify layers 
        if 'NAME=HDF5' in line:
            print ('line',line)
            fpPartsL = line.split('/')
            hdfgrid = fpPartsL[len(fpPartsL)-1]
            hdffolder = fpPartsL[len(fpPartsL)-2]
            formatD[hdfgrid] = {}
            #Set the EPSG code
            formatD[hdfgrid]['epsg'] = str(EPSG)
            #Set the tags
            formatD[hdfgrid]['<columns>'] = '<values>'
            formatD[hdfgrid]['</columns>'] = '</values>'
            
            #source
            source = fpPartsL[4]
            formatD[hdfgrid]['source'] = source
            
            #product and version
            product, version = source.split('.')
            formatD[hdfgrid]['product'] = product
            formatD[hdfgrid]['version'] = version
            
            #hdfgrid and hdffolder
            formatD[hdfgrid]['hdfgrid'] = hdfgrid
            formatD[hdfgrid]['hdffolder'] = hdffolder
            
            #searchstring for identifying additional metadata further down
            formatD[hdfgrid]['searchcstring'] = '%s_%s' %( hdffolder, hdfgrid)
            
            #band
            band  = hdfgrid.replace('_','-')
            band = band.replace('_retrieved','')
            if band == 'soil-moisture' and 'SPL3SMP' in product:
                if hdffolder == 'Soil_Moisture_Retrieval_Data_AM':
                    band = '%(b)s-am' %{'b':band}
                elif hdffolder == 'Soil_Moisture_Retrieval_Data_PM':
                    band = '%(b)s-pm' %{'b':band}
                else:
                    print (product, hdffolder)
                    KOLLAFELET
            formatD[hdfgrid]['band'] = formatD[hdfgrid]['prefix'] = band
            if len(formatD[hdfgrid]['band']) > 32:
                print ('Band name too long',formatD[hdfgrid]['band'])
                EXITHERE
      
            #title
            title = hdfgrid.replace('_',' ')
            formatD[hdfgrid]['title'] = title
            
            #folder
            folder = hdffolder.lower()
            folder = folder.replace('_retrieval','')
            folder = folder.replace('_data','')
            folder = folder.replace('_polar','')
            folder = folder.replace('_global','')
            formatD[hdfgrid]['folder'] = folder.replace('_','-')
            if len(formatD[hdfgrid]['folder']) > 32:
                print ('folder name too long',formatD[hdfgrid]['folder'])
                EXITHERE
            #suffix 
            formatD[hdfgrid]['suffix'] = version
            
            #fileext
            formatD[hdfgrid]['fileext'] = 'tif'
            
            #compid
            compid = '%(f)s_%(b)s' %{'f':formatD[hdfgrid]['folder'].lower(), 'b':formatD[hdfgrid]['band'].lower()}
            formatD[hdfgrid]['compid'] = compid
            
            #region
            if 'polar' in hdffolder.lower():
                region = 'polar'
            else:
                region = 'global'
            formatD[hdfgrid]['region'] = region
            
            #dataunit
            formatD[hdfgrid]['dataunit'] = 'index'
            
            #scalefac
            formatD[hdfgrid]['scalefac'] = '1'
            
            #offsetadd
            formatD[hdfgrid]['offsetadd'] = '0'
            
            #retrieve defaulted to N
            formatD[hdfgrid]['retrieve'] = 'N'
            
            #cellnull
            formatD[hdfgrid]['cellnull'] = '-1'
            
            #label
            formatD[hdfgrid]['label'] = 'label?'
            
            #measure
            formatD[hdfgrid]['measure'] = 'R'
            
        elif 'DESC=' in line:
            #identify celltype for the description (always follows the NAME field   
            partsL = line.split('(')
            celltype = partsL[1].split(')')[0]
            formatD[hdfgrid]['celltype'] = celltypeD[celltype]
            
    #With all the grids identified, loop again to find the unit, cellnull and long name (label)
    for key in formatD:
        for line in content:
            if formatD[key]['searchcstring'] in line:
                if '__FillValue' in line:
                    if line.split('=')[1] != '?':
                        formatD[key]['cellnull'] = line.split('=')[1]
                elif '_units' in  line:
                    formatD[key]['dataunit'] = line.split('=')[1]
                elif '_long_name' in  line:
                    label = line.split('=')[1]
                    #remove comma
                    label = label.replace(',','')
                    label = label.replace('&apos;s','s')
                    formatD[key]['label'] = label

    #Print the columns with the column tag
    print (",".join(orderL ))
    for key in formatD:
        itemL = []
        for item in orderL:
            itemL.append(formatD[key][item])
        printstr = '            %(v)s' %{'v':"','".join(itemL )}
        printstr = printstr.replace ("<values>',","<values>")
        printstr = printstr.replace (",'</values>","</values>")
        printstr = printstr.replace ("<columns>,","<columns>")
        printstr = printstr.replace (",</columns>","</columns>")
        #Print the values with the values tag
        print (printstr )
        
    print ('Copy and paste the tags to the xml file used for defining smap templates.')
    print ('You might need to edit special characters, and you MUST set the grids you want to retrieve manually.')

    


if __name__ == "__main__":
    #srcFPN = ['/Users/thomasgumbricht/Dropbox/projects/geoimagine/USERS/karttur/SMAP/gdalinfo-excerpts/SPL3FTP.txt',4310]
    #srcFPN = ['/Users/thomasgumbricht/Dropbox/projects/geoimagine/USERS/karttur/SMAP/gdalinfo-excerpts/SPL3FTP-E.txt',4310]
    #srcFPN = ['/Users/thomasgumbricht/Dropbox/projects/geoimagine/USERS/karttur/SMAP/gdalinfo-excerpts/SPL3SMP.txt',4310]
    #srcFPN = ['/Users/thomasgumbricht/Dropbox/projects/geoimagine/USERS/karttur/SMAP/gdalinfo-excerpts/SPL3SMP-E.txt',4310]
    #srcFPN = '/Users/thomasgumbricht/Dropbox/projects/geoimagine/USERS/karttur/SMAP/gdalinfo-excerpts/SPL3SMAP.txt'
    #srcFPN = '/Users/thomasgumbricht/Dropbox/projects/geoimagine/USERS/karttur/SMAP/gdalinfo-excerpts/SPL3SMA.txt'
    srcFPN = ['/Users/thomasgumbricht/Dropbox/projects/geoimagine/USERS/karttur/SMAP/gdalinfo-excerpts/SPL3SMP-E.txt',4310]
    
    ReadSubdatasets(srcFPN)
 
            