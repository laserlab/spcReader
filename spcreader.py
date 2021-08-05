# Read bh .spc files.
# This program is a Python library to read SPC files produced by Becker & Hickl
# SPCM software. SPC files are time-tag data from FIFO, containing micro time.
# macro time and detector channels for each individual photon.
#
# Author: Carl, Tim
# 6/21/2021


import warnings
import struct
import os
import numpy as np
import pandas as pd


class readspc:
    
    def __init__(self, arg):
        '''
        input: path to file
        '''
        self.filename = os.fspath(arg)
        #with open(self.filename, 'rb') as fh:
        #    self._fromfile(fh)

    # read lowtime part from buffer
    def _lowtimeread(self,buffer):
        lowtimebuff = buffer[0:3]                   # update lowertime buffer
        channelbuff = buffer[-1] & 0b00011111       # update channel buffer, as an int
        GAPbuff = buffer[-1] & (1 << 5)             # update gap buffer from reading the GAP bit
        return lowtimebuff, channelbuff, GAPbuff
    
    # read hightime part
    def _hightimeread(self,buffer):
        high24_47 = buffer[0:3]   # extract bits for time 24-47
        high48_53 = int.to_bytes(buffer[-1] & 0b00111111,length=1,byteorder='little')   # extract bits for time 48-53
        return high24_47 + high48_53   # set hightime
    
    def topandas(self):
        '''
        Return pandas data frame in the following format:
        | timestamp | chl1 | chl2 | ... | chl20(ref) |
        '''
        times,channels = self.datatoarray()   # extract time stamps and channels arrays
        channeltype = np.unique(channels)    # list all channel types
        columnname = []   # a list of pd.DataFrame column names
        for i in channeltype:   #  get channel names 
            columnname.append('chl '+ str(i))

        # write data into pandas
        data = pd.DataFrame(index=times,columns=columnname)
        for i,ctype in enumerate(channeltype):
            data.loc[:,columnname[i]] = (channels == ctype)

        # remove same elements 
        if data.index.has_duplicates:   #  there's multi clicks at the same time
            if len(channeltype) < 3:    # dula clicks is easy to deal with
                dupindex = data.index[data.index.duplicated()]
                data = data[~data.index.duplicated(keep='first')]
                data.loc[dupindex] = True

            else:  #  multi clicks more than two
                #  extract multi clicks rows
                dupindex = data.index[data.index.duplicated()]
                multi = data.loc[dupindex]
                data = data[~data.index.duplicated(keep='first')]   

                #  loop over rows and channels
                for idx in dupindex:
                    for column in columnname:
                        if multi.loc[:,column].to_numpy().any():
                            data.loc[idx,column] = True
        return data



    def toascii(self, path):
        '''
        This extracts all data: timestamps and channels and write them into a ASCII text file.
        path -> output file name and path
        return: a .txt file
        
        '''
        filestring = ''
        
        with open(self.filename, 'rb') as fh:
            recordbytes = 4

            content = fh.read()
            leng = len(content)
            descriptor = content[:4] # read the descriptor
            bytestring = '{0:08b}'.format(descriptor[-1])

            if bytestring == '11000001':  # check if it's from DPC-230 FIFO
                self.device_type = "DPC-230"
                self.data_type = "pre-processed"
            elif bytestring == '11000101':   # check if it's raw data
                self.data_type = "raw"
                self.device_type = "DPC-230"
                warnings.warn('Data is not pre-processed. It should be used for test only.')
            else:
                raise ValueError('This is not FIFO data from DPC-230, currently cannot be read, translation aborts.')

            self.step = int.from_bytes(descriptor[:3], 'big',signed=False) # Calculate time step in fs

             # loop over all data
            for i in range(1, leng//4):
                data = content[i*4:(i+1)*4] # loop over every 4 bytes
                bitt = data[-1] >> 6

                if bitt == 1:    # it's hightime
                    hightime = self._hightimeread(data)

                elif bitt == 0:   # it's lowtime
                    lowtime, channel, gap = self._lowtimeread(data)
                    if gap != 0:
                        warnings.warn('There is a gap at %0.0f th record, data before this entry is lost.'%(i))
                        
                    timestamp = int.from_bytes(lowtime + hightime, 'little',signed=False)
                    
                    filestring += '{} {:2d}\r\n'.format(timestamp,channel)        # write it into the string
                else:
                    raise ValueError('Data format is wrong at %0.0f th record'%(i))
                    
            with open(path, 'w') as f:
                f.write(filestring)
                
                
                
    def datatoarray(self,chl=0):
        '''
        output data into two arrays;
        
        chl -> the channel number (1-20)
        output -> abs timestamp (fs) for the given channel number, if no input, will return all timestamps contained.
        
        channel dictionary:
        
        channel 1,2  -> CFD 1,2 at TDC1
        channel 11,12 -> CFD 3,4 at TDC2
        channel 3-10  -> LVTTL 1-8 at TDC1
        channel 13-20 -> LVTTL 9-16 at TDC2
        '''
        gaplist = []
        stamplist = []
        channellist = []
                
        with open(self.filename, 'rb') as fh:
            recordbytes = 4
            
            content = fh.read()
            self.leng = len(content)
            descriptor = content[:4]   # read the descriptor
            bytestring = '{0:08b}'.format(descriptor[-1])   # read indicator bits

            # check if it's from DPC-230 FIFO
            if bytestring == '11000001':  
                self.device_type = "DPC-230"
                self.data_type = "pre-processed"
                
            elif bytestring == '11000101':   # check if it's raw data
                self.data_type = "raw"
                self.device_type = "DPC-230"
                warnings.warn('Data is not pre-pocessed. It should be used for test only.')
                
            else:
                raise ValueError('This is not FIFO data from DPC-230, currently cannot be read, translation aborts.')

            self.step = int.from_bytes(descriptor[:3], 'big',signed=False)  # Calculate time unit in fs


            # loop over all data  self.leng//4
            for i in range(1, self.leng//4):
                data = content[i*4:(i+1)*4] # loop over every 4 bytes
                bitt = data[-1] >> 6   # high / low time indicator
                
                if bitt == 1:   # it's high time
                    hightime = self._hightimeread(data)
                # raise ValueError('First data record must be high time part.') 

                elif bitt == 0:
                    lowtime, channel, gap = self._lowtimeread(data)
                    timestamp = int.from_bytes(lowtime + hightime, 'little',signed=False) * self.step
                    
                    if gap != 0:      # there is a gap
                        warnings.warn('There is a gap at %0.0f th record, data before this entry is lost.'%(i))
                        gaplist.append(timestamp)
                        gaplist.append(channel)
                    stamplist.append(timestamp)
                    channellist.append(channel)

                else:
                    raise ValueError('Data format is wrong at %0.0f th record'%(i))
                    
            channelarray = np.array(channellist)
            stamparray = np.array(stamplist) 
            gaparray = np.array(gaplist)
                
        if chl<0 or chl>20:
            raise ValueError('no such channel')
            
        if chl==0:
            output = stamparray, channelarray
                
        else:
            output = stamparray[channelarray==chl]
                
        return output

        
