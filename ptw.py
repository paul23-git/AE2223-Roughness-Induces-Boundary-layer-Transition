#! /usr/bin/env python
# Read ptw-file 
import os.path
from scipy import interpolate
import numpy as np
import Data_Reduction



#from PyQt4 import QtGui

class imarray_backend(object):
	def __new__(subtype,data,dtype=None):
		subarr=np.empty(data,dtype=dtype)
		return subarr
	
	def __init__(self, xoff, yoff, toff, lenx, leny, lent):
		self.offsets = (xoff, yoff, toff)
		self.lengths = (lenx, leny, lent)
		self.imarray = np.zeros([lenx, leny, lent])
	
	def __setitem__(self, key, value):
		newkey = (self.__slice_ext_to_int(k, i) for i,k in enumerate(key))
		return self.imarray.__setitem__(tuple(newkey), value)
	
	def __len__(self):
		return len(imarray)
	
	def __slice_ext_to_int(self, s, i):
		if isinstance(s, slice):
			start = s.start - self.offsets[i] if s.start is not None else None
			stop = s.stop - self.offsets[i] if s.start is not None else None
			return slice(start, stop, s.step)
		else:
			return s - self.offsets[i]
	
	def getRaw(self):
		return self.imarray	
	
	
	def __getitem__(self, key):
		newkey = (self.__slice_ext_to_int(k, i) for i,k in enumerate(key))
		return self.imarray.__getitem__(tuple(newkey))
	
	@property
	def shape(self):
		return self.imarray.shape
	


class ptw_file(object):
	class Deferred(object):
		def isIntialized(self):
			return self._initizalized
		def __init__(self, initializer):
			self._initializer = initializer
			self._value = None
			self._initialized = False
		def __call__(self):
			if not self._initialized:
				self._value = self._initializer()
				self._initialized = True
			return self._value
	
	def __init__(self, filename=None, copyFrom=None, slice = None, qfilename = None): #Reads the file header and obtains info on the file
		#maintenance variables
		self._py_temp = self.Deferred(self.python_convert)
		self._ml_temp = self.Deferred(self.matlab_convert_helper)
		self._py_delta_temp =self.Deferred(self.calcDeltaTPy)
		self._ml_delta_temp =self.Deferred(self.calcDeltaTML)
		self._ml_py_diff = self.Deferred(self.calcMLPyDiff)
		self._ml_q = self.Deferred(self.calcQML)
		self._data = None
		self._offsets = (0,0,0)
		self.qfilename = qfilename
		if filename is not None:
			self.loadFile(filename)
		elif copyFrom is not None and slice is not None:
			self.fname = copyFrom.fname
			self.MainHeaderSize = copyFrom.MainHeaderSize
			self.FrameHeaderSize = copyFrom.FrameHeaderSize
			self.nframes = copyFrom.nframes
			
			self.maxlut = copyFrom.maxlut
			self.minlut = copyFrom.minlut
			self.specialscale = copyFrom.specialscale
			self.scaleunit = copyFrom.scaleunit
			self.scalevalue = copyFrom.scalevalue
			self.cols = copyFrom.cols
			self.rows = copyFrom.rows
			self.bitres = copyFrom.bitres
			self.frameperiod = copyFrom.frameperiod
			self.integration = copyFrom.integration
			self.comment = copyFrom.comment
			self.calibration = copyFrom.calibration
			self.FrameSize = copyFrom.FrameSize
			xtup, ytup, timetup = slice
			#self._data = imarray_backend(xtup[0], ytup[0], timetup[0], xtup[1]-xtup[0], ytup[1] - ytup[0], timetup[1] - timetup[0])
			self._offsets = [x[0] + x[1][0] for x in zip(self._offsets, slice)]#(xtup[0], ytup[0], timetup[0])
			self.data = copyFrom.data[xtup[0]:xtup[1], ytup[0]:ytup[1], timetup[0]:timetup[1]]
			self.depth = timetup[1] - timetup[0]
			
	def sliceSelf(self, slice):
		xtup, ytup, timetup = slice
		self.depth = timetup[1] - timetup[0]
		self._offsets = tuple(x[0] + x[1][0] for x in zip(self._offsets, slice))
		self._data = self._data[xtup[0]:xtup[1], ytup[0]:ytup[1], timetup[0]:timetup[1]]
		#self.offsets = [xtup[0], ytup[0], timetup[0]]

	def loadFile(self, filename):
		self.fname = filename 
		fileobj = open(self.fname, mode='rb')
		print("--- Reading Header data ---")
		fileobj.seek(11)
		self.MainHeaderSize = float(np.fromfile(fileobj,'Int32',1))
		self.FrameHeaderSize = float(np.fromfile(fileobj,'Int32',1))
		#print 'MainHeaderSize:',self.MainHeaderSize
		#print 'FrameHeaderSize:',self.FrameHeaderSize

		fileobj.seek(27,0)
		self.nframes = int(np.fromfile(fileobj,'Int32',1))
		#print('nframes:',self.nframes)

		fileobj.seek(245,0);
		self.minlut=float(np.fromfile(fileobj,'Int16',1))
		self.maxlut=float(np.fromfile(fileobj,'Int16',1))
		#print 'minlut:',minlut
		#print 'maxlut:',minlut

		#if(s.m_maxlut==-1)
		#    s.m_maxlut=2^16-1;
		#end; %if
		fileobj.seek(277,0);
		self.specialscale=float(np.fromfile(fileobj,'UInt16',1))
		#print 'specialscale:',specialscale
		self.scaleunit=np.fromfile(fileobj,'c',10)
		#print 'scaleunit:',scaleunit
		self.scalevalue=np.fromfile(fileobj,'Float32',17)
		#print 'scalevalue:', scalevalue
		#if(s.m_specialscale==0)
		#    s.m_unit='dl';                           % [dl T rad]
		#else
		#    s.m_unit=scaleunit;                      % [dl T rad]
		#end; %if

		fileobj.seek(377,0)
		self.cols=int(np.fromfile(fileobj,'UInt16',1)) # Columns
		self.rows=int(np.fromfile(fileobj,'UInt16',1)) # Rows
		#if s.m_rows==0 
		#   s.m_rows=128;
		#end;%if
		#if s.m_cols==0 
		#   s.m_cols=128;
		#end;%if
		self.bitres=float(np.fromfile(fileobj,'UInt16',1)) # bit resolution
		#print 'cols',self.cols
		#print 'rows',self.rows

		fileobj.seek(403,0);
		self.frameperiod = float(np.fromfile(fileobj,'Float32',1))# frame rate
		self.integration =  float(np.fromfile(fileobj,'Float32',1)) # integration time
		#print('frameperiod',self.frameperiod)
		#print ('integration',self.integration)

		fileobj.seek(563,0)
		self.comment=np.fromfile(fileobj,'c',1000)
		#print('comment',comment)

		fileobj.seek(1563,0)
		self.calibration=np.fromfile(fileobj,'c',100) # calibration file name
		#print ('calibration',calibration)

		fileobj.seek(int(self.MainHeaderSize),0) #skip main header
		fileobj.seek(int(self.FrameHeaderSize),1); #skip frame header
		self.FrameSize = self.FrameHeaderSize + self.cols * self.rows * 2;
		
		fileobj.close()
	
	def read(self,framepointer): #Reads the _data from the file
		print("--- Reading frames _data ---")
		fileobj = open(self.fname, mode='rb')
		self._data = np.zeros([self.rows,self.cols,len(framepointer)])
		#self._data = imarray_backend(0,0,0, self.rows, self.cols, len(framepointer))
		self.depth = len(framepointer)
		self._offsets = (self.offsets[0], self.offsets[1], self.offsets[2] + framepointer[0])
		
		for n in range(len(framepointer)):
			fpointer = framepointer[n]
			fileobj.seek(int(self.MainHeaderSize),0) 
			fileobj.seek(int(fpointer * self.FrameSize),1)
			fileobj.seek(int(self.FrameHeaderSize),1)
			data_dum = np.fromfile(fileobj,'UInt16',self.cols * self.rows)
			data_dum = np.reshape(data_dum, (self.rows,self.cols))
			self._data[:,:,n] = data_dum
		fileobj.close()
		return self._data
	@staticmethod
	def matlab_convert(Im, tint, Tcam, eps):
		if(np.round(tint*1.e6) == 400):
		    print("--- Converting (Matlab 400 microsec) ---")
		    elambda=np.array([2498026.51787585,1784.71556672801,6.67171636375203e-12,3215.59874368947,493144.349437419])
		elif(math.round(tint*1.e6) == 200):
			print("--- Converting (Matlab 200 microsec) ---")
			elambda=np.array([436445.608980990,1408.69064523959,0.394428724501193,0.000889315394332730,186268.678717702])
		else:
			raise ValueError("Incorrect time integral")
		
		v1 =elambda[4]/(np.exp(elambda[1]/Tcam)-elambda[2])
		
		Tm=elambda[1]/np.log((eps*elambda[0]/((Im)-elambda[3]-elambda[4]/(np.exp(elambda[1]/Tcam)-elambda[2])))+elambda[2])
		print("--- Converting done ---")
		return Tm
	
	def matlab_convert_helper(self):
		return ptw_file.matlab_convert(self.dataRaw(),self.integration, 22.4+273.15, 0.86)


	def readSlice(self, slice):
		print("--- Reading frames data ---")
		if slice is not None:
			xtup, ytup, timetup = slice
		else:
			xtup = (0, int(self.rows))
			ytup = (0, int(self.cols))
			timetup = (0, int(self.nframes))
		framepointer = range(*timetup)
		fileobj = open(self.fname, mode='rb')
		l = len(framepointer)
		self._data = np.zeros([self.rows,self.cols,l])
		#self._data = imarray_backend(0,0,0, self.rows, self.cols, len(framepointer))
		
		for n in range(l):
			fpointer = framepointer[n]
			fileobj.seek(int(self.MainHeaderSize),0) 
			fileobj.seek(int(fpointer * self.FrameSize),1)
			fileobj.seek(int(self.FrameHeaderSize),1)
			data_dum = np.fromfile(fileobj,'UInt16',self.cols * self.rows)
			data_dum = np.reshape(data_dum, (self.rows,self.cols))
			self._data[:,:,n] = data_dum
		fileobj.close()
		
		
		#slicing
		self.depth = l		
		self._offsets = tuple(x[0] + x[1][0] for x in zip(self._offsets, (xtup, ytup, timetup)))
		self._data = self._data[xtup[0]:xtup[1], ytup[0]:ytup[1], :]

		
		
		return self._data

	def calcDeltaTPy(self):
		print("--- Calculating dT python _data ---")
		t = np.zeros([self.py_temp.shape[0], self.py_temp.shape[1], 1])
		v  = np.concatenate((t, self.py_temp[:,:,1:] - self.py_temp[:,:,:-1]), axis=2 )
		
		print("--- dT done ---")
		return v
	

	def calcDeltaTML(self):
		print("--- Calculating dT matlab _data ---")
		t = np.zeros([self.ml_temp.shape[0], self.ml_temp.shape[1], 1])
		v  = np.concatenate((t, self.ml_temp[:,:,1:] - self.ml_temp[:,:,:-1]), axis=2 )
		
		print("--- dT done ---")
		return v
	
	def calcMLPyDiff(self):
		print("--- Calculating matlab-python difference ---")
		v = self.py_temp - self.ml_temp
		print("--- Diff done ---")
		return v
  
 
	def calcQML(self):
		print("--- Calculating q ---")
		print(self.qfilename)
		if os.path.exists(self.qfilename):
			try:
				#"H:/AE2223/AE2223/3cm_LE/cylinder_r_4_h_2_100bar_run1_Q.dat"
				print(self.qfilename)
				qFileObj = open(self.qfilename, mode = "rb")
				readoff = np.fromfile(qFileObj, dtype = np.int32, count=3)
				readshape = np.fromfile(qFileObj, dtype = np.int32, count=3)
				newoff = np.array(self.offsets, dtype = np.int32)
				newshape = np.array(self.data.shape, dtype = np.int32)
				print(readshape, readoff, newshape, newoff)
				if np.array_equal(readoff, newoff) and np.array_equal(readshape, newshape):
					print("tst-----------=================")
					v = self.loadQ(qFileObj, readshape)
					qFileObj.close()
					return v
			except EOFError as e:
				qFileObj.close()
		return self.calcQML_Backup()
				
			
	def calcQML_Backup(self):
		print ("%i iterations" % int(self.data.shape[0] * self.data.shape[1] * self.data.shape[2] * (self.data.shape[2]+1)/2))
		v = np.zeros(self.ml_delta_temp.shape)
		T = 590
		P = 28 * 100000
		R = 287
		c = 1005
		k = 1.4
		q = 20
	
	    #compute rho
		rho = P/(R*T)  
		for x in range((self.ml_delta_temp.shape[0])):
			print ("calculating column: ", x)
			for y in range((self.ml_delta_temp.shape[1])):
				#print ("calculating point", x, y)
				v_tmp = Data_Reduction.data_reduction_constant_time(self.ml_delta_temp[x,y,:], 
																self.integration,
																rho, c, k)
				v[x,y,:] = v_tmp
		print(v.shape)
		print("--- q done ---")
		return v

	def python_convert(self): #Convert the raw _data to temperatures
		print("--- Converting (python way) ---")
		dl_tab = [5700,6200,6683,7205,8363,9668,11110,12700];
		t_tab = [-5,0,5,10,20,30,40,50];
		tck = interpolate.splrep(dl_tab,t_tab,s=0)
		temp = np.zeros(self.data.shape)
		for n in range(0,self.data.shape[2]):
			temp_dum = interpolate.splev(self._data[:,:,n].flatten(0),tck,der=0)
			temp[:,:,n] = np.reshape(temp_dum,(self._data.shape[0],self._data.shape[1]))
		temp += 273.15
		#self.calcDeltaT()
		print("--- Converting done ---")
		return temp
	
	#def slice(self, xrange, yrange, timrange):
	
	def createSlice(self, slice):
		return ptw_file(copyFrom=self, slice=slice)
	

	def loadQ(self, qFileObj, shape):
		size = np.prod(shape)
		q_dum = np.fromfile(qFileObj, dtype=np.float64, count=size)
		print("after loading", q_dum.shape, q_dum)
		shape = tuple(shape.tolist())
		if size != q_dum.shape[0]:
			raise EOFError("Too file too small " + str(q_dum.shape[0]))
		temp = np.zeros(shape)
		temp = np.reshape(q_dum, shape)
		return temp
		
		
	def saveQ(self):
		filename = self.qfilename
		print(filename, "writing")
		fileobj = open(filename, mode='wb')
		off = np.array(self.offsets)
		
		shape = np.array(self.data.shape)
		print(off.dtype, shape.dtype)
		print(shape)
		off.tofile(fileobj)
		shape.tofile(fileobj)
		self.ml_q.tofile(fileobj)
		fileobj.close()
		
	
	def dataRaw(self):
		return self._data
	@property
	def data(self):
		return self._data
	@data.setter
	def data(self, numpy_array):
		self._data = numpy_array

	@property
	def offsets(self):
		return self._offsets
	@property
	def shape(self):
		return self._data.shape
	
	def _change_data(self, numpy_array):
		self._data = numpy_array
	
	@property
	def py_temp(self):
		if self._data is not None:
			return self._py_temp()
		else:
			return None

	@property
	def ml_temp(self):
		if self._data is not None:
			return self._ml_temp()
		else:
			return None
	
	@property
	def py_delta_temp(self):
		if self._data is not None:
			return self._py_delta_temp()
		else:
			return None

	@property
	def ml_delta_temp(self):
		if self._data is not None:
			return self._ml_delta_temp()
		else:
			return None
	@property
	def ml_py_diff(self):
		if self._data is not None:
			return self._ml_py_diff()
		else:
			return None
  
	@property
	def ml_q(self):
		if self._data is not None:
			return self._ml_q()
		else:
			return None
