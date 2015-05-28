from collections import namedtuple
import os.path
import itertools
import six
import numpy as np
import ptw
import Data_Reduction
import operator
import weakref
import Stanton

class measurement(object):
    """
    Measurement object contains all data relevant to a measurement
    """ 
    
    def __init__(self, filepath, shape, size, height, pressure, LE, measurement_slice=None, point = (0,0), scale = 1, undisturbed = None, Ktresh_hold = 2):
        """
        Initializes a measurement object, keeps track of all measurements and other data
        @param filepath: filepath
        @param shape: roughness shape
        @type shape: string 
        @param size: roughness size (mm)
        @param height: roughness height (mm)
        @param pressure: total pressure (bar)
        @param measurement_slice: slice containing measurement
        """ 
        self.shape = shape
        self.size = size
        self._pressure = 0.28*pressure *100000
        self._filepath = filepath
        self.height = height
        self.slice = measurement_slice
        self._undisturbed = undisturbed
        self.leading_edge_distance = LE
        self._gamma = 1.4
        self._R = 287.05 
        self._M = 7.5 #mach
        self._T0 = 579#775 #kelvin
        self._mu0 = 0.00001827
        self._cp = 1005
        self._C = 120
        self.Tw = 290
        
        self._chord = 1 #meters
        self.point = point #abs px
        self.scale = scale #px/mm

        self._a = np.sqrt(self.gamma*self.R*self.T0) #speed of sound
        
        self.possibleCalculations = (lambda:self.data.data, 
                        lambda:self.data.ml_temp, 
                        lambda:self.data.ml_q,
                        lambda:self.data.ml_relq,
                        lambda:self.data.ml_st,
                        lambda:self.data.ml_K,
                        lambda:self.data.ml_binK)
        self.possibleCalculationNames = ("Raw Data", "Matlab Temp", "Matlab q","Norm q value","Stanton Number","K value","ThreshHolded K")
        self._data = None
        self.reynolds = None
        self.st_lam = None
        self.st_turb = None
        self.Ktresh_hold = Ktresh_hold = 2
        
        
    def calcReML(self):
        m = self
        T = 61.#self.ml_temp*(1+(m.gamma-1)/2*m.M**2)**-1
        P = m.static_pressure
        C=120.
        rho = P/(m.R*T) #density
        T0 = 579.#m.T0
        mu = m.mu0*(T0+C)/(T+C)*(T/T0)**1.5
        v = mu/rho #dynamic viscosity
        Rex = m.V*m.chord/v
        x0 = self.data.cols - self.slice[1][0]
        x1 = self.data.cols - self.slice[1][1] 
        print("=====", self.data.cols, x1, x0, self.point[0], self.leading_edge_distance -(self.data.cols - x1 - self.point[0])/m.scale)
        xarr = np.array(range(x1, x0) ) + 0.5
        return Rex*((xarr)/(m.scale*1000)  )
        # calculate and insert other flow & measurement parameters here!
    @property
    def cp(self):
        return self._cp   
    @property
    def mu0(self):
        return self._mu0
    @property
    def C(self):
        return self._C
    @property 
    def a(self):
        return self._a
    @property
    def chord(self):
        return self._chord
           
    @property
    def gamma(self):
        return self._gamma
    @property
    def pressure(self):
        return self._pressure/100000
    @property
    def R(self):
        return self._R
    @property
    def M(self):
        return self._M
    @property
    def T0(self):
        return self._T0
        #recalculate other flow parameters here!
    @property
    def V(self):
        return self.M*self.a
    @property
    def static_pressure(self):
        return self.pressure_pascal()*(1+(self.gamma-1)/2*self.M**2)**(-self.gamma/(self.gamma-1))
    
    def pressure_bar(self):
        return self._pressure/100000
   
    def pressure_pascal(self):
        return self._pressure
    
    def rho_total_upstream(self):
        return self.pressure_pascal()/(self.R*self.T0)
    
    
    @property
    def offsets(self):
        return (self.data.offsets[1], self.data.offsets[0], self.data.offsets[2]) 
    
    @property
    def filepath(self):
        return self._filepath
    @filepath.setter
    def filepath(self, newpath):
        self._filepath = newpath
    
    @property
    def undisturbed(self):
        return self._undisturbed
    @undisturbed.setter
    def undisturbed(self, undisturbed):
        self._undisturbed = undisturbed
        if self.data is not None:
            self.data.makeUndisturbedData(self._undisturbed)
    
    @property
    def data(self):
        return self._data
    def load(self):
        if os.path.exists(self.filepath):
            self._data = ptw.ptw_file(self.filepath, qfilename = self.qFilename(), measurement=weakref.ref(self))
            if self.data is not None:
                self.reynolds = self.calcReML()
                self.st_lam = Stanton.stanton_laminar(self, self.reynolds)
                self.st_turb = Stanton.stanton_turbulent(self, self.reynolds)
        else:
            print("File not found", self.filepath)
            self._data = None
    def unload(self):
        self._data = None
        Data_Reduction.data_reduction_constant_time.preCalcedPeriods = {}
    def readSlice(self):
        if self.data is not None:
            self.data.readSlice(self.slice, self.undisturbed)
    def qFilename(self):
        fpath, p = os.path.split(self.filepath)
        fname, ext = os.path.splitext(p)
        return os.path.join(fpath, fname + "_Q.dat")
    def saveQ(self):
        if self.data is not None:
            self.data.saveQ()
    def isLoaded(self):
        return self.data is not None

    def to_absolute_pos(self, relative_pos):
        return tuple(map(operator.add,relative_pos, self.offsets))
    
    def to_relative_pos(self, absolute_pos):
        return tuple(map(operator.sub,absolute_pos, self.offsets))
    
    def to_point_pos(self, absolute_pos):
        return tuple(map(operator.sub,absolute_pos, self.point))
        

    def __repr__(self):
        return "measure" + str((self.filepath, str(self.leading_edge_distance)))
    
    
    def __str__(self):
        #LE_str = str()
        return "measure" + str((self.filepath, str(self.leading_edge_distance)))


class simple_measurement(object):
    """
    Measurement object contains all data relevant to a measurement
    """ 
    
    def __init__(self, measurement, data = None, get_lambda = None):
        """
        Initializes a measurement object, keeps track of all measurements and other data
        @param filepath: filepath
        @param shape: roughness shape
        @type shape: string 
        @param size: roughness size (mm)
        @param height: roughness height (mm)
        @param pressure: total pressure (bar)
        @param measurement_slice: slice containing measurement
        """ 
        self.shape = measurement.shape
        self.size = measurement.size
        self._pressure = measurement._pressure *100000
        self._filepath = measurement.filepath
        self.height = measurement.height
        self.leading_edge_distance = measurement.leading_edge_distance
        self._gamma = 1.4
        self._R = 287.05 
        self._M = 7.5 #mach
        self._T0 = 775 #kelvin
        self._mu0 = 0.00001827
        self._C = 120
        self._cp = 1005
        self._chord = 1 #meters
        self.point = measurement.point #abs px
        self.scale = measurement.scale #px/mm
        self._a = np.sqrt(self.gamma*self.R*self.T0) #speed of sound
        self.offsets = measurement.offsets[:-1]
        self.ind = -1
        self.Ktresh_hold = measurement.Ktresh_hold
        
        if data is None:
            data = measurement.data
        
        if get_lambda is not None:
            self._data = get_lambda(data)
        else:
            self._data = data[:,:,-1]

        self.reynolds = measurement.reynolds
        self.st_lam = measurement.st_lam
        self.st_turn = measurement.st_turb
        
        
        # calculate and insert other flow & measurement parameters here!
    

        
        
    
    
    @property
    def cp(self):
        return self._cp       
    @property
    def gamma(self):
        return self._gamma
    @property
    def pressure(self):
        return self._pressure/100000
    @property
    def R(self):
        return self._R
    @property
    def M(self):
        return self._M
    @property
    def T0(self):
        return self._T0
        #recalculate other flow parameters here!
    @property
    def V(self):
        return self.M*self.a
    @property
    def static_pressure(self):
        return self.pressure_pascal()*(1+(self.gamma-1)/2*self.M**2)**(-self.gamma/(self.gamma-1))
    
    def pressure_bar(self):
        return pressure/100000
    
    def pressure_pascal(self):
        return pressure
    
    @property
    def filepath(self):
        return self._filepath
    @filepath.setter
    def filepath(self, newpath):
        self._filepath = newpath
    
    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, newdata):
        self._data = newdata
    
    def to_absolute_pos(self, relative_pos):
        return tuple(map(operator.add,relative_pos, self.offsets))
    
    def to_relative_pos(self, absolute_pos):
        return tuple(map(operator.sub,absolute_pos, self.offsets))
    
    def to_point_pos(self, absolute_pos):
        return tuple(map(operator.sub,absolute_pos, self.point))
        
    
    def __repr__(self):
        return "measure" + str((self.filepath, str(self.leading_edge_distance)))
    
    
    def __str__(self):
        #LE_str = str()
        return "measure" + str((self.filepath, str(self.leading_edge_distance)))



class all_measurements(object):
    """
    List to keep track of all measurements
    """ 
    
    def __init__(self, filepaths = ("", "")):
        """
        Initializes a measurement object, keeps track of all measurements and other data
        @param loc30mm: path to 30 mm LE data
        @param loc60mm: path to 60 mm LE data 
        """ 
        #self.measurements_conn = sqlite3.connect(:memory:)
        self.KeyTypes = namedtuple("KeyTypes", ["shape","size", "height", "pressure", "LE"])
        tmp = (list([]) for _ in range(len(self.KeyTypes._fields)))
        self.keys = self.KeyTypes(*tmp )
        self.measurements = np.empty([0]*len(self.keys), dtype = list) #(4,5,3,3)
        self.filepaths = filepaths
        
    
    def __getitem__(self, key):
        """
        item access, format: ["shape","size", "height", "pressure", "LE"]
        """
        ind = []
        if isinstance(key, str):
            ind.append(next(j for j,v in enumerate(self.keys.shape) if v == key))
        else:
            for i,k in enumerate(key):
                if isinstance(k, slice):
                    start = next(j for j,v in enumerate(self.keys[i]) if v == k.start) if k.start is not None else None
                    stop = next(j for j,v in enumerate(self.keys[i]) if v == k.stop) if k.stop is not None else None
                    step = k.step
                    newk = slice(start, stop, step)
                    ind.append(newk)
                else: 
                    try:
                        print(k, self.keys[i])
                        tst = next(j for j,v in enumerate(self.keys[i]) if v == k)
                    except StopIteration as _:
                        print("Error: trying to open inexisting index", k)
                        raise
                    ind.append(tst)
        return self.measurements[tuple(ind)]
    
    def __iter__(self):
        return self.measurements.__iter__()
    
    def __len__(self):
        return sum(len(x.tolist()) for x in np.nditer(self.measurements, flags=("refs_ok",)))
    
    def _newShape(self, newshape):
        if newshape != self.measurements.shape:
            tlist = self.measurements.tolist()
            
        self.measurements = np.empty(newshape, dtype = list)

        for (sh,si,h,p,le),v in np.ndenumerate( self.measurements ):
            try:
                self.measurements[sh,si,h,p,le] = tlist[sh][si][h][p][le]
            except IndexError:
                self.measurements[sh,si,h,p,le] = []
        
    
    def add_measurement(self, measurement):

        li = [measurement.shape, measurement.size, measurement.height, measurement.pressure/0.28, measurement.leading_edge_distance]
        measurement.filepath = self.get_path(measurement.leading_edge_distance) + measurement.filepath
        ind = [0]*len(li)
        for i,k in enumerate(li):
            curkey = self.keys[i]
            if k not in curkey:
                ind[i] = len(curkey)
                curkey.append(k)
            else:
                ind[i] = next(j for j,v in enumerate(curkey) if v == k)

        size = [max(self.measurements.shape[i], v+1) for i,v in enumerate(ind)]
        self._newShape(size)

        
        ind = tuple(ind)

        if isinstance(self.measurements[ind], list): 
            self.measurements[ind].append(measurement)
        else:
            self.measurements[ind] = [measurement]
 
    def get_measurements(self, shape = slice(None, None, None), 
                         size = slice(None, None, None), 
                         height = slice(None, None, None), 
                         pressure = slice(None, None, None),  
                         LE = slice(None, None, None),                          
                         fname = None):
        
        """
        Item access
        @param shape: shape slice
        @param size: size slice
        @param height: height slice
        @param pressure: pressure slice
        @param LE: leading edge slice
        @param fname: filepath. If fname is given, loads the specific filenames and ignores other params  
        @return: List of measurements fitting criteria
        """
        ret = []
        items = self.__getitem__((shape, size, height, pressure, LE ))
        if fname is not None:
            for idx,mlist in np.ndenumerate( items ):
                try:
                   ret.extend(m for m in mlist if os.path.split(m.filepath)[1] == fname)
                except TypeError:
                    if os.path.split(mlist.filepath)[1] == fname:
                        ret.append(mlist)
        else:
            try:
                return list(itertools.chain(*[mlist for idx,mlist in np.ndenumerate( items ) if len(mlist) > 0]))    
            except TypeError:
                return items
            
        return ret
    
    def get_path(self, leading_edge_size):
        if leading_edge_size == 30:
            return self.filepaths[0]
        elif leading_edge_size == 60:
            return self.filepaths[1]
        else:
            raise ValueError("Incorrect leading edge")