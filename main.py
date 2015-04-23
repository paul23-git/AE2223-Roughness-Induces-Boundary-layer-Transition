
from collections import namedtuple
import os.path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import gspread
import Data_Reduction
import Measurements
import ptw

class interface(object):
    def __init__(self, figure, measurement, actual_data = None):
        font = {'family' : 'sans serif',
            'size'   : 12}
        fontsmall = {'family' : 'sans serif',
            'size'   : 10}
        mpl.rcParams["axes.labelsize"] = 10
        mpl.rcParams["xtick.labelsize"] = 10
        mpl.rcParams["ytick.labelsize"] = 10
        mpl.rcParams["legend.fontsize"] = 10
        mpl.rcParams["legend.labelspacing"] = 0.3         
        mpl.rcParams["legend.borderpad"] = 0.3
        
        if (actual_data is not None):
            self.actual_data = actual_data
        else:
            self.actual_data = measurement.data.data
        
        self.display = figure
        self.current_time = 0#measurement.offsets[2]
        self.measurement = measurement
        self.gs = mpl.gridspec.GridSpec(4,4, width_ratios=[50,1,20,10], height_ratios=[30,20,20,20])
        
        
        self._xRange = (0, 0+self.actual_data.shape[0])
        self._yRange = (0, 0+self.actual_data.shape[1])
        self._timeRange = (0, 0+self.actual_data.shape[2])

        
        self.axes = self.display.add_subplot(self.gs[0:2, 0])
        self.axColour = self.display.add_subplot(self.gs[0:2,1])
        self.axText = self.display.add_subplot(self.gs[0,-1], axisbg='lightgoldenrodyellow')
        self.axRadio = self.display.add_subplot(self.gs[1,-1], axisbg='lightgoldenrodyellow')
        self.axDetail = self.display.add_subplot(self.gs[2, :-1])
        self.axColumn = self.display.add_subplot(self.gs[0:2,2])
        self.axTime = self.display.add_subplot(self.gs[3, :-1])
        self.axSaveButton = self.display.add_subplot(self.gs[3,-1], axisbg='lightgoldenrodyellow')
        #self.gs.tight_layout(self.display)
        self.img = self.axes.imshow(self.actual_data[:,:,self.current_time], interpolation='none', 
                                    extent=(0,self.actual_data.shape[1],
                                            self.actual_data.shape[0], 0
                                            ))
        min_temp = np.amin(self.actual_data)
        max_temp = np.amax(self.actual_data)
        self.img.set_clim(min_temp, max_temp)
        #self.maxtime = len(self.actual_data[0,0,:])-1
        self.axes.set_autoscalex_on(False)
        self.axes.set_autoscaley_on(False)
        
        ticksize = np.max((1, int(round((int(np.ceil(max_temp))-int(min_temp))/10.0))))
        self.cbar = self.display.colorbar(self.img, self.axColour, ticks=range(int(min_temp), int(np.ceil(max_temp)),ticksize))
        #print(self.cbar.ax.get_xticklabels(minor=True))
        #self.cbar.ax.set_xticklabels(np.round(self.cbar.ax.get_xticklabels()))
        
        self.axDetail.set_title("Row Temperature", **font)
        self.axDetail.set_ylabel("T*")
        self.axDetail.set_xlabel("horizontal position")
        self.axDetail.set_xlim([self._yRange[0], self._yRange[1]-1])
        self.axDetail.set_ylim([min_temp, max_temp])
        self.axDetail.set_autoscalex_on(False)
        self.axDetail.set_autoscaley_on(False)
        self.axTime.set_title("Time profile point", **font)
        self.axTime.set_ylabel("T*")
        self.axTime.set_xlabel("time")
        self.axTime.set_xlim([self._timeRange[0], self._timeRange[1]-1])
        self.axTime.set_ylim([min_temp, max_temp])
        self.axTime.set_autoscalex_on(False)
        self.axTime.set_autoscaley_on(False)
        self.axColumn.set_title("Column Temperature", **font)
        self.axColumn.set_ylabel("T*")
        self.axColumn.set_xlabel("vertical position")
        self.axColumn.set_xlim([self._xRange[0], self._xRange[1]-1])
        self.axColumn.set_ylim([min_temp, max_temp])
        self.axColumn.set_autoscalex_on(False)
        self.axColumn.set_autoscaley_on(False)
        
        self.detail_pos = (np.nan, np.nan)
        self.hor_line = None
        self.ver_line = None
        self.detail_line = None
        self.time_line = None
        self.detail_ver_line = None
        self.column_line = None
        self.column_ver_line = None
        dat_slice = self.actual_data[self._xRange[0]:self._xRange[1], self._yRange[0]:self._yRange[1], self._timeRange[0]:self._timeRange[1]]
        min_temp = np.amin(self.actual_data[self._xRange[0]:self._xRange[1], self._yRange[0]:self._yRange[1], self._timeRange[0]:self._timeRange[1]])
        max_temp = np.amax(self.actual_data[self._xRange[0]:self._xRange[1], self._yRange[0]:self._yRange[1], self._timeRange[0]:self._timeRange[1]])
        self.time_ver_line, = self.axTime.plot([0, 0], [min_temp, max_temp], '-', color='black',linewidth=1)
        #self.
        
        bmean = np.zeros(self._timeRange[1] - self._timeRange[0])
        
        for t in range(self._timeRange[0], self._timeRange[1]):
            bmean[t-self._timeRange[0]] = np.mean(self.actual_data[self._xRange[0]:self._xRange[1],
                                                                self._yRange[0]:self._yRange[1],
                                                                t])
            
        trange = range(*self._timeRange)
        self.totmean_line, = self.axTime.plot(trange, bmean, color='b', linewidth=1, label="Image mean")

        self.boxmean_line, = self.axTime.plot(trange, bmean, color='r', linewidth=1, label="Box mean")
        self.boxmean_line2, = self.axTime.plot(trange, bmean, color='r', linewidth=1.5)
        self.axTime.legend()
        
        self.mean = np.mean(self.actual_data[:,:,self.current_time])
        self.boxmean = np.mean(self.actual_data[self._xRange[0]:self._xRange[1],self._yRange[0]:self._yRange[1],self.current_time])
        self.rowmean = 0
        self.val = 0
        self.axText.set_autoscalex_on(False)
        self.axText.set_autoscaley_on(False)
        self.axText.set_xlim([0,1])
        self.axText.set_ylim([0,1])
        self.axText.get_xaxis().set_visible(False)
        self.axText.get_yaxis().set_visible(False)
        txt = "(%i, %i, %i) \n" % self.measurement.data.offsets
        txt += "(%.0f, %.0f) \n" % self.detail_pos
        txt += "total mean \n    %.2f \n" % self.mean
        txt += "box mean \n    %.2f \n" % self.boxmean
        txt += "row mean \n    %.2f \n" % self.rowmean
        txt += "Value: %.2f \n" % self.val
        self.Txt = self.axText.text(0.01,1-0.01,txt,verticalalignment='top', **fontsmall)
                
        
        #axcolor = 'lightgoldenrodyellow'
        self.radioLabels = self.measurement.possibleCalculationNames
        self.radioFuncs = self.measurement.possibleCalculations
        self.radio = mpl.widgets.RadioButtons(self.axRadio, self.radioLabels)
        for txt in self.radio.labels:
            txt.set_size(10)
        self.radio.on_clicked(self.ChangeDataSetButton)
        self.saveButton = mpl.widgets.Button(self.axSaveButton, "save\nQ")
        self.saveButton.on_clicked(self._saveQFunc)
        
        self.bid = self.display.canvas.mpl_connect('key_press_event', self.on_press)
        self.cid = self.display.canvas.mpl_connect('button_press_event', self.on_mouse_press)
        
    def ChangeDataSetButton(self, label):
        i = self.radioLabels.index(label)
        self.ChangeDataSet(self.radioFuncs[i]())

    
    def ChangeDataSet(self, newset):
        self.actual_data = newset
        self.rescaleOnSlice(self._xRange, self._yRange, self._timeRange)
        if self.detail_line is not None:
            self._updateDetailLine(self.detail_pos[1], self.current_time)
        
        if self.time_line is not None:
            self._updateTimeLine(self.detail_pos[0],self.detail_pos[1])
        if self.column_line is not None:
            self._updateColumnLine(self.detail_pos[0], self.current_time)
        
        self.img.set_data(self.actual_data[:,:,self.current_time])
        self._updateBoxMeanLine()
        self._updateMeanLine()
        self.display.canvas.draw()
        self._updateLabels()
        
            
        
    def _updateLabels(self):
        self.mean = np.mean(self.actual_data[:,:,self.current_time])
        self.boxmean = np.mean(self.actual_data[self._xRange[0]:self._xRange[1],self._yRange[0]:self._yRange[1],self.current_time])
        if self.detail_line is not None:
            self.rowmean = np.mean(self.actual_data[self.detail_pos[1], : , self.current_time])
            self.val = self.actual_data[self.detail_pos[1], self.detail_pos[0], self.current_time]
        else:
            self.rowmean = 0
            self.val = 0
        txt = "(%.0f, %.0f, %.0f) \n" % self.measurement.data.offsets
        txt += "pos (%.0f, %.0f) \n" % self.detail_pos
        txt += "time %i \n" % int(self.current_time)
        txt += "total mean \n    %.2f \n" % self.mean
        txt += "box mean \n    %.2f \n" % self.boxmean
        txt += "row mean \n    %.2f \n" % self.rowmean
        txt += "Value: %.2f \n" % self.val
        self.Txt.set_text(txt);
            
        
    def update(self, time):
        self._updateTime(time)
    
        
    def rescaleOnSlice(self, xRange, yRange, timeRange):
        self._xRange = xRange
        self._yRange = yRange
        self._timeRange = timeRange
        min_temp = np.amin(self.actual_data[xRange[0]:xRange[1], yRange[0]:yRange[1], timeRange[0]:timeRange[1]])
        max_temp = np.amax(self.actual_data[xRange[0]:xRange[1], yRange[0]:yRange[1], timeRange[0]:timeRange[1]])
        self._reScaleTemperature(min_temp, max_temp)
        self._updateBoxMeanLine()
        self._updateLabels()
        
        
    def _updatePos(self, xnew, ynew):
        
        self._updateTimeLine(xnew, ynew)
        if ynew != self.detail_pos[1]:
            self._updateColumnverLine(ynew)
            self._updateDetailLine(ynew, self.current_time)
        
        if xnew != self.detail_pos[0]:
            self._updateDetailverLine(xnew)
            self._updateColumnLine(xnew, self.current_time)
        
        
        
        
        if xnew != self.detail_pos[0]:
            self._updateVerLine(xnew)
        if ynew != self.detail_pos[1]:
            self._updateHorLine(ynew)

            
        self.detail_pos = (xnew, ynew)
        self._updateLabels()
        self.display.canvas.draw()

    def _updateTime(self, tnew):
        tnew = np.clip(tnew, self._timeRange[0], self._timeRange[1]-1)
        if tnew != self.current_time:
            self.img.set_data(self.actual_data[:,:,tnew])

            self._updateTimeVerLine(tnew)
            if self.detail_line is not None:
                self._updateDetailLine(self.detail_pos[1], tnew)
            if self.column_line is not None:
                self._updateColumnLine(self.detail_pos[0], tnew)
            self.current_time = tnew
            self._updateLabels()
            self.display.canvas.draw()
            
    
    def _updateBoxMeanLine(self):
        bmean = np.zeros(len(self.actual_data[0,0,:]))
        for t in range(len(self.actual_data[0,0,:])):
            bmean[t] = np.mean(self.actual_data[self._xRange[0]:self._xRange[1],
                            self._yRange[0]:self._yRange[1], t])
        self.boxmean_line.set_ydata(bmean)
        self.boxmean_line.set_xdata(range(*self._timeRange))
        self.boxmean_line2.set_ydata(bmean[self._timeRange[0]:self._timeRange[1]])
        self.boxmean_line2.set_xdata(range(*self._timeRange))
    def _updateMeanLine(self):
        tmean = np.zeros(len(self.actual_data[0,0,:]))
        for t in range(len(self.actual_data[0,0,:])):
            tmean[t] = np.mean(self.actual_data[:,:, t])
        self.totmean_line.set_xdata(range(*self._timeRange ))
        self.totmean_line.set_ydata(tmean)
    def _updateTimeVerLine(self, tnew):
        if self.time_ver_line is not None: 
            self.time_ver_line.set_xdata([tnew, tnew])
        else:
            ymin, ymax = self.axTime.get_ylim()
            self.time_ver_line, = self.axTime.plot([tnew, tnew], [ymin, ymax], '-', color='black',linewidth=1)
    def _updateTimeLine(self, xnew, ynew):
        xrange = range(self._timeRange[0], self._timeRange[1])
        if self.time_line is not None:
            self.time_line.set_xdata(xrange)
            self.time_line.set_ydata(self.actual_data[ynew, xnew , :])
        else:
            yrange = self.actual_data[ynew, xnew , :]
            
            self.time_line, = self.axTime.plot(xrange, yrange, color='g', linewidth=1.5, label="Point")
            self.axTime.legend()
            
    def _updateDetailverLine(self, xnew):
        if self.detail_ver_line is not None:
            self.detail_ver_line.set_xdata([xnew, xnew])
        else:
            ymin, ymax = self.axDetail.get_ylim()
            self.detail_ver_line, = self.axDetail.plot([xnew, xnew], [ymin, ymax], '-', color='black',linewidth=1)
    def _updateDetailLine(self, ynew, tnew):
        xrange = range(*self._yRange )
        if self.detail_line is not None:
            self.detail_line.set_xdata(xrange)
            self.detail_line.set_ydata(self.actual_data[ynew, : , tnew])
        else:
            yrange = self.actual_data[ynew, : , tnew]
            self.detail_line, = self.axDetail.plot(xrange, yrange)
            
    def _updateColumnverLine(self, ynew):
        if self.column_ver_line is not None:
            self.column_ver_line.set_xdata([ynew, ynew])
        else:
            ymin, ymax = self.axColumn.get_ylim()
            self.column_ver_line, = self.axColumn.plot([ynew, ynew], [ymin, ymax], '-', color='black',linewidth=1)
    def _updateColumnLine(self, xnew, tnew):
        yrange = range(*self._xRange )
        if self.column_line is not None:
            self.column_line.set_xdata(yrange)
            self.column_line.set_ydata(self.actual_data[:, xnew , tnew])
        else:
            xrange = self.actual_data[:, xnew , tnew]
            self.column_line, = self.axColumn.plot(yrange, xrange)
    def _updateVerLine(self, xnew):
        if self.ver_line is not None:
            self.ver_line.set_xdata([xnew, xnew])
        else:
            ymin, ymax = self.axes.get_ylim()
            self.ver_line, = self.axes.plot([xnew, xnew], [ymin, ymax], '-', color='black',linewidth=1)
    def _updateHorLine(self, ynew):
        if self.hor_line is not None:
            self.hor_line.set_ydata([ynew, ynew])
        else:
            xmin, xmax = self.axes.get_xlim()
            self.hor_line, = self.axes.plot([xmin, xmax], [ynew, ynew], '-', color='black',linewidth=1)
    def _reScaleTemperature(self, min_temp, max_temp):
        self.img.set_clim(min_temp, max_temp)
        self.axDetail.set_ylim([min_temp, max_temp])
        self.axColumn.set_ylim([min_temp, max_temp])
        self.axTime.set_ylim([min_temp, max_temp])
        if self.detail_ver_line is not None:
            self.detail_ver_line.set_ydata([min_temp, max_temp])
        if self.column_ver_line is not None:
            self.column_ver_line.set_ydata([min_temp, max_temp])
        
        self.time_ver_line.set_ydata([min_temp, max_temp])
        ticksize = np.max((1, int(round((int(np.ceil(max_temp))-int(min_temp))/10.0))))
        self.cbar.set_ticks(range(int(min_temp), int(np.ceil(max_temp)),ticksize))
        self._updateLabels()

    def _saveQFunc(self, event):
        self.measurement.saveQ()
        
    
    def on_press(self, event):
        t = None
        if (event.key == "right"):
            t = self.current_time + 1
        if (event.key == "left"):
            t = self.current_time - 1 
        
        if (t != None):
            self.update(t)
    
    def on_mouse_press_main(self, event):
        self._updatePos(round(event.xdata), round(event.ydata))
    
    def on_mouse_press_detail(self, event):
        self._updatePos(round(event.xdata), self.detail_pos[1])

    def on_mouse_press_column(self, event):
        self._updatePos(self.detail_pos[0], round(event.xdata))
        
    def on_mouse_press_time(self, event):
        self._updateTime(round(event.xdata))
        
        
    def createSubFigure(self, axes):
        ftemp = plt.figure()
        ax = ftemp.add_subplot("111")
        if axes == self.axes:
            ax.imshow(self.actual_data[:,:,self.current_time], interpolation='none')
        elif axes == self.axDetail and self.detail_line is not None:
            xrange = range(*self._yRange)
            yrange = self.actual_data[self.detail_pos[1], : , self.current_time]
            ax.plot(xrange, yrange)
        elif axes == self.axColumn and self.column_line is not None:
            xrange = range(*self._xRange)
            yrange = self.actual_data[:, self.detail_pos[0] , self.current_time]
            ax.plot(xrange,yrange)
        elif axes == self.axTime:
            trange = range(*self._timeRange)
            if self.time_line is not None:
                yrange = self.actual_data[self.detail_pos[1], self.detail_pos[0] , :]
                ax.plot(trange,yrange, color='g', linewidth=1.5, label="Point")
            bmean = np.zeros(self._timeRange[1] - self._timeRange[0])
            
            for t in range(self._timeRange[0], self._timeRange[1]):
                bmean[t-self._timeRange[0]] = np.mean(self.actual_data[self._xRange[0]:self._xRange[1],
                                                                self._yRange[0]:self._yRange[1],
                                                                t])
            
            tmean = np.zeros(self.actual_data.shape[2])
            for t in range(0, self.actual_data.shape[2]):
                bmean[t] = np.mean(self.actual_data[:,:,t])
            
            ax.plot(trange, tmean, color='b', linewidth=1, label="Image mean")
    
            ax.plot(trange, bmean, color='r', linewidth=1, label="Box mean")
            ax.plot(trange, bmean, color='r', linewidth=1.5)
            ax.legend()
        ftemp.axes.append(axes)
        plt.show()
        
    def on_mouse_press(self, event):
        if event.button == 1:
            if event.inaxes == self.axes: self.on_mouse_press_main(event)
            elif event.inaxes == self.axDetail and self.detail_line is not None: self.on_mouse_press_detail(event)
            elif event.inaxes == self.axTime: self.on_mouse_press_time(event)
            elif event.inaxes == self.axColumn and self.column_line is not None: self.on_mouse_press_column(event)
        else:
            if event.inaxes == self.axes or \
                event.inaxes == self.axDetail or \
                event.inaxes == self.axTime or \
                event.inaxes == self.axColumn: 
                self.createSubFigure(event.inaxes)

        

def convNameDataDetail(name):
    tok = name.split("_")
    j = next(j for j,v in enumerate(tok) if v == "r" or v == "l" or v[-3:] == "bar" or v[:3] == "run")
    shape = "_".join(tok[:j])
    i = 1
    l = len(tok)
    DocsData = namedtuple("DocsData", ["shape","height", "size", "pressure"])
    size = 0;
    height = 0;
    pressure = 0;
    while i<l:
        curtok = tok[i]
        if curtok == "h" and i < (l-1):
            height = int(tok[i+1])
            i+=1
        elif (curtok == "r" or curtok == "l") and i < (l-1):
            size = int(tok[i+1])
            i+=1
        elif curtok[-3:] == "bar":
            pressure = int(curtok[:-3])
        i+=1
    dat = DocsData(shape=shape, size=size, height=height, pressure=pressure) #Measurements.measurement(fname, shape, height, size, pressure, LE)
    return dat


def convNameData(fname, LE):
    k = fname.rfind(".")
    if k > 0:
        name = fname[:k]
    else:
        name = fname
    #unique typos
    s = "half_spherer"
    if name.startswith(s):
        name = "half_sphere_" + name[len(s)-1:]
    
    if name:
        dat = convNameDataDetail(name)
        return Measurements.measurement(fname, shape = dat.shape, height = dat.height, size = dat.size, pressure = dat.pressure, LE = LE)
    raise ValueError("Badly formatted name %s" % name)

    

def loadAllMeasurementsGoogleDocs(allmeasurements, login, sheetkey, Verbose = False):
    print(" --- Loading measurement from google docs --- ")
    sh = login.open_by_key(sheetkey)
    worksheets = sh.worksheets()
    for ws in worksheets:
        if Verbose:
            print("Sheet", ws.title)
        try:
            LE = int(ws.title)
        except ValueError as e:
            print("Undefined sheet", ws.title)
        else:
            listoflists = ws.get_all_values()
            for row in listoflists[1:]:
                if row and row[0]:
                    name = row[0]
                    try:
                        m = convNameData(name, LE)
                    except ValueError as e:
                        print ("Badly formatted filename", fname)
                    else:
                        allmeasurements.add_measurement(m)
                        if len(row) >= 7:
                            try:
                                slice = ([int(v) for v in row[1:3]], [int(v) for v in row[3:5]], [int(v) for v in row[5:7]])
                                m.slice = slice
                            except ValueError as E:
                                pass
                        if Verbose:
                            print("Loading measurement", m)


def findMeasurementDataFromFilename(all_measurements, fname, leading_edge = slice(None, None, None)):
    k = fname.rfind(".")
    if k > 0:
        name = fname[:k]
    else:
        name = fname
    if name:
        try:
            dat = convNameDataDetail(name)
        except ValueError as e:
            return None
        else:
            return all_measurements.get_measurements(shape = dat.shape, size = dat.size, height = dat.height, pressure = dat.pressure, LE = leading_edge, fname = fname)
            
    return None               
                
                
    #Measurements.measurements_list.add(Measurements.measurement())


def main_loadgoogle(allMeasurements):
    gc = gspread.login("D07TAS", "hypersonicroughness")
        
    
    

    loadAllMeasurementsGoogleDocs(allMeasurements, gc, "1Xw_EXTmFHbKhSj4OGKRff0ClR-_QSO_V_YLUTaOD_GM")
    print("---------------------")

def main_load_data(test_measurements):
    for m in test_measurements:
        m.load()
        m.readSlice()
def main_show_measurements(test_measurements):
    datalist = []
    displist = []
    for m in test_measurements:
        fmain = plt.figure()
        fmain.suptitle(m.filepath)
        disp = interface(fmain, m, m.data.data)
        displist.append(disp)
    plt.show()
    return (datalist, displist)

def main_calculate_and_save_q(test_measurements):
    for m in test_measurements:
        print("-----",m,"-----")
        m.saveQ()
def main_calculate_and_save_q_all_memory_efficient(test_measurements):
    print(len(test_measurements))
    for m in test_measurements:
        print("-----",m,"-----")
        m.load()
        m.readSlice()
        m.data.ml_q
        print("continuing")
        m.saveQ()
        print("saved")
        m.unload()
        print("unload")
    
def main():



    
    filename = "cylinder_r_2_h_1_60bar_run1.ptw"
    allMeasurements = Measurements.all_measurements(("C:\Users\Roeland\Documents\GitHub\AE2223-Roughness-Induces-Boundary-layer-Transition/AE2223/AE2223/3cm_LE/","C:\Users\Roeland\Documents\GitHub\AE2223-Roughness-Induces-Boundary-layer-Transition/AE2223/AE2223/6cm_LE/"))
    main_loadgoogle(allMeasurements)
    #m = convNameData(filename, 30)
    #m.slice = ((90, 180), (120, 300), (29, 48))
    #allMeasurements.add_measurement(m)
    print(len(allMeasurements))
    print("loaded google docs")
    #test_measurements = findMeasurementDataFromFilename(allMeasurements, filename)
    test_measurements = allMeasurements.get_measurements(shape = "square", pressure = 100, height =  2, size = 2, LE=30) 
    """
    Item access
    @param shape: shape slice
    @param size: size slice
    @param height: height slice
    @param pressure: pressure slice
    @param LE: leading edge slice (mm)
    @param fname: filepath. If fname is given, loads the specific filenames and ignores other params  
    @return: List of measurements fitting criteria
    """
    print(test_measurements)
    

    main_load_data(test_measurements)
    #main_calculate_and_save_q_all_memory_efficient(test_measurements)
    main_show_measurements(test_measurements)
    

    print("end")
    
if (__name__ == "__main__"):
    main()
