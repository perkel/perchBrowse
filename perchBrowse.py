#!/usr/bin/env python
"""
Version 0.0.1. PerchBrowse. To browse perch torque data

David Perkel 23 January 2016
"""
try:
    import os
    import math
    import sys
    import wx
    import matplotlib
    matplotlib.use('WXAgg')
    from matplotlib.figure import Figure
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_wxagg import \
        FigureCanvasWxAgg as FigCanvas
        #        NavigationToolbar2WxAgg as NavigationToolbar
    from matplotlib.ticker import AutoMinorLocator
    from matplotlib.widgets import SpanSelector

    import numpy as np
    from numpy import mean, sqrt, square
    #from matplotlib import *
#    from findEvents import findEvents

except ImportError, target:
    print 'Software configuration error: %s.' % str(target)
    print 'Please email %s for assistance.' % "perkel@uw.edu"
    raw_input('Press enter to close this window.')
    sys.exit(1)

def rms(vals):
    temp = (np.array(vals, dtype = np.int32)**2).mean()
    if temp < 0.0:
        print "in rms: temp = ", temp
    return np.sqrt(temp)

def align_time(time, ChunkSize, OverlapFraction):
    result = time[::int(ChunkSize * OverlapFraction)]
    return result
    

def rms_new(data, ChunkSize, OverlapFraction):
    list_length = len(data)
    result = []
    start_point = 0
    
    while start_point < list_length:
        chunk = data[start_point:int(start_point + ChunkSize)]
        chunk_rms = sqrt(mean(square(chunk)))
        result.append(chunk_rms)
        start_point += int(ChunkSize - (ChunkSize * OverlapFraction))
    return result
    
def myMax(array):
    # Takes an array of values and returns an array consisting of the max value
    # and the index of that maximum value
    maxVal = array[0]
    maxIdx = 0;
    for i in range(0,len(array)):
        if array[i] > maxVal:
            maxVal = array[i]
            maxIdx = i
    return [maxVal, maxIdx]

## Global Parameters
NFFT = 512       # the length of the windowing segments
noverlap = 0
fileList = [] # stores full path now (BMW 10/6/15)
debug = True
initialPath = "~/Dropbox/Perching Bird project/Data"


class GraphFrame(wx.Frame):
    # The main frame of the application
    title = 'Perch Data Viewer'
    
    def __init__(self):
        # width appears to be defined by control panel
        wx.Frame.__init__(self, None, -1, self.title, size=wx.Size(1000, 600))
        
        self.paused = False
        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()
        self.thisFile = 0
        self.currFileName = ''
        self.nFiles = 0
        self.chunkDur = 0.1 # chunk duration for segmentation (in s)
        self.fs = 0.0 # sampling rate in Hz
        self.dt = 0.0 # reciprocal of fs
        self.data = [] # waveform data from wav file
        self.times = None
        self.spec = []
        self.FFResults = []
        self.FFVals = []
        self.currFolder = os.getcwd()
        self.defaultChunkSize = 10
        self.defaultOverlapFraction = 0.5
        self.defaultMinOffDuration = 0.26
        
    def create_menu(self):
        self.menubar = wx.MenuBar()
        
        menu_file = wx.Menu()
        m_About= menu_file.Append(wx.ID_ABOUT, "&About"," Information about this program")
        self.Bind(wx.EVT_MENU, self.OnAbout, m_About)
        m_open = menu_file.Append(-1, "&Open wav file\tCtrl-O", "Open wav file")
        self.Bind(wx.EVT_MENU, self.on_file_open, m_open)
        menu_file.AppendSeparator()
        m_save = menu_file.Append(-1, "&Save plot\tCtrl-S", "Save plot to file")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_save)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
                
        self.menubar.Append(menu_file, "&File")
        
        menu_analysis = wx.Menu()
        m_Browse = menu_analysis.Append(-1, "&Browse files\tCtrl-B", "Browse spectrograms of wave files")
        self.Bind(wx.EVT_MENU, self.on_browse, m_Browse)
                
        m_Segment = menu_analysis.Append(-1, "&Segment file", "Segment file into syllables")
        self.Bind(wx.EVT_MENU, self.on_segment, m_Segment)
       
        self.menubar.Append(menu_analysis, "&Analysis")

        menu_view = wx.Menu()
        m_Inspect = menu_view.Append(-1, "&Inspect widgets", "Inspect widgets on the panel")
        self.Bind(wx.EVT_MENU, self.on_Inspect, m_Inspect)
        self.menubar.Append(menu_view, "&View")
        self.SetMenuBar(self.menubar)

    def create_main_panel(self):
        self.panel = wx.Panel(self)

        self.init_plot()
        self.canvas = FigCanvas(self.panel, -1, self.fig)
        self.scroll_range = 400
        self.canvas.SetScrollbar(wx.HORIZONTAL, 0, 5,
                                 self.scroll_range)   
        self.canvas.Bind(wx.EVT_SCROLLWIN, self.OnScrollEvt)

        nSizeNumBox = 45    # n pix across for number input text field
        self.hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        self.open_button = wx.Button(self.panel, -1, "Open", style=wx.ALIGN_BOTTOM)
        self.Bind(wx.EVT_BUTTON, self.on_file_open, self.open_button)
        self.browse_button = wx.Button(self.panel, -1, "Browse", style=wx.ALIGN_BOTTOM)
        self.Bind(wx.EVT_BUTTON, self.on_browse, self.browse_button)
        self.next_button = wx.Button(self.panel, -1, "Next>>", style=wx.ALIGN_BOTTOM)
        self.Bind(wx.EVT_BUTTON, self.on_next_file, self.next_button)
        self.prev_button = wx.Button(self.panel, -1, "<<Previous", style=wx.ALIGN_BOTTOM)
        self.Bind(wx.EVT_BUTTON, self.on_prev_file, self.prev_button)
        self.segment_button = wx.Button(self.panel, -1, "Segment", style=wx.ALIGN_BOTTOM)
        self.Bind(wx.EVT_BUTTON, self.on_segment, self.segment_button)
        self.mark_presence_button = wx.Button(self.panel, -1, "Mark bird's presence", style=wx.ALIGN_BOTTOM)
        self.Bind(wx.EVT_BUTTON, self.on_mark_presence, self.mark_presence_button)
        self.save_data_button = wx.Button(self.panel, -1, "Save Data as txt", style=wx.ALIGN_BOTTOM)
        self.Bind(wx.EVT_BUTTON, self.on_save_data, self.save_data_button)
        self.quit_button = wx.Button(self.panel, -1, "Quit", style=wx.ALIGN_BOTTOM)
        self.Bind(wx.EVT_BUTTON, self.on_quit, self.quit_button)
        

        self.file_num = wx.TextCtrl(self.panel, -1, '', size=(nSizeNumBox,-1))
        self.num_files = wx.TextCtrl(self.panel, -1, '', size=(nSizeNumBox,-1))
        self.jump_button = wx.Button(self.panel, -1, "Jump")
        self.Bind(wx.EVT_BUTTON, self.on_jump, self.jump_button)
        self.hbox3.AddSpacer((40,10))
        self.show_plots = wx.CheckBox(self.panel, -1, 'Plot', size=(nSizeNumBox,-1))
        self.auto_segment = wx.CheckBox(self.panel, -1, 'Auto segment', size=(nSizeNumBox*2,-1))
        #self.auto_calcFF = wx.CheckBox(self.panel, -1, 'Auto CalcFF', size=(nSizeNumBox*2,-1))

        self.ChunkSize = wx.TextCtrl(self.panel, -1, '10', size=(nSizeNumBox,-1))
        self.OverlapFraction = wx.TextCtrl(self.panel, -1, '0.5', size=(nSizeNumBox,-1))
        self.noise_threshold = wx.TextCtrl(self.panel, -1, '0.1', size=(nSizeNumBox,-1))
        self.MinOffDuration = wx.TextCtrl(self.panel, -1, '0.26', size=(nSizeNumBox,-1))
         
        self.panel.Bind(wx.EVT_KEY_DOWN, self.onKeyPress)
        self.panel.SetFocus()
    
        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)      
        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
         # holds controls for file_num, junp_button, num_files
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, flag=wx.LEFT | wx.TOP | wx.GROW)        
        self.vbox.Add(self.hbox3, 0, flag=wx.ALIGN_LEFT | wx.TOP)
        self.vbox.Add(self.hbox1, 0, flag=wx.ALIGN_LEFT | wx.TOP)
        self.hbox1.Add(self.open_button, 0, flag=wx.ALIGN_LEFT | wx.BOTTOM)
        self.hbox1.Add(self.browse_button, 0, flag=wx.ALIGN_LEFT | wx.BOTTOM)
        self.hbox1.Add(self.prev_button, 0, flag=wx.ALIGN_RIGHT)
        self.hbox1.Add(self.next_button, 0, flag=wx.ALIGN_RIGHT)
        self.hbox1.Add(self.segment_button, 0, flag=wx.ALIGN_RIGHT)
        self.hbox1.Add(self.mark_presence_button, 0, flag=wx.ALIGN_RIGHT)
        self.hbox1.Add(self.save_data_button, 0, flag=wx.ALIGN_RIGHT)
        self.hbox1.Add(self.quit_button, 0, flag=wx.ALIGN_RIGHT)
        

        self.hbox3.Add(wx.StaticText(self.panel, -1, 'File num (1st is 0): '), 0)
        self.hbox3.Add(self.file_num, 0)
        self.hbox3.AddSpacer((40,10))
        self.hbox3.Add(self.jump_button, 0)
        self.hbox3.AddSpacer((40,10))
        self.hbox3.Add(wx.StaticText(self.panel, -1, '   Num files: '), 0)
        self.hbox3.Add(self.num_files, 0)
        self.hbox3.AddSpacer((40,10))
        self.hbox3.Add(self.show_plots, 0)
        self.hbox3.Add(self.auto_segment, 0)
        self.hbox3.AddSpacer((40,10))
        self.hbox3.Add(wx.StaticText(self.panel, -1, '  RMS Chunk Size:  '), 0)
        self.hbox3.Add(self.ChunkSize, 0)
        self.hbox3.Add(wx.StaticText(self.panel, -1, '  RMS Overlap %:  '), 0)
        self.hbox3.Add(self.OverlapFraction, 0)
        self.hbox3.Add(wx.StaticText(self.panel, -1, '  Noise Threshold:  '), 0)
        self.hbox3.Add(self.noise_threshold, 0)
        self.hbox3.Add(wx.StaticText(self.panel, -1, '  Minimum Off Duration:  '), 0)
        self.hbox3.Add(self.MinOffDuration, 0)


        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)
        self.vbox.Add(self.hbox2, 0, flag=wx.ALIGN_LEFT | wx.TOP)

    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()
    
    def clear_axes(self):
        self.axes.clear()
        self.axes.set_title('Perch angle', size=12)
        self.axes.xaxis.set_minor_locator(AutoMinorLocator())
        self.axes.set_ylabel("Angle (deg)")
        self.axes.tick_params(which='minor',axis='x', color='w')
    
    def clear_axes_two(self):
        self.axes2.clear()
        self.axes2.set_xlabel("Time (s)")
        self.axes2.set_ylabel("RMS (deg)")
        self.axes2.set_title('RMS', size=12)

    def init_plot(self):
        self.dpi = 100
        self.fig = Figure((6.0, 6.0), dpi=self.dpi, facecolor='white')
        self.axes = self.fig.add_subplot(2, 1,1)
        self.axes.set_axis_bgcolor('black')
        self.axes.set_title('Perch angle', size=12)
        self.axes.xaxis.set_minor_locator(AutoMinorLocator())
        self.axes.set_ylabel("Angle (deg)")
        self.axes.tick_params(which='minor',axis='x', color='w')
        
        self.axes2 = self.fig.add_subplot(2, 1, 2)
        self.axes2.set_xlabel("Time (s)")
        self.axes2.set_ylabel("RMS (deg)")
        self.axes2.set_axis_bgcolor('black')
        self.axes2.set_title('RMS', size=12)
        self.axes2.xaxis.set_minor_locator(AutoMinorLocator())
        self.axes2.tick_params(which='minor',axis='x', color='w')
        plt.setp(self.axes.get_xticklabels(), fontsize=12)
        plt.setp(self.axes.get_yticklabels(), fontsize=12)
        plt.setp(self.axes2.get_xticklabels(), fontsize=12)
        plt.setp(self.axes2.get_yticklabels(), fontsize=12)

    def updateStatusBar(self):
        #textString = "Directory: " + os.getcwd() + "; " + "File: " + self.currFileName
        textString = "File: " + self.currFileName
        self.statusbar.SetStatusText(textString)

    def on_file_open(self, event): # open file
        file_choices = "NPY (*.npy)|*.npy"

        dlg = wx.FileDialog(
            self, 
            message="Open npy file...",
            defaultDir=self.currFolder,
            defaultFile="*.npy",
            wildcard=file_choices,
            style=wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.currFileName = dlg.GetPath()
            self.currFolder = os.path.dirname(self.currFileName)
            self.open_file(self.currFileName)
            self.nFiles = 0 # this refers to number you can browse thru
            self.file_num.SetValue('')
            self.num_files.SetValue('')


    def onselect(self, xmin, xmax):
        print "xmin, xmax: ", xmin, xmax
        
        #Create an array consisting of RMS values of data (so this function will also update Chunk Size and Overlap Fraction)
        rms_data = np.array(rms_new(self.data, self.defaultChunkSize, self.defaultOverlapFraction))
        #Create corresponding time array
        rms_time = np.array(align_time(self.times, self.defaultChunkSize, self.defaultOverlapFraction))
        indmin, indmax = np.searchsorted(rms_time, (xmin, xmax))
        indmax = min(len(rms_time) - 1, indmax)
        print "indmin, indmax: ", indmin, indmax
    
        thisx = rms_time[indmin:indmax]
        #thisy = self.axes2.get_ydata()
        thisy = rms_data[indmin:indmax]
        print "thisy min, max: ", thisy.min(), thisy.max()
        self.axes2.set_xlim(xmin, xmax)
        self.axes2.set_ylim(thisy.min(), thisy.max())
        self.canvas.draw()

    def plot(self):
        if debug:
            print "in plot. len data is ", len(self.data)
        self.clear_axes()
        self.axes.plot(self.times, self.data, 'r')
        
        self.clear_axes_two()
        self.axes2.plot(self.times, np.absolute(self.data), 'b.')
        
        self.time_hist_min = self.data.min()
        self.time_hist_max = self.data.max()
        self.canvas.draw()
        # set useblit True on gtkagg for enhanced performance
        self.span = SpanSelector(self.axes, self.onselect, 'horizontal', useblit=True,
                    rectprops=dict(alpha=0.5, facecolor='red'))



    def open_file(self, path):
        tempData = np.load(path) # Load data from file
        self.times = tempData[0]
        self.data = tempData[1]
        self.dt = self.times[1] - self.times[0]
        self.fs = 1.0/self.dt
        #if debug:
        #    print "data loaded"
        
        # make histogram of all points (to look for likely baseline)
        #fig = plt.figure()
        #ax = fig.add_subplot(111)

        numBins = 360
        hRange = (-180.5, 179.5)
        hv = np.histogram(self.data, numBins, hRange)
        pkIndex = np.argmax(hv[0])
        self.data -= np.arange(-180, 180, 1.0)[pkIndex]

        # Plot waveform
        self.plot()
        self.updateStatusBar()
        self.syllable_times = []
        if self.auto_segment.GetValue():
            self.on_segment(0)
        else:
            print "not calling on_segment"

    def on_next_file(self, event):
        if self.nFiles < 1:
            self.WarnDlg('You need to click "Browse" before you can click "Next>>"')
            return
        self.thisFile += 1
        if self.thisFile >= self.nFiles:
            self.thisFile = 0
            print '\a' # beep
        self.currFileName = fileList[self.thisFile]
        self.open_file(self.currFileName)
        self.file_num.SetValue(str(self.thisFile))

    def on_prev_file(self, event):
        if self.nFiles < 1:
            self.WarnDlg('You need to click "Browse" before you can click "<<Previous"')
            return
        if self.thisFile == 0:
            print "self.nFiles = ", self.nFiles
            print "len fileList = ", len(fileList)
            self.thisFile = self.nFiles - 1
            print '\a' # beep
        else:
            self.thisFile -= 1
        self.currFileName = fileList[self.thisFile]
        self.open_file(self.currFileName)
        self.file_num.SetValue(str(self.thisFile))
             
    def on_jump(self, event):
        if self.nFiles < 1:
            self.WarnDlg('You need to click "Browse" before you can click "Jump"')
            return
        self.thisFile = int(self.file_num.GetValue())
        if self.thisFile >= self.nFiles:
            self.thisFile = self.nFiles-1
        elif self.thisFile < 0:
            self.thisFile = 0
        self.open_file(fileList[self.thisFile])
        self.file_num.SetValue(str(self.thisFile))

    def on_segment(self, event):

        ChunkSize = float(self.ChunkSize.GetValue())
        OverlapFraction = float(self.OverlapFraction.GetValue())
        self.clear_axes_two()
        self.axes2.plot(align_time(self.times, ChunkSize, OverlapFraction), rms_new(self.data, ChunkSize, OverlapFraction), 'b.')
        self.canvas.draw()
        
    def OnScrollEvt(self, event):
        window_size_1 = 50
        axes_start = 0 + event.GetPosition()
        axes_end = axes_start + window_size_1 
        
        self.clear_axes()
        #self.axes.autoscale(enable = False, axis = 'y')
        self.axes.plot(self.times[axes_start:axes_end], self.data[axes_start:axes_end], 'r')
        

        self.clear_axes_two()
        #self.axes2.autoscale(enable = False, axis = 'y')
        self.axes2.plot(self.times[axes_start:axes_end], np.absolute(self.data[axes_start:axes_end]), 'b')
        
        self.canvas.draw()


    def getFileList(self):
        # load directory list
        global fileList

        fileList = []
        dlg = wx.DirDialog(
            self, message="Choose a folder", defaultPath=initialPath, 
            )
        if dlg.ShowModal() != wx.ID_OK:
            return 0 # as if zero files
        dirname = dlg.GetPath()
        self.currFolder = dirname

        for file in os.listdir(dirname): # arbitrary order
            if file.endswith(".npy"):
                fileList.append(os.path.join(dirname,file))

        fileList.sort()
        nf = len(fileList)
        self.num_files.SetValue(str(nf))
        print "found ", nf, " npy files"
        return nf
            
    def on_browse(self, event):
        self.nFiles = self.getFileList()
        #print "on_browse. nf = ", self.nFiles
        if self.nFiles > 0:
            # open first file
            self.thisFile = 0
            self.currFileName = fileList[self.thisFile]
            self.open_file(self.currFileName)
            self.file_num.SetValue(str(self.thisFile))
        # display buttons for easy browsing next and previous through the list
        else:
            print "no files found"
            return
        self.updateStatusBar()

    def includeOrExcludeDlg(self):
        print "in dlg"
        dlg = wx.MessageDialog(self, 
             "Exclude this file?",
             "Exclude", wx.YES_NO|wx.CANCEL|wx.ICON_QUESTION)
        result = dlg.ShowModal()
        print result
        dlg.Destroy()
        return result
 
    def on_save_plot(self, event):
        file_choices = "PNG (*.png)|*.png"
        
        dlg = wx.FileDialog(
            self, 
            message="Save plot as...",
            defaultDir=os.getcwd(),
            defaultFile="plot.png",
            wildcard=file_choices,
            style=wx.SAVE)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.canvas.print_figure(path, dpi=self.dpi)
            self.flash_status_message("Saved to %s" % path)
            self.updateStatusBar()

    def onKeyPress(self, event):
        keycode = event.GetKeyCode()
        print keycode
        #if keycode == ord('A'):
        #    print "you pressed the spacebar!"
        #    sound_file = "notation1.wav"
        #    sound=wx.Sound(sound_file)
        #    print(sound_file)
        #    sound.Play(wx.SOUND_ASYNC)
        event.Skip()

    def OnAbout(self,e):
        # Create a message dialog box
        dlg = wx.MessageDialog(self, " A simple WAV file spectrogram viewer\nDavid J. Perkel\nUniversity of Washington\nperkel@uw.edu", "About wavBrowse", wx.OK)
        dlg.ShowModal() # Shows it
        dlg.Destroy() # finally destroy it when finished.

    def on_Inspect(self, e):
        import wx.lib.inspection
        wx.lib.inspection.InspectionTool().Show()

    def on_exit(self, event):
        self.Destroy()
    
    def on_save_data(self, event):
        fname= self.currFileName + '_data.txt'
        np.savetxt(fname, np.transpose([self.times, self.data]), delimiter=',')
    
    def on_mark_presence(self, event):
        #Obtain values from TextCtrls for Noise Threshold, Chunk Size, Overlap Fraction, and Minimum Off Duration        
        noise_threshold = float(self.noise_threshold.GetValue())
        ChunkSize = float(self.ChunkSize.GetValue())
        OverlapFraction = float(self.OverlapFraction.GetValue())
        MinOffDuration = float(self.MinOffDuration.GetValue())
        
        #Create an array consisting of RMS values of data (so this function will also update Chunk Size and Overlap Fraction)
        rms_data = np.array(rms_new(self.data, ChunkSize, OverlapFraction))
        #Create corresponding time array
        rms_time = np.array(align_time(self.times, ChunkSize, OverlapFraction))

        #Create new array consisting of RMS values above the "Noise Threshold" Value
        thresholded_data = rms_data[rms_data > noise_threshold]
        #Create corresponding time array by targeting the corresponding indices of the data values
        valid_indices = np.where(rms_data > noise_threshold)
        thresholded_time = rms_time[valid_indices]
        
        #Find "On Times"
        #Create an array consisting of the first thresholded time (this is the first on time and the function below will not catch it)
        on_times = np.array([thresholded_time[0]])
        #Iterate through threhsolded times and also find the value of the previous "thresholded time"
        for time in np.nditer(thresholded_time):
            current_index_tuple = list(np.where(thresholded_time == time))
            current_index = int(current_index_tuple[0])
            previous_index = current_index - 1
            value_previous = thresholded_time[previous_index]
            #Calculate the difference between the current time and the previous time
            if time - value_previous > MinOffDuration :
                #if this difference is greater than the Minimum Off Duration currently set, add this to the array of "On Times"
                on_times = np.append(on_times, float(time))
                
        #Find "Off Times"
        off_times = np.array([])
        #Iterate through threhsolded times
        for time in np.nditer(thresholded_time):
            current_index_tuple = list(np.where(thresholded_time == time))
            current_index = int(current_index_tuple[0])
            #If this is not the last time in the array, find the value of the next "thresholded time" (this avoids an IndexError)
            if current_index < thresholded_time.shape[0] - 1:
                next_index = current_index + 1
                value_next = thresholded_time[next_index]
                #Calculate the difference between the previous time and the current time
                if value_next - time > MinOffDuration :
                    #If this difference is greater than the Minimum Off Duration currently set, add this to the array of "Off Times"
                    off_times = np.append(off_times, float(time))
        #Add the final thresholded time to the array of Off Times (this is the last off time and the function above will not catch it)
        off_times = np.append(off_times, thresholded_time[-1])
        
        #Pair the on and off times in a list
        zipped = zip(on_times, off_times)
        
        #Find all the times during which the bird is on the perch
        on_perch = []
        #Iterate through times
        for time in np.nditer(self.times):
            time = float(time)
            #Iterate through each on and off time pair
            for on_time, off_time in zipped:
                #If the time is between the on and off time, append to on_perch
                if time > on_time and time < off_time:
                    on_perch.append(time)        
        
        #Find the max data point on the y axis (subtract one from the ymax used for self.axes for formatting reasons)
        ymax = max([max(rms_data) for data in rms_data]) - 1
        ymax2 = max([max(rms_data) for data in rms_data])
        
        #Create an array full of the ymax value at the times that the bird is on the perch to plot on self.axes
        y_values = np.full(len(on_perch), ymax, dtype=np.int)
        #Create an array full of the ymax2 value at the times that the bird is on the perch to plot on self.axes2
        y_values2 = np.full(len(on_perch), ymax2, dtype=np.int)
        
        #Clear self.axes2
        self.clear_axes_two()  
        #Plot only the "thresholded" data
        self.axes2.plot(thresholded_time, thresholded_data, 'g.')
        #Plot a horizontal line on self.axes2 at the times during which the bird is on the perch
        self.axes2.plot(on_perch, y_values2, 'y.')
        self.axes2.set_xlim(self.axes.get_xlim())
        
        #Plot a horizontal line on self.axes at the times during which the bird is on the perch
        self.clear_axes()
        self.axes.plot(on_perch, y_values, 'y.')
        self.axes.plot(self.times, self.data, 'r')        
                        
        #Draw everything plotted above
        self.canvas.draw()
        
        
        #ADD "WOULD YOU LIKE TO SAVE INFO AS A TEXT FILE"        
        
        #Save as text file
        #Make a name for the file using the current file name
        fname2 = self.currFileName + '_onoff.txt'
        
        #Open the file in "write" mode (overwrite anything currently on that file so that the user can repeat this function to find the proper values and only final correct info will be in txt file)
        fname2_write = open(fname2, "w")  
        #Write the Duration of the Recording to the txt file
        fname2_write.write("Duration of Recording:  ")
        fname2_write.write(str(float(self.times[-1])))
        #Close file in "write" mode and reopen in "append" mode (in "write" mode, saving np arrays will overwrite everything currently written)
        fname2_write.close()
        fname2_append = open(fname2,'a')
        #Write the Number of On Bouts, the On Times, and the Off Times to the txt file
        fname2_append.write("\n\nNumber of On Bouts:  ")
        fname2_append.write(str(on_times.shape[0]))
        fname2_append.write("\n\nOn Times: \n")
        np.savetxt(fname2_append, on_times, delimiter=',')
        fname2_append.write("\nOff Times: \n")
        np.savetxt(fname2_append, off_times, delimiter=',')
        #Close the txt file
        fname2_append.close()

        #print thresholded_time
        #print off_times
        
    def on_quit(self, event):
        self.Destroy()
    
    def flash_status_message(self, msg, flash_len_ms=1500):
        self.statusbar.SetStatusText(msg)
        self.timeroff = wx.Timer(self)
        self.Bind(
            wx.EVT_TIMER, 
            self.on_flash_status_off, 
            self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)
    
    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')

    def ErrorDlg(self, s):
        dlg = wx.MessageDialog(self, s, 'Error', wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

    def WarnDlg(self, s):
        dlg = wx.MessageDialog(self, s, 'Warning', wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

if __name__ == '__main__':
    app = wx.App(False)
    app.frame = GraphFrame()
    app.frame.Show()
    app.MainLoop()