import tkinter as tk
from tkinter.filedialog import asksaveasfile
from glob import glob
import sys

import plotter as PL

from threading import Thread

import numpy as np

import time
# import datetime
from datetime import datetime

from collections import deque

from fft_funcs import extracter

import sensors as SS

# from reader import DataReader as Reader
from reader import DataReader,FakeReader,CoreReader,StatusError

datatypes = [int,int,float]
sensors = []
# sensors = {}

DEBUG = True

nt_plot = (500,3)
nt_save = None


default_reader = CoreReader(0)

default_list_init = '   Select   '

default_baud = 115200
# default_baud = 9600  # 115200

# print('MAIN: ',default_reader)
start = time.time()

class ExperimentApp(object):
    def __init__(self, root,reader=default_reader):
        self.root = root
        self.root.title("2.20 Lab 2024 -- Wave-Body Interaction")

        self.threads = {}

        # Needs to be a mutable variable
        self._sensor_reader = {'reader':reader}
        # print(self._sensor_reader)

        self.record_bag = []
        self.recording = False

        self.data_labels = ['Body Position','PTO Position','Current']
        self.data_labels_short = ['z{body} [cm]','z{PTO} [cm]','I [mA]']

        # self.set_layout()
        self.layout = Layout(self)
        self.canvas = self.layout.canvas

        self.report(self._sensor_reader)

        # self.update_plots()

        self.power_buffer = deque(maxlen=100)  # Store last 100 power calculations
        self.fft_buffer = deque(maxlen=1000)  # (t, x, I)
        self.last_fft_time = time.time()
        self.fft_interval = 1.0  # seconds

        return

    @property
    def sensor_reader(self):
        return self._sensor_reader['reader']

    def connect_selected_reader(self,dev=None):  # reader
        if dev is None:
            dev = self.layout.device.get()

        if dev == default_list_init:
            self._sensor_reader['reader'] = default_reader

        if dev == 'fake':
            reader = FakeReader(nt_plot)
        else:
            baud = self.layout.devicebaud.get()
            reader = DataReader(nt_plot,
                                SS.SensorReader(dev,datatypes,baudrate=baud))
                                    # timeout=0,writeTimeout=0))
            self.report(reader.reader.serial)
            reader.reader.serial.flush()
            # reader.connect()
        self.report(f'Reader ID: {hex(id(reader))}')

        reader.set_scales([sensor.scaling_factor for sensor in sensors])

        self.connect_reader(reader)
        return

    def connect_reader(self,reader):
        self._sensor_reader['reader'] = reader

        self.report(f'self._sensor_reader ID: {id(self._sensor_reader["reader"])}')
        # SS.SensorReader("/dev/cu.usbmodem14201",[int,int],baudrate=115200)

        # if not self.sensor_reader.is_ready:
        #     self.sensor_reader.connect()
        self.sensor_reader.connect()

        # self.threads[self.sensor_reader.port] = Thread()
        return

    def report(self,message,newline=True):
        pos = '1.0'  # line 1, col 0
        # pos = tk.END   # insert at the end
        msg = message
        if newline:
            try:
                msg = message+'\n'
                # self.layout.stream_box.insert(pos,message+'\n')
            except Exception:
                msg = str(message)+'\n'
                # self.layout.stream_box.insert(pos,str(message)+'\n')
        self.layout.stream_box.insert(pos,msg)
        return

    def get_readings(self):
        if self.sensor_reader.is_ready:
            # self.report('in update_app')
            # self.report(f'in update_app: self._sensor_reader ID: {id(self._sensor_reader["reader"])}')
            # self.report(f'in update_app: self.sensor_reader ID: {id(self.sensor_reader)}')
            self.sensor_reader.read()
            vals = self.sensor_reader.scaled_values
            # vals = self.sensor_reader.values

            if self.layout.show_raw_data.get():
                raw_vals = self.sensor_reader.values
                vals[:,[0,1]] = raw_vals[:,[0,1]]

            t = self.sensor_reader.read_times

            # print(vals[-1])
            self.report(vals[-1])

            # if DEBUG and isinstance(self.sensor_reader,DataReader):
            #     self.report(self.sensor_reader.reader.serial.readline().strip(),newline=False)

            # print("Times and Position:", t, vals)

            return t,vals
        else:
            raise StatusError('Reader not ready')

    def update_app(self):

        delay = 10

        try:
            t,vals = self.get_readings()
            # Assuming:
            # - vals[:,0] is body position
            # - vals[:,2] is current (current)

            for ti, xi, Ii in zip(t, vals[:,0], vals[:,2]):
                self.fft_buffer.append((ti, xi, Ii))

            # Every second, run enhanced FFT analysis
            now = time.time()
            if now - self.last_fft_time >= self.fft_interval and now - start >= 10:
                self.last_fft_time = now
                self.run_fft_power_analysis()
        except StatusError:
            pass
        else:
            if self.power_buffer:
                power_times, power_vals = zip(*self.power_buffer)
                power_times = np.array(power_times)
                power_vals = np.array(power_vals)
            else:
                power_times = np.array([])
                power_vals = np.array([])
            # Update plots
            self.canvas.update_live(t,(power_times, power_vals), vals, self.plot_bools,
                                    raw_data=self.layout.show_raw_data.get())
            # self.canvas.update_live(vals,plot_bools)

            # # Save readings
            # if self.recording:
            #     self.record(t,vals)

        self.root.after(delay,self.update_app)

        return

    # def record(self):
    #     self.sensor_reader.read_times

    def run_fft_power_analysis(self):
        if len(self.fft_buffer) < 20:
            return

        times, positions, currents = zip(*self.fft_buffer)
        # print(self.fft_buffer)
        times = np.array(times)
        positions = np.array(positions) / 100.0  # Convert cm to meters
        currents = np.array(currents) / 1000.0   # Convert mA to A
        # times -= times[0]  # Normalize time

        if len(times) >= 100:
            # Check time intervals
            time_window = times[-100:]
            time_intervals = np.diff(time_window)
            min_interval = np.min(time_intervals)
            max_interval = np.max(time_intervals)
            avg_interval = np.mean(time_intervals)
            
            if min_interval < 0.001:  # Less than 1ms between samples is suspicious
                self.report(f"Warning: Very small time interval detected: {min_interval:.6f}s")
                return
                
            if max_interval / min_interval > 10:  # Suspicious time interval variation
                self.report(f"Warning: Large time interval variation - min: {min_interval:.3f}s, max: {max_interval:.3f}s, avg: {avg_interval:.3f}s")
                return  # Skip this FFT calculation to avoid potential spikes
            
            # Create evenly spaced time points for interpolation
            desired_sample_rate = 50  # Hz
            desired_interval = 1.0 / desired_sample_rate
            t_uniform = np.linspace(time_window[0], time_window[-1], 
                                  int((time_window[-1] - time_window[0]) / desired_interval))
            
            # Interpolate position data to uniform time points
            from scipy.interpolate import interp1d
            position_interp = interp1d(time_window, positions[-100:], kind='cubic', bounds_error=False)
            positions_uniform = position_interp(t_uniform)
            
            # For current, use linear interpolation since it might be noisier
            current_interp = interp1d(time_window, currents[-100:], kind='linear', bounds_error=False)
            currents_uniform = current_interp(t_uniform)
            
            # Use interpolated data for calculations
            mean_current = np.mean(currents_uniform)
            norm_times = t_uniform - t_uniform[0]
            
            # Debug interpolation quality
            if np.any(np.isnan(positions_uniform)):
                self.report("Warning: Interpolation produced NaN values")
                return
                
            # Apply Hanning window to reduce spectral leakage
            window = np.hanning(len(positions_uniform))
            positions_windowed = positions_uniform * window
            
            # Compensate for window amplitude reduction
            # Hanning window reduces amplitude by ~0.5, so multiply by 2
            window_correction = .4/.294
            
            amplitudes, frequencies = extracter(positions_windowed, norm_times, k=1)
            amplitudes = [amp * window_correction for amp in amplitudes]  # Apply correction

            if not amplitudes or not frequencies:
                self.report("FFT: No significant frequency detected.")
                return

            # Apply frequency stability check
            if hasattr(self, 'last_freq') and self.last_freq is not None:
                freq_change = abs(frequencies[0] - self.last_freq)
                if freq_change > 1:  # More than 1 Hz change
                    self.report(f"Warning: Large frequency change: {freq_change:.2f} Hz")
                    return
            self.last_freq = frequencies[0]

            damping = 1.5 * mean_current**2
            power = sum(0.5*damping*(amp**2)*(2*np.pi*f)**2 for amp, f in zip(amplitudes, frequencies))

            # Store the absolute time for the power value
            self.power_buffer.append((self.fft_buffer[-1][0], power))

            self.report(f"Power: {power:.4f} W | Freq: {frequencies[0]:.2f} Hz | Amp: {amplitudes[0]:.4f} m | Damping: {damping} N/m/s")
            # return power

    def toggle_recording(self):
        # print('In toggle_recording start',self.recording,self.layout)
        self.recording = not self.recording
        # print('In toggle_recording end',self.recording,self.layout)
        last_read_time = self.sensor_reader.read_times[-1]
        if self.recording:
            self.record_start_time = last_read_time
            self.record_bag = []
        else:
            self.record_end_time = last_read_time
            t,vals = self.get_readings()
            self.record_bag = np.column_stack((t,vals))
        return

    # def update_plots(self):

    #     if self.sensor_reader.is_ready:
    #         # self.report('in update_plots')
    #         # self.report(f'in update_plots: self._sensor_reader ID: {id(self._sensor_reader["reader"])}')
    #         # self.report(f'in update_plots: self.sensor_reader ID: {id(self.sensor_reader)}')
    #         self.sensor_reader.read()
    #         vals = self.sensor_reader.scaled_values
    #         # vals = self.sensor_reader.values

    #         t = self.sensor_reader.read_times

    #         # print(vals[-1])
    #         self.report(vals[-1])

    #         # if DEBUG and isinstance(self.sensor_reader,DataReader):
    #         #     self.report(self.sensor_reader.reader.serial.readline().strip(),newline=False)

    #         plot_bools = self.plot_bools
    #         # self.canvas.update_live(vals,plot_bools)
    #         self.canvas.update_live(t,vals,plot_bools)

    #         if self.recording:
    #             self.

    #     self.root.after(40,self.update_plots)

    #     return

    @property
    def plot_bools(self):
        plot_bools = [bool(self.layout.plot_pot1_line.get()),
                      bool(self.layout.plot_pot2_line.get()),
                      bool(self.layout.plot_current_line.get())]
        return plot_bools

    def save_data(self,f,fmt='%10.2f'):
        t = self.record_bag[:,0]
        tmask = ~np.isnan(t)
        A = np.array(self.record_bag[tmask])  # All non-NaN measurements
        print(A)

        t = A[:,0]
        twindow = np.logical_and(t >= self.record_start_time,
                                 t <= self.record_end_time)
        # print('In app.save_data. ',self.record_bag.shape,A.shape)
        # print(twindow)
        A = A[twindow]
        A[:,0] -= A[0,0]  # initialize zero-time
        # print(A)
        # np.savetxt(f,A[twindow,[True]+self.plot_bools],fmt=fmt)
        # np.savetxt(f,A[:,[True]+self.plot_bools][twindow],fmt=fmt)
        data_flag = tuple([True]+self.plot_bools)
        # print(data_flag,A[:,data_flag])

        field_size = 10  # default header field size
        try:
            field_size = int(fmt[1:].split('.')[0])
        except ValueError:
            pass

        ncols = np.count_nonzero(self.plot_bools)
        data_col_labels = [label for label,mask in
                           zip(self.data_labels_short,self.plot_bools)
                           if mask]

        # +1 here is to account for the delimiter
        fmt_str = '{:>'+str(field_size+1)+'}'
        data_header = (ncols*fmt_str).format(*data_col_labels)

        # -2 here is to account for 2 characters for the comment sign
        fmt_str = '{:>'+str(field_size-2)+'}'
        print(field_size,str(field_size),fmt_str)
        header = fmt_str.format('t [s] ') + data_header + \
            '\n' + (ncols*(field_size+1))*'-'

        np.savetxt(f,A[:,data_flag],fmt=fmt,header=header)
        # np.savetxt(f,A[:,[True]+self.plot_bools],fmt=fmt,header=header)
        return


class Layout(object):
    def __init__(self,app):
        self.app = app
        self.root = app.root
        self.reader = app._sensor_reader
        self.report = self.app.report
        self.create()
        self.set()
        # self.set_grid()
        self.update()
        return

    def reader_func_eval(self,fname,*args,**kwargs):
        # print('reader_func_eval: ',fname,*args)
        return getattr(self.app.sensor_reader,fname)(*args,**kwargs)

    def reader_func_factory(self,fname,*args,**kwargs):
        # print('reader_func_factory: ',fname,*args)

        # if 'widget' in kwargs:
        #     args = args + (float(kwargs.pop('widget').get()),)

        def func():
            # foo = tuple(args)
            # if 'widget' in kwargs:
            #     foo = foo + (float(kwargs.get('widget').get()),)

            # getattr(self.app.sensor_reader,fname)(*foo)
            f = getattr(self.app.sensor_reader,fname)

            # Adjust label if necessary
            # i = foo[0]
            i = args[0]
            sensor = sensors[i]

            # if fname == 'set_as_zero_offset_i' and i in [0,1]:
            #     zero_offset = self.app.sensor_reader.zero_offset[i]
            #     key = self.sensor_labels[i][0]
            #     status_label = getattr(self,f'{key}_status_label')
            #     amp = sensor.get_centered_range()
            #     txt = (f'Zero set at {100*zero_offset:.0f}% pot position. '
            #            f'Available +/- amplitude is {amp:.1f} cm.')
            #     status_label.config(text=txt)

            if fname == 'set_as_zero_offset_i':
                f(i)
                if i in [0,1]:
                    self.set_status_label_text(i)
                    # zero_offset = self.app.sensor_reader.zero_offset[i]
                    # sensor.center_point = zero_offset
                    # key = self.sensor_labels[i][0]
                    # status_label = getattr(self,f'{key}_status_label')
                    # # amp = sensor.get_centered_range(zero_offset)
                    # # txt = (f'Zero set at {100*zero_offset:.0f}% pot position. '
                    # #        f'Available +/- amplitude is {amp:.1f} cm.')
                    # txt = sensor.report_range()
                    # status_label.config(text=txt)

            elif fname == 'set_scale_i':
                val = float(kwargs.get('widget').get())
                # Entry is a correction to the default scaling factor
                val *= sensor.scaling_factor
                f(i,val)

            return  # getattr(self.app.sensor_reader,fname)(*foo) #,**kwargs)
        return func

    def create(self):
        root = self.root
        reader = self.app.sensor_reader
        self.widgets = {}  # dict for minor widgets

        self.canvas = PL.MplCanvas(root)
        self.canvas.setup()

        self.canvas.draw()

        # ------------------  Device Settings --------------------

        device_frame = tk.Frame(root)

        self.widgets['device_frame'] = device_frame

        r_circ = 10
        self.devstatus = tk.Canvas(device_frame,height=3*r_circ,width=3*r_circ)

        circ_args = (r_circ,)*2 + (2*r_circ,)*2
        # circ_args = (0.0,)*2 + (r_circ,)*2

        self.devstatus_i = self.devstatus.create_oval(*circ_args,fill='blue')

        device_address_label = tk.Label(device_frame, text="Select device")
        self.widgets['device_address_label'] = device_address_label

        # devices = glob('/dev/cu*usb*')
        if sys.platform.startswith('win'):
            devices = ['COM%s' % i for i in [5,6,7,8,9]]
        elif sys.platform.startswith('linux') or sys.platform.startswith('cygwin'):
            # this excludes your current terminal "/dev/tty"
            devices = glob('/dev/tty[A-Za-z]*')
        elif sys.platform.startswith('darwin'):
            # devices = glob('/dev/tty.*')
            devices = glob('/dev/cu*usb*')
        else:
            raise EnvironmentError('Unsupported platform')

        if DEBUG:
            devices.append('fake')

        self.device = tk.StringVar(device_frame)
        self.device.set(default_list_init)  # default value
        self.device.set(devices[0])  # default value
        self.devlist = tk.OptionMenu(device_frame, self.device,*devices)

        # for i,dev in enumerate(devices):
        #     self.devlist.insert (i+1,dev)

        self.devconnectButton = tk.Button(device_frame,text='Connect',
                                          command=self.device_connection)

        bauds = [9600,115200]
        self.devicebaud = tk.IntVar(device_frame)
        self.devicebaud.set(default_baud)  # default value
        self.baudlist = tk.OptionMenu(device_frame,self.devicebaud,*bauds)

        # self.button_quit = tk.Button(master=root, text="Quit",
        #                              command=root.destroy)

        # ------------------  Plot Settings --------------------

        plot_checkbox_frame = tk.Frame(root)
        self.widgets['plot_checkbox_frame'] = plot_checkbox_frame

        self.time_window_var = tk.StringVar(value=str(self.canvas.time_window))
        self.time_window_entry = tk.Entry(plot_checkbox_frame,width=6,
                                          textvariable=self.time_window_var)
        self.time_window_button = tk.Button(plot_checkbox_frame,text='Set',
                                            command=self.set_time_window)

        self.show_raw_data = tk.IntVar(value=0)
        self.raw_data_checkbox = tk.Checkbutton(plot_checkbox_frame,text='Raw data',
                                                variable=self.show_raw_data)

        default_checkbox_state = 1  # True  # 'selected'
        self.plot_pot1_line = tk.IntVar(value=default_checkbox_state)
        self.plot_pot2_line = tk.IntVar(value=0)
        self.plot_current_line = tk.IntVar(value=default_checkbox_state)

        self.pot1_checkbox = tk.Checkbutton(plot_checkbox_frame,text=self.app.data_labels[0], #'Plot Pot 1',
                                            variable=self.plot_pot1_line)
        self.pot2_checkbox = tk.Checkbutton(plot_checkbox_frame,text=self.app.data_labels[1], #'Plot Pot 2',
                                            variable=self.plot_pot2_line)

        self.current_checkbox = tk.Checkbutton(plot_checkbox_frame,text=self.app.data_labels[2], #'Plot current',
                                             variable=self.plot_current_line)

        self.reset_plot1 = tk.Button(plot_checkbox_frame,text='Reset Plot 1',
                                     command=lambda: self.canvas.initialize_position_plot(init=False))

        self.reset_plot2 = tk.Button(plot_checkbox_frame,text='Reset Plot 2',
                                     command=self.canvas.initialize_current_plot)


        self.record_button_actions = {True:'Stop',False:'Start'}
        self.record_button_colors = {True:'red',False:'green'}
        action = self.record_button_actions[self.app.recording]
        clr = self.record_button_colors[self.app.recording]
        rec_butt_text = f'{action} Recording'
        self.record_button_var = tk.StringVar(value=rec_butt_text)

        self.record_button = tk.Button(plot_checkbox_frame,text=rec_butt_text,
                                       bg=clr,command=self.toggle_recording)


        self.save_button = tk.Button(plot_checkbox_frame, text= "Save Data",
                                     command=self.save_file)

        # ------------------  Calibration Settings --------------------

        calibration_frame = tk.Frame(root,relief='sunken')

        # self.calibration = {''}

        # calibration_frame_pot1 = tk.Frame(calibration_frame)

        self.widgets['calibration_frame'] = calibration_frame
        # self.widgets['calibration_frame_pot1'] = calibration_frame_pot1

        # self.sensor_labels = {'pot1':'Body Position','pot2':'PTO position',
        #                       'current':'Current'}

        self.sensor_labels = [('pot1','Body Position'),('pot2','PTO position'),
                              ('current','Current')]

        # for key,label in self.sensor_labels.items():
        for i,(key,label) in enumerate(self.sensor_labels):
            # print(i,key,label)
            sub_calibration_frame = tk.Frame(calibration_frame,bd=2,
                                             highlightbackground='darkgray',
                                             highlightthickness=1)

            button = tk.Button(sub_calibration_frame,
                               command=self.reader_func_factory('set_as_zero_offset_i',i),
                               # command=lambda: self.reader_func_eval('set_as_zero_offset_i',i),
                               # command=lambda: reader.set_as_zero_offset_i(i),
                               text='Set Zero')

            # setvar = tk.DoubleVar(value=1.0)
            setvar = tk.StringVar(value=str(1.0))
            entry = tk.Entry(sub_calibration_frame,textvariable=setvar,width=8)

            cal_button = tk.Button(sub_calibration_frame,
                                   command=self.reader_func_factory('set_scale_i',i,widget=entry),
                                   # command=lambda: self.reader_func_eval('set_scale_i',i,float(entry.get())),
                                   # command=lambda: reader.set_scale_i(i,entry.get()),
                                   text='Set Calibration Factor')

            lbl = tk.Label(sub_calibration_frame,text=label)

            self.widgets[f'calibration_frame_{key}'] = sub_calibration_frame
            setattr(self,f'{key}_label',lbl)
            setattr(self,f'{key}_zero_button',button)
            setattr(self,f'{key}_calibration_var',setvar)
            setattr(self,f'{key}_calibration_entry',entry)
            setattr(self,f'{key}_calibration_button',cal_button)

            if key != 'current':
                # sensor = sensors[i]

                # zero_offset = self.app.sensor_reader.zero_offset[i]
                # sensor.center_point = zero_offset
                # txt = sensor.report_range()
                status_label = tk.Label(sub_calibration_frame,anchor='w',
                                        justify=tk.LEFT,width=23) #,text=txt)
                setattr(self,f'{key}_status_label',status_label)
                self.set_status_label_text(i)

            if key == 'pot1':
                if self.pot1_zero_button is not button:
                    print('ERROROR!!!! ZERO BUTTON ID DOES NOT MATCH')

                if self.pot1_calibration_button is not cal_button:
                    print('ERROROR!!!! CAL BUTTON ID DOES NOT MATCH')

                if self.pot1_calibration_entry is not entry:
                    print('ERROROR!!!! ENTRY ID DOES NOT MATCH')

        # self.pot1_zero_button = tk.Button(calibration_frame,
        #                                   # command=set_zero,
        #                                   text='Pot 1: Set Zero')
        # self.pot2_zero_button = tk.Button(calibration_frame,
        #                                   # command=set_zero,
        #                                   text='Pot 2: Set Zero')

        # self.pot1_calibration_var = tk.DoubleVar(value=1.0)
        # self.pot1_calibration_entry = tk.Entry(calibration_frame_pot1,
        #                                        textvariable=self.pot1_calibration_var)

        # self.pot1_calibration_button = tk.Button(calibration_frame_pot1,
        #                                   # command=set_calibration,
        #                                   text='Pot 1: Set Calibration Factor')

        # self.slider_update = tk.Scale(root, from_=1, to=5,orient=tk.HORIZONTAL,
        #                               command=self.canvas.update_frequency,
        #                               label="Frequency [Hz]")

        # stream_row = tk.Frame(root, relief='sunken')
        self.stream_box = tk.Text(root, height=7, width=200,fg='blue')
        # self.stream_box.insert(tk.END,'asdfasdfas\n')

        return

    def set(self):

        # ------------------  Device --------------------
        self.widgets['device_frame'].pack(side=tk.TOP)
        self.devstatus.pack(side=tk.LEFT,fill=tk.Y)

        self.widgets['device_address_label'].pack(side=tk.LEFT)
        self.devlist.pack(side=tk.LEFT)
        self.baudlist.pack(side=tk.LEFT)

        self.devconnectButton.pack(side=tk.LEFT)

        # ------------------  Calibration --------------------

        self.widgets['calibration_frame'].pack(side=tk.TOP)

        self.widgets['calibration_frame_pot1'].pack(side=tk.LEFT)

        self.pot1_zero_button.pack(side=tk.TOP)
        self.pot1_calibration_entry.pack(side=tk.TOP)
        self.pot1_calibration_button.pack(side=tk.TOP)
        # self.pot2_zero_button.pack(side=tk.TOP)

        # for key,label in self.sensor_labels.items():
        for key,label in self.sensor_labels:

            self.widgets[f'calibration_frame_{key}'].pack(side=tk.LEFT)
            getattr(self,f'{key}_label').pack(side=tk.TOP)
            getattr(self,f'{key}_zero_button').pack(side=tk.LEFT)
            getattr(self,f'{key}_calibration_entry').pack(side=tk.LEFT)
            getattr(self,f'{key}_calibration_button').pack(side=tk.LEFT)

            if key != 'current':
                getattr(self,f'{key}_status_label').pack(side=tk.BOTTOM)

        # ------------------  Plot settings --------------------

        self.widgets['plot_checkbox_frame'].pack(side=tk.TOP)

        self.raw_data_checkbox.pack(side=tk.LEFT)

        self.time_window_entry.pack(side=tk.LEFT)
        self.time_window_button.pack(side=tk.LEFT)

        self.pot1_checkbox.pack(side=tk.LEFT)
        self.pot2_checkbox.pack(side=tk.LEFT)
        self.current_checkbox.pack(side=tk.LEFT)

        self.reset_plot1.pack(side=tk.LEFT)
        self.reset_plot2.pack(side=tk.LEFT)

        self.record_button.pack(side=tk.LEFT)
        self.save_button.pack(side=tk.LEFT)

        # Packing order is important. Widgets are processed sequentially and if there
        # is no space left, because the window is too small, they are not displayed.
        # The canvas is rather flexible in its size, so we pack it last which makes
        # sure the UI controls are displayed as long as possible.

        # self.button_quit.pack(side=tk.BOTTOM)
        # self.slider_update.pack(side=tk.BOTTOM)

        self.stream_box.pack(side=tk.BOTTOM)
        # self.widgets['device_frame'].pack(side=tk.BOTTOM)
        # toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        return

    def set_grid(self):

        # ------------------  Root elements --------------------

        self.widgets['device_frame'].grid(row=0, column=0)
        self.canvas.get_tk_widget().grid(row=1, column=0)
        self.stream_box.grid(row=2, column=0)

        # ------------------  device_frame elements --------------------
        # self.devstatus.grid(row=0, column=0, rowspan=2)
        self.devstatus.grid(row=1, column=0)

        self.widgets['device_address_label'].grid(row=0, column=1)
        self.devlist.grid(row=1, column=1)

        self.baudlist.grid(row=1, column=2)

        self.devconnectButton.grid(row=1, column=3)

        self.pot1_checkbox.grid(row=0, column=1)
        self.pot2_checkbox.grid(row=0, column=2)
        self.current_checkbox.grid(row=0, column=3)

        self.reset_plot1.grid(row=0, column=0)
        self.reset_plot2.grid(row=0, column=4)

        return

    def device_connection(self):
        reader = self.reader.get('reader')
        # if status == 'disconnected'
        if self.device.get() == default_list_init:
            self.reader['reader'] = default_reader
            self.report('Setting default reader')
            # raise ValueError('Select a device first')
            return
        # self.app.connect_selected_reader()
        if reader is default_reader:
            # reader not connected yet, but there is something to connect to
            self.report(f'Connecting from default reader to {self.device.get()}')
            self.app.connect_selected_reader()
            # reader.start()
            # self.devconnectButton.config(text='Disconnect',fg='red')
        else:
            # reader.stop()
            self.report(reader)
            # reader.stop_thread = True
            reader.disconnect()
            self.reader['reader'] = default_reader
            self.report('Stopping and setting default reader')
        return

    def toggle_recording(self):
        self.app.toggle_recording()
        # button_colors = {True:'red',False:'green'}
        action = self.record_button_actions[self.app.recording]
        label = f'{action} Recording'
        clr = self.record_button_colors[self.app.recording]
        # print('In toggle_recording layout: ',label)
        #
        # NOTE: Setting background color (e.g. `bg=clr`)
        # does not seem to be working on macOS
        # https://stackoverflow.com/a/1530913
        self.record_button.config(text=label,bg=clr)

        if self.app.recording:
            self.devconnectButton.config(state='disabled')
        else:
            self.devconnectButton.config(state='normal')
        return

    def save_file(self):
        filetypes = [("All Files","*.*"),("Text Documents","*.txt")]
        f = asksaveasfile(mode='w',initialfile = 'Results.txt',
                          defaultextension=".txt",
                          filetypes=filetypes)
        if f is None:
            # asksaveasfile return `None` if dialog closed with "cancel".
            return
        self.app.save_data(f)
        f.close()
        return

    def set_time_window(self):
        tw = float(self.time_window_entry.get())
        self.canvas.time_window = tw
        self.canvas.update_xaxes()
        self.canvas.draw()
        return

    def set_status_label_text(self,i):
        sensor = sensors[i]
        key = self.sensor_labels[i][0]

        try:
            zero_offset = self.app.sensor_reader.zero_offset[i]
        except TypeError:
            txt = 'foo'
        else:
            sensor.center_point = zero_offset
            txt = sensor.report_range()

        status_label = getattr(self,f'{key}_status_label')
        status_label.config(text=txt)
        return

    def update(self):
        # print('boo: ' + datetime.now().strftime('%H:%M:%S.%f'))

        reader = self.reader.get('reader')

        # if self.serial is not None:
        if self.device.get() != default_list_init:
            # connected =  self.serial.isOpen()
            # connected = self.slider_update.get() == 2
            connected = reader.is_ready
            if connected:
                self.devconnectButton.config(text='Disconnect',fg='red')
                self.devstatus.itemconfig(self.devstatus_i,fill='green')
                self.devlist.config(state='disabled')
                self.baudlist.config(state='disabled')
            else:
            # elif self.slider_update.get() == 3:
                self.devconnectButton.config(text='Connect',fg='green')
                self.devstatus.itemconfig(self.devstatus_i,fill='red')
                self.devlist.config(state='normal')
                self.baudlist.config(state='normal')
        else:
            self.devconnectButton.config(text='Connect')
            self.devstatus.itemconfig(self.devstatus_i,fill='gray')
            self.devlist.config(state='normal')

        self.root.after(10,self.update)


devconnectButton_params = dict(text='Disconnect',fg='red')

# def Entry2func(entry,entry_var,func):
#     new


class Circ2Lin(object):
    def __init__(self,diameter,max_turns,center_point=0.0,label=None,unit='cm'):
        self.diameter = diameter
        self.max_turns = max_turns
        self._center_point = center_point
        self.label = label
        self.unit = unit
        return

    @property
    def scaling_factor(self):
        return self.get_range()

    def get_range(self):
        return get_range(self.diameter,self.max_turns)

    def get_centered_range(self,center_point=None):
        if center_point is None:
            center_point = self.center_point
        return get_centered_range(center_point,self.diameter,self.max_turns)

    @property
    def center_point(self):
        return self._center_point

    @center_point.setter
    def center_point(self,v):
        if not (0.0 <= v <= 1.0):
            raise ValueError('Zero point outside [0,1] interval. ',v)
        self._center_point = v
        return

    def report_range(self):  # ,center_point=None):
        zero_offset = self.center_point
        amp = self.get_centered_range()
        txt = (f'Zero set at {100*zero_offset:.0f}% pot position. \n'
               f'Available +/- amplitude is {amp:.1f} {self.unit}.')
        return txt


class Lin(object):
    def __init__(self,scaling_factor=1.0,center_point=0.0,label=None):
        self.scaling_factor = scaling_factor


def get_range(diameter,max_turns):
    pi = np.pi  # 3.14159
    return diameter*pi*max_turns


def get_centered_range(center_point,diameter,max_turns):
    if not (0.0 <= center_point <= 1.0):
        raise ValueError('Zero point outside [0,1] interval. ',center_point)

    rng = get_range(diameter,max_turns)
    return min(center_point*rng,(1-center_point)*rng)


if __name__ == "__main__":
    import toml
    import sys

    with open('params.toml','r') as f:
        params = toml.load(f)

    for key,ps in params['sensor'].items():
        if key == 'pot1':  # Time of flight sensor
            item = Lin(ps.get('scaling_factor',1.0), label=key)
        elif key == 'pot2':  # Still using rotary sensor
            try:
                item = Circ2Lin(ps['diameter'],ps['max_turns'],label=key)
            except KeyError:
                item = Lin(ps.get('scaling_factor',1.0),label=key)
        else:  # Current sensor
            item = Lin(ps.get('scaling_factor',1.0),label=key)

        # sensors[key] = item
        sensors.append(item)

    print(sensors)

    root = tk.Tk()

    app = ExperimentApp(root)

    # main_thread = Thread(target=root.mainloop)

    # main_thread = Thread(target=app.update_plots)
    # main_thread = Thread(target=app.update_app)
    # main_thread.start()

    app.update_app()
    # app.update_plots()

    print('Starting GUI...')
    root.mainloop()
