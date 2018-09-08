# outputfile_converter.py - converts gprMax merged output HDF5 file to RD3, DZT, DT1 and IPRB files
#
# Author: Dimitrios Angelis
# Copyright: 2017-2018
# Last modified: 05/09/2018


# Libraries ============================================================================================================
import os
import sys
import struct
import bitstruct
import datetime
import h5py as h5
import numpy as np
import tkinter as tk
import matplotlib.pyplot as plt                                                 # v. 2.2.3

from scipy import signal
from bitstring import Bits
from tkinter import filedialog
from tkinter import messagebox


# Select file ==========================================================================================================
dlg = tk.Tk()
dlg.withdraw()
file_path = filedialog.askopenfilename(title='Select gprMax File',
                                       filetypes=(('HDF5 Files', '*.out'), ('all files', '*.*')))

if not file_path:
    erd = messagebox.showerror('Error!', 'No input file')
    dlg.destroy()
    sys.exit(0)

dlg.destroy()


# Read data from HDF5 file =============================================================================================
fid = h5.File(file_path, 'r')
dataex = np.array(fid.get('/rxs/rx1/Ex'))
dataey = np.array(fid.get('/rxs/rx1/Ey'))
dataez = np.array(fid.get('/rxs/rx1/Ez'))
fid.close()


# Field check ==========================================================================================================
if dataey.all == 0 and dataez.all == 0:
    data = dataex
elif dataex.all == 0 and dataez.all == 0:
    data = dataey
elif dataex.all == 0 and dataey.all == 0:
    data = dataez
else:
    maxex = np.max(dataex)
    maxey = np.max(dataey)
    maxez = np.max(dataez)
    if (maxex > maxey) and (maxex > maxez):
        data = dataex
    elif (maxey > maxex) and (maxey > maxez):
        data = dataey
    elif (maxez > maxex) and (maxez > maxey):
        data = dataez


# Additional information ===============================================================================================
class entry_info_window:
    def __init__(self, window):
        window.title('Additional Information')
        window.resizable(False, False)

        label1 = tk.Label(window, text='Waveform Centre Frequency (MHz)')
        label1.grid(row=1, column=1)
        label2 = tk.Label(window, text='Source-Receiver Distance (m)')
        label2.grid(row=2, column=1)
        label3 = tk.Label(window, text='Trace Interval / Step (m)')
        label3.grid(row=3, column=1)

        self.entry1 = tk.StringVar()
        self.entry2 = tk.StringVar()
        self.entry3 = tk.StringVar()
        box1 = tk.Entry(window, textvariable=self.entry1)
        box1.grid(row=1, column=2)
        box2 = tk.Entry(window, textvariable=self.entry2)
        box2.grid(row=2, column=2)
        box3 = tk.Entry(window, textvariable=self.entry3)
        box3.grid(row=3, column=2)

        self.ok_btn()
        self.cancel_btn()

    def ok_btn(self):
        self.ok_btn = tk.Button(window, text='OK', bg='green', fg='white', command=self.ok)
        self.ok_btn.grid(row=5, column=1)
        self.ok_btn.config(height=1, width=10)

    def cancel_btn(self):
        self.cancel_btn = tk.Button(window, text='Cancel', bg='red', fg='white', command=self.cancel)
        self.cancel_btn.grid(row=5, column=2)
        self.cancel_btn.config(height=1, width=10)

    def ok(self):
        self.answers = [self.entry1.get(), self.entry2.get(), self.entry3.get()]
        window.destroy()

    def assign(self):
        return self.answers

    def cancel(self):
        self.answers = []
        window.destroy()


while True:
    if __name__ == '__main__':
        window = tk.Tk()
        Entry_Info_Window = entry_info_window(window)
        window.mainloop()
        answers = Entry_Info_Window.assign()

    if not answers:
        sys.exit(0)
    else:
        try:
            answers = [float(i) for i in answers]
        except ValueError:
            answers = []

        if answers == [] or answers[0] <= 0 or answers[1] < 0 or answers[2] <= 0:
            continue
        else:
            break


# Create header with the basic information =============================================================================
class structure:
    pass


HDR = structure()
HDR.fname = os.path.basename(file_path)                                         # File Name
HDR.pname = os.path.dirname(file_path)                                          # Path
HDR.fname, HDR.fext = os.path.splitext(HDR.fname)                               # File extension
HDR.centre_freq = answers[0]                                                    # Centre frequency (MHz)
HDR.ant_sep = answers[1]                                                        # Antenna seperation / Tx-Rx distance (m)
HDR.trac_int = answers[2]

fid = h5.File(file_path, 'r')
HDR.samp_int = fid.attrs['dt'] * 10**9                                          # Sampling interval (ns)
fid.close()

HDR.samp_freq = (1 / HDR.samp_int) * 10**3                                      # Sampling frequency (MHz)
HDR.num_samp, HDR.num_trac = np.shape(data)                                     # Number of samples and traces
HDR.time_window = HDR.num_samp * HDR.samp_int                                   # Time window (ns)
HDR.antenna = ('gprMax ' + str(HDR.centre_freq) + 'MHz')                        # Antenna name


# **************************************************** Optional ********************************************************
# Resample to 1024 samples =============================================================================================
# I usually perform this step for either 512 or 1024 samples (line 162) because many GPR processing software cannot load
# files with more samples.
data = signal.resample(data, 1024)                                              # <------- 1024 samples

HDR.num_samp = len(data)                                                        # New number of samples
HDR.samp_int = HDR.time_window / HDR.num_samp                                   # New sampling interval (ns)
HDR.samp_freq = (1 / HDR.samp_int) * 10**3                                      # New sampling frequency (MHz)


# **************************************************** Optional ********************************************************
# Data scale ===========================================================================================================
data = data * 32767.5 / np.max(np.abs(data))                       #signal * ((1 - 1 / 2^bitrate) * 32768) / max(signal)
data[data >= 32767] = 32767
data[data <= -32768] = -32768
data = np.round(data, decimals=0)


# Plots ================================================================================================================
x = np.linspace(0, (HDR.num_trac - 1) * HDR.trac_int, HDR.num_trac)             # Distance of each trace (m)
t = np.linspace(HDR.samp_int, (HDR.num_samp * HDR.samp_int), HDR.num_samp)      # TIme of each sample (ns)

# Bscan
f1 = plt.figure(num='Bscan', facecolor='white', edgecolor='white')
plt.imshow(data, extent=[np.min(x), np.max(x), np.max(t), np.min(t)], aspect='auto',
           cmap=plt.cm.bone, vmin=np.min(data), vmax=np.max(data), interpolation='nearest')

plt.xlabel('Distance (m)')
plt.ylabel('Time (ns)')
plt.title(HDR.fname, pad=40)
ax1 = f1.gca()
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position('top')
plt.draw()


# Frequency spectrum
def nextpow2(args):
    return int(np.ceil(np.log2(np.abs(args))))


m = 2**nextpow2(HDR.num_samp)

amp = np.array(np.fft.fft(data, m, axis=0))
amp = (np.abs(amp[: int(m / 2), :] / m)) * 2
amp = np.mean(amp, axis=1)

freq = np.array(np.fft.fftfreq(m, HDR.samp_int * 10**-9) * 10**-6)
freq = freq[: int(m / 2)]

f2 = plt.figure(num='Frequency Spectrum', facecolor='white', edgecolor='white')
plt.fill(freq, amp, 'black')
plt.ylabel('Amplitude')
plt.xlabel('Frequency (MHz)')
plt.axis([np.min(freq), np.max(freq), np.min(amp), np.max(amp)])
plt.draw()
plt.show(block=False)


# Export option ========================================================================================================
class entry_export_window:
    def __init__(self, window):
        window.title('Choose GPR Format')
        window.resizable(False, False)

        label1 = tk.Label(window, text='1 = RD3,  2 = DZT,  3 = DT1,  4 = IPRB')
        label1.grid(row=1, column=1)

        self.entry1 = tk.StringVar()

        box1 = tk.Entry(window, textvariable=self.entry1)
        box1.grid(row=1, column=2)

        self.ok_btn()
        self.cancel_btn()

    def ok_btn(self):
        self.ok_btn = tk.Button(window, text='OK', bg='green', fg='white', command=self.ok)
        self.ok_btn.grid(row=5, column=1)
        self.ok_btn.config(height=1, width=10)

    def cancel_btn(self):
        self.cancel_btn = tk.Button(window, text='Cancel', bg='red', fg='white', command=self.cancel)
        self.cancel_btn.grid(row=5, column=2)
        self.cancel_btn.config(height=1, width=10)

    def ok(self):
        self.answer = [self.entry1.get()]
        window.destroy()

    def assign(self):
        return self.answer

    def cancel(self):
        self.answer = []
        window.quit()


while True:
    if __name__ == '__main__':
        window = tk.Tk()
        Entry_Export_Window = entry_export_window(window)
        window.mainloop()
        gpr_format = Entry_Export_Window.assign()


    if not gpr_format:
        sys.exit(0)
    else:
        try:
            gpr_format = [int(i) for i in gpr_format]
        except ValueError:
            gpr_format = []

        if gpr_format == [] or gpr_format[0] != 1 and gpr_format[0] !=2 and gpr_format[0] != 3 and gpr_format[0] != 4:
            continue
        else:
            gpr_format = gpr_format[0]
            break


# RAD / RD3, Mala GeoScience ===========================================================================================
# Rad is the header file. In this file is all the important information such as the number of samples, traces,
# measurement intervals can be found.
# Rd3 is the data file. This file contains only the data (amplitude values) in a binary form.

if gpr_format == 1:
    # Header structure
    HDR.fname = HDR.fname                                                       # File name
    HDR.num_samp = HDR.num_samp                                                 # Number of samples
    HDR.samp_freq = HDR.samp_freq                                               # Sampling frequency (MHz)
    HDR.frequency_steps = 1                                                     # Frequency steps
    HDR.signal_pos = 0                                                          # Signal position
    HDR.raw_signal_pos = 0                                                      # Raw signal position
    HDR.distance_flag = 1                                                       # Distance flag: 0 time interval, 1 distance interval
    HDR.time_flag = 0                                                           # Time flag    : 0 distance interval, 1 time interval
    HDR.program_flag = 0                                                        # Program flag
    HDR.external_flag = 0                                                       # External flag
    HDR.trac_int_sec = 0                                                        # Trace interval in seconds(only if Time flag = 1)
    HDR.trac_int = HDR.trac_int                                                 # Trace interval in meters (only if Distance flag = 1)
    HDR.operator = 'Unknown'                                                    # Operator
    HDR.customer = 'Unknown'                                                    # Customer
    HDR.site = 'gprMax'                                                         # Site
    HDR.antenna = HDR.antenna                                                   # Antenna name
    HDR.antenna_orientation = 'NOT VALID FIELD'                                 # Antenna orientation
    HDR.ant_sep = HDR.ant_sep                                                   # Antenna seperation / Tx-Rx distance (m)
    HDR.comment = '----'                                                        # Comment
    HDR.time_window = HDR.time_window                                           # Time window (ns)
    HDR.stacks = 1                                                              # Stacks
    HDR.stack_exponent = 0                                                      # Stack exponent
    HDR.stacking_time = 0                                                       # Stacking Time
    HDR.num_trac = HDR.num_trac                                                 # Number of traces
    HDR.stop_pos = HDR.num_trac * HDR.trac_int                                  # Stop position (m)
    HDR.system_calibration = 0
    HDR.start_pos = 0                                                           # Start position (m)
    HDR.short_flag = 1
    HDR.intermediate_flag = 0
    HDR.long_flag = 0
    HDR.preprocessing = 0
    HDR.high = 0
    HDR.low = 0
    HDR.fixed_increment = 0
    HDR.fixed_moves_up = 0
    HDR.fixed_moves_down = 1
    HDR.fixed_position = 0
    HDR.wheel_calibration = 0
    HDR.positive_direction = 1

    with open(HDR.fname + '.rad', 'w') as fid:
        fid.write('%s%i\n' % ('SAMPLES:', HDR.num_samp))
        fid.write('%s%0.6f\n' % ('FREQUENCY:', HDR.samp_freq))
        fid.write('%s%i\n' % ('FREQUENCY STEPS:', HDR.frequency_steps))
        fid.write('%s%0.6f\n' % ('SIGNAL POSITION:', HDR.signal_pos))
        fid.write('%s%i\n' % ('RAW SIGNAL POSITION:', HDR.raw_signal_pos))
        fid.write('%s%i\n' % ('DISTANCE FLAG:', HDR.distance_flag))
        fid.write('%s%i\n' % ('TIME FLAG:', HDR.time_flag))
        fid.write('%s%i\n' % ('PROGRAM FLAG:', HDR.program_flag))
        fid.write('%s%i\n' % ('EXTERNAL FLAG:', HDR.external_flag))
        fid.write('%s%0.6f\n' % ('TIME INTERVAL:', HDR.trac_int_sec))
        fid.write('%s%0.6f\n' % ('DISTANCE INTERVAL:', HDR.trac_int))
        fid.write('%s%s\n' % ('OPERATOR:', HDR.operator))
        fid.write('%s%s\n' % ('CUSTOMER:', HDR.customer))
        fid.write('%s%s\n' % ('SITE:', HDR.site))
        fid.write('%s%s\n' % ('ANTENNAS:', HDR.antenna))
        fid.write('%s%s\n' % ('ANTENNA ORIENTATION:', HDR.antenna_orientation))
        fid.write('%s%0.6f\n' % ('ANTENNA SEPARATION:', HDR.ant_sep))
        fid.write('%s%s\n' % ('COMMENT:', HDR.comment))
        fid.write('%s%0.6f\n' % ('TIMEWINDOW:', HDR.time_window))
        fid.write('%s%i\n' % ('STACKS:', HDR.stacks))
        fid.write('%s%i\n' % ('STACK EXPONENT:', HDR.stack_exponent))
        fid.write('%s%0.6f\n' % ('STACKING TIME:', HDR.stacking_time))
        fid.write('%s%i\n' % ('LAST TRACE:', HDR.num_trac))
        fid.write('%s%0.6f\n' % ('STOP POSITION:', HDR.stop_pos))
        fid.write('%s%0.6f\n' % ('SYSTEM CALIBRATION:', HDR.system_calibration))
        fid.write('%s%0.6f\n' % ('START POSITION:', HDR.start_pos))
        fid.write('%s%i\n' % ('SHORT FLAG:', HDR.short_flag))
        fid.write('%s%i\n' % ('INTERMEDIATE FLAG:', HDR.intermediate_flag))
        fid.write('%s%i\n' % ('LONG FLAG:', HDR.long_flag))
        fid.write('%s%i\n' % ('PREPROCESSING:', HDR.preprocessing))
        fid.write('%s%i\n' % ('HIGH:', HDR.high))
        fid.write('%s%i\n' % ('LOW:', HDR.low))
        fid.write('%s%0.6f\n' % ('FIXED INCREMENT:', HDR.fixed_increment))
        fid.write('%s%i\n' % ('FIXED MOVES UP:', HDR.fixed_moves_up))
        fid.write('%s%i\n' % ('FIXED MOVES DOWN:', HDR.fixed_moves_down))
        fid.write('%s%0.6f\n' % ('FIXED POSITION:', HDR.fixed_position))
        fid.write('%s%0.6f\n' % ('WHEEL CALIBRATION:', HDR.wheel_calibration))
        fid.write('%s%i\n' % ('POSITIVE DIRECTION:', HDR.positive_direction))
    fid.close()

    with open(HDR.fname + '.rd3', 'wb') as fid:
        data.T.astype('short').tofile(fid)
    fid.close()


# DZT, Geophysical Survey Systems Inc. (GSSI) ==========================================================================
# Dzt is a binary file that consists of the file header with all the  important information (number of samples, traces,
# channels, etc.) followed by the data section.
# All the information is contained in this file except the TxRx distance (antenna separation). There is a possibility
# that the official GSSI software has stored this information and by using the antenna name presents the correct one.
# All the other software does not detect the TxRx distance.

elif gpr_format == 2:
    # Header structure
    HDR.fname = HDR.fname                                                       # File name
    HDR.tag = 255                                                               # Header = 255
    HDR.data_offset = 1024                                                      # Offset to data from the beginning of file
    HDR.num_samp = HDR.num_samp                                                 # Number of samples
    HDR.data_format = 16                                                        # Bits per data word (8, 16, 32)
    HDR.binary_offset = 32768                                                   # Binary offset, 8 bit = 128, 16 bit = 32768
    HDR.scans_per_second = 0                                                    # Scans per second
    HDR.scans_per_meter = 1 / HDR.trac_int                                      # Scans per meter
    HDR.meters_per_mark = 0                                                     # Meters per mark
    HDR.zero_time_adjustment = 0                                                # Time zero adjustment (ns)
    HDR.time_window = HDR.time_window                                           # Time window (with no corrections i.e zero time)
    HDR.scans_per_pass = 0                                                      # Scan per pass for 2D files

    HDR.createdate_sec = int(0 / 2)                                             # Structure, date created
    HDR.createdate_min = 0
    HDR.createdate_hour = 0
    HDR.createdate_day = 0
    HDR.createdate_month = 0
    HDR.createdate_year = 1980 - 1980

    date_time = datetime.datetime.now()
    HDR.modifydate_sec = int(np.round(date_time.second / 2, decimals=0))        # Structure, date modified
    HDR.modifydate_min = date_time.minute
    HDR.modifydate_hour = date_time.hour
    HDR.modifydate_day = date_time.day
    HDR.modifydate_month = date_time.month
    HDR.modifydate_year = date_time.year - 1980

    HDR.offset_to_range_gain = 0                                                # Offset to range gain
    HDR.size_of_range_gain = 0                                                  # Size of range gain
    HDR.offset_to_text = 0                                                      # Offset to text
    HDR.size_of_text = 0                                                        # Size of text
    HDR.offset_to_proc_his = 0                                                  # Offset to processing history
    HDR.size_of_proc_his = 0                                                    # Size of processing history
    HDR.num_channels = 1                                                        # Number of channels
    HDR.dielectric_constant = 8                                                 # Dielectric constant (8 is a random number)
    HDR.top_position = 0                                                        # Top position

    c = 299792458
    v = (c / np.sqrt(HDR.dielectric_constant)) * 10**-9
    HDR.range_depth = v * (HDR.time_window / 2)                                 # Range depth (m)

    if len(HDR.antenna) == 14:                                                  # Antenna name
        HDR.antenna = HDR.antenna
    elif len(HDR.antenna) < 14:
        HDR.antenna = HDR.antenna.ljust(14)
    elif len(HDR.antenna) > 14:
        HDR.antenna = HDR.antenna[:14]

    HDR.channel_mask = 0                                                        # Channel mask

    if len(HDR.fname) == 12:                                                    # Raw file name
        HDR.raw_file_name = HDR.fname
    elif len(HDR.fname) < 12:
        HDR.raw_file_name = HDR.fname.ljust(12)
    elif len(HDR.fname) > 12:
        HDR.raw_file_name = HDR.fname[:12]

    # DZT file
    with open(HDR.fname + '.dzt', 'wb') as fid:
        fid.write(struct.pack('<H', HDR.tag))
        fid.write(struct.pack('<H', HDR.data_offset))
        fid.write(struct.pack('<H', HDR.num_samp))
        fid.write(struct.pack('<H', HDR.data_format))
        fid.write(struct.pack('<H', HDR.binary_offset))
        fid.write(struct.pack('<f', HDR.scans_per_second))
        fid.write(struct.pack('<f', HDR.scans_per_meter))
        fid.write(struct.pack('<f', HDR.meters_per_mark))
        fid.write(struct.pack('<f', HDR.zero_time_adjustment))
        fid.write(struct.pack('<f', HDR.time_window))
        fid.write(struct.pack('<H', HDR.scans_per_pass))

        createdate_sec = Bits(uint=HDR.createdate_sec, length=5)
        createdate_min = Bits(uint=HDR.createdate_min, length=6)
        createdate_hour = Bits(uint=HDR.createdate_hour, length=5)
        createdate_day = Bits(uint=HDR.createdate_day, length=5)
        createdate_month = Bits(uint=HDR.createdate_month, length=4)
        createdate_year = Bits(uint=HDR.createdate_year, length=7)
        bits1 = Bits().join([createdate_year, createdate_month, createdate_day,
                             createdate_hour, createdate_min, createdate_sec])
        createdate = bits1.tobytes()
        fid.write(bitstruct.pack('>r32<', createdate))

        modifydate_sec = Bits(uint=HDR.modifydate_sec, length=5)
        modifydate_min = Bits(uint=HDR.modifydate_min, length=6)
        modifydate_hour = Bits(uint=HDR.modifydate_hour, length=5)
        modifydate_day = Bits(uint=HDR.modifydate_day, length=5)
        modifydate_month = Bits(uint=HDR.modifydate_month, length=4)
        modifydate_year = Bits(uint=HDR.modifydate_year, length=7)
        bits2 = Bits().join([modifydate_year, modifydate_month, modifydate_day,
                             modifydate_hour, modifydate_min, modifydate_sec])
        modifydate = bits2.tobytes()
        fid.write(bitstruct.pack('>r32<', modifydate))

        # fid.seek(40)
        fid.write(struct.pack('<H', HDR.offset_to_range_gain))
        fid.write(struct.pack('<H', HDR.size_of_range_gain))
        fid.write(struct.pack('<H', HDR.offset_to_text))
        fid.write(struct.pack('<H', HDR.size_of_text))
        fid.write(struct.pack('<H', HDR.offset_to_proc_his))
        fid.write(struct.pack('<H', HDR.size_of_proc_his))
        fid.write(struct.pack('<H', HDR.num_channels))
        fid.write(struct.pack('<f', HDR.dielectric_constant))
        fid.write(struct.pack('<f', HDR.top_position))
        fid.write(struct.pack('<f', HDR.range_depth))

        fid.seek(98)
        fid.write(HDR.antenna.encode('utf-8'))
        fid.write(struct.pack('<H', HDR.channel_mask))
        fid.write(HDR.raw_file_name.encode('utf-8'))

        fid.seek(HDR.data_offset, 0)
        data = data + 2**15
        data.T.astype('ushort').tofile(fid)
    fid.close()

# HD / DT1, Sensors & Software Inc. ====================================================================================
# Hd is the header file. In this file all the important information such as the number of samples, traces, stacks, etc.
# can be found.
# Dt1 is the data file written in binary form. This file contains as many records as there are traces. Each record
# consists of a header and a data section. This means that also in this file there are stored information such as the
#  number of samples, traces, etc.
elif gpr_format == 3:
    # Header structure of HD
    HDR.fname = HDR.fname                                                       # File name
    HDR.file_tag = 1234                                                         # File tag = 1234
    HDR.system = 'gprMax'                                                       # The system the data collected with

    date_time = datetime.datetime.now()
    HDR.date = (str(date_time.year) + '-' +
                str(date_time.month) + '-' + str(date_time.day))                # Date

    HDR.num_trac = HDR.num_trac                                                 # Number of traces
    HDR.num_samp = HDR.num_samp                                                 # Number of samples
    HDR.time_zero_point = 0                                                     # Time zero point
    HDR.time_window = HDR.time_window                                           # Total time window (ns)
    HDR.start_position = 0                                                      # Start position (m)
    HDR.final_position = (HDR.num_trac - 1) * HDR.trac_int                      # Stop position (m)
    HDR.trac_int = HDR.trac_int                                                 # Trace interval (m)
    HDR.pos_units = 'm'                                                         # Position units
    HDR.nominal_freq = HDR.centre_freq                                          # Nominal freq. / Centre freq. (MHz)
    HDR.ant_sep = HDR.ant_sep                                                   # Antenna seperation / Tx-Rx distance (m)
    HDR.pulser_voltage = 0                                                      # Pulser voltage (V)
    HDR.stacks = 1                                                              # Number of stacks
    HDR.survey_mode = 'Reflection'                                              # Survey mode
    HDR.odometer = 0                                                            # Odometer Cal (t/m)
    HDR.stacking_type = 'F1'                                                    # Stacking type
    HDR.dvl_serial = '0000-0000-0000'                                           # DVL serial
    HDR.console_serial = '000000000000'                                         # Console serial
    HDR.tx_serial = '0000-0000-0000'                                            # Transmitter serial
    HDR.rx_serial = '0000-0000-0000'                                            # Receiver Serial

    # Header structure of DT1
    HDR.num_each_trac = np.linspace(1, HDR.num_trac, HDR.num_trac)              # Number of each trace 1, 2, 3, ... num_trac
    HDR.position = np.linspace(0, (HDR.num_trac - 1) *
                               HDR.trac_int, HDR.num_trac)                      # Position of each trace (m)

    HDR.num_samp_each_trac = np.zeros([HDR.num_trac]) + HDR.num_samp            # Number of samples of each trace
    HDR.elevation = np.zeros([HDR.num_trac])                                    # Elevation / topography of each trace
    HDR.not_used1 = np.zeros([HDR.num_trac])                                    # Not used
    HDR.bytes = np.zeros([HDR.num_trac]) + 2                                    # Always 2 for Rev 3 firmware
    HDR.time_window_each_trac = np.zeros([HDR.num_trac]) + HDR.time_window      # Time window of each trace (ns)
    HDR.stacks_each_trac = np.ones([HDR.num_trac])                              # Number of stacks each trace
    HDR.not_used2 = np.zeros([HDR.num_trac])                                    # Not used
    HDR.rsv_gps_x = np.zeros([HDR.num_trac])                                    # Reserved for GPS X poition (double*8 number)
    HDR.rsv_gps_y = np.zeros([HDR.num_trac])                                    # Reserved for GPS Y position (double*8 number)
    HDR.rsv_gps_z = np.zeros([HDR.num_trac])                                    # Reserved for GPS Z position (double*8 number)
    HDR.rsv_rx_x = np.zeros([HDR.num_trac])                                     # Reserved for receiver x position
    HDR.rsv_rx_y = np.zeros([HDR.num_trac])                                     # Reserved for receiver y position
    HDR.rsv_rx_z = np.zeros([HDR.num_trac])                                     # Reserved for receiver z position
    HDR.rsv_tx_x = np.zeros([HDR.num_trac])                                     # Reserved for transmitter x position
    HDR.rsv_tx_y = np.zeros([HDR.num_trac])                                     # Reserved for transmitter y position
    HDR.rsv_tx_z = np.zeros([HDR.num_trac])                                     # Reserved for transmitter z position
    HDR.time_zero = np.zeros([HDR.num_trac])                                    # Time zero adjustment where: point(x) = point(x + adjustment)
    HDR.zero_flag = np.zeros([HDR.num_trac])                                    # 0 = data ok, 1 = zero data
    HDR.num_channels = np.zeros([HDR.num_trac])                                 # Number of channels
    HDR.time = np.zeros([HDR.num_trac])                                         # Time of day data collected in seconds past midnight
    HDR.comment_flag = np.zeros([HDR.num_trac])                                 # Comment flag
    HDR.comment = 'No comment'.ljust(24)                                        # Comment

    # HD file
    with open(HDR.fname + '.hd', 'w') as fid:
        fid.write('%i\n' % (HDR.file_tag))
        fid.write('%s%s\n' % ('Data Collected with ', HDR.system))
        fid.write('%s\n' % (HDR.date))
        fid.write('%s%i\n' % ('NUMBER OF TRACES   = ', HDR.num_trac))
        fid.write('%s%i\n' % ('NUMBER OF PTS/TRC  = ', HDR.num_samp))
        fid.write('%s%i\n' % ('TIMEZERO AT POINT  = ', HDR.time_zero_point))
        fid.write('%s%0.6f\n' % ('TOTAL TIME WINDOW  = ', HDR.time_window))
        fid.write('%s%0.6f\n' % ('STARTING POSITION  = ', HDR.start_position))
        fid.write('%s%0.6f\n' % ('FINAL POSITION     = ', HDR.final_position))
        fid.write('%s%0.6f\n' % ('STEP SIZE USED     = ', HDR.trac_int))
        fid.write('%s%s\n' % ('POSITION UNITS     = ', HDR.pos_units))
        fid.write('%s%0.6f\n' % ('NOMINAL FREQUENCY  = ', HDR.nominal_freq))
        fid.write('%s%0.6f\n' % ('ANTENNA SEPARATION = ', HDR.ant_sep))
        fid.write('%s%0.6f\n' % ('PULSER VOLTAGE (V) = ', HDR.pulser_voltage))
        fid.write('%s%i\n' % ('NUMBER OF STACKS   = ', HDR.stacks))
        fid.write('%s%s\n' % ('SURVEY MODE        = ', HDR.survey_mode))
        fid.write('%s%0.6f\n' % ('ODOMETER CAL (t/m) = ', HDR.odometer))
        fid.write('%s%s\n' % ('STACKING TYPE      = ', HDR.stacking_type))
        fid.write('%s%s\n' % ('DVL Serial#        = ', HDR.dvl_serial))
        fid.write('%s%s\n' % ('Console Serial#    = ', HDR.console_serial))
        fid.write('%s%s\n' % ('Transmitter Serial#= ', HDR.tx_serial))
        fid.write('%s%s\n' % ('Receiver Serial#   = ', HDR.rx_serial))
    fid.close()

    # DT1 file
    with open(HDR.fname + '.dt1', 'wb') as fid:
        for i in range(0, HDR.num_trac):
            fid.write(struct.pack('<f', HDR.num_each_trac[i]))
            fid.write(struct.pack('<f', HDR.position[i]))
            fid.write(struct.pack('<f', HDR.num_samp_each_trac[i]))
            fid.write(struct.pack('<f', HDR.elevation[i]))
            fid.write(struct.pack('<f', HDR.not_used1[i]))
            fid.write(struct.pack('<f', HDR.bytes[i]))
            fid.write(struct.pack('<f', HDR.time_window_each_trac[i]))
            fid.write(struct.pack('<f', HDR.stacks_each_trac[i]))
            fid.write(struct.pack('<f', HDR.not_used2[i]))
            fid.write(struct.pack('<d', HDR.rsv_gps_x[i]))
            fid.write(struct.pack('<d', HDR.rsv_gps_y[i]))
            fid.write(struct.pack('<d', HDR.rsv_gps_z[i]))
            fid.write(struct.pack('<f', HDR.rsv_rx_x[i]))
            fid.write(struct.pack('<f', HDR.rsv_rx_y[i]))
            fid.write(struct.pack('<f', HDR.rsv_rx_z[i]))
            fid.write(struct.pack('<f', HDR.rsv_tx_x[i]))
            fid.write(struct.pack('<f', HDR.rsv_tx_y[i]))
            fid.write(struct.pack('<f', HDR.rsv_tx_z[i]))
            fid.write(struct.pack('<f', HDR.time_zero[i]))
            fid.write(struct.pack('<f', HDR.zero_flag[i]))
            fid.write(struct.pack('<f', HDR.num_channels[i]))
            fid.write(struct.pack('<f', HDR.time[i]))
            fid.write(struct.pack('<f', HDR.comment_flag[i]))
            fid.write(HDR.comment.encode('utf-8'))
            data[:, i].astype('short').tofile(fid)
    fid.close()

# IPRH / IPRB, Impulse Radar ===========================================================================================
# IPRH is the header file. In this file is all the important information such as the number of samples, traces,
#  measurement intervals can be found.
# IPRB is the data file. This file contains only the data (amplitude values) in a binary form.
elif gpr_format == 4:
    # Header structure
    HDR.fname = HDR.fname                                                       # File name
    HDR.hdr_version = 20                                                        # Header version
    HDR.data_format = 16                                                        # Data format 16 or 32 bit

    date_time = datetime.datetime.now()
    HDR.date = (str(date_time.year) + '-' +
                str(date_time.month) + '-' + str(date_time.day))                # Date

    HDR.start_time = '00:00:00'                                                 # Measurement start time
    HDR.stop_time = '00:00:00'                                                  # Measurement end time
    HDR.antenna = (str(HDR.centre_freq) + ' MHz')                               # Antenna frequency (MHz)
    HDR.ant_sep = HDR.ant_sep                                                   # Antenna seperation / Tx-Rx distance (m)
    HDR.num_samp = HDR.num_samp                                                 # Number of samples
    HDR.signal_pos = 0                                                          # Signal position
    HDR.clipped_samps = 0                                                       # Clipped samples
    HDR.runs = 0                                                                # Number of runs
    HDR.stacks = 1                                                              # Maximum number of stacks
    HDR.auto_stacks = 1                                                         # Autostacks (1 = On)
    HDR.samp_freq = HDR.samp_freq                                               # Sampling frequency (MHz)
    HDR.time_window = HDR.time_window                                           # Total time window (ns)
    HDR.num_trac = HDR.num_trac                                                 # Number of traces
    HDR.trig_source = 'wheel'                                                   # Trig source (wheel or time)
    HDR.trac_int_sec = 0                                                        # Trace interval if trig source is time (sec)
    HDR.trac_int_met = HDR.trac_int                                             # Trace interval if trig source is wheel (m)
    HDR.user_trac_int = HDR.trac_int                                            # User defined trace interval if trig source is wheel (m)
    HDR.stop_pos = HDR.num_trac * HDR.trac_int                                  # Stop position (meters or seconds) -> num_trac * trac_int
    HDR.wheel_name = 'Cart'                                                     # Wheel name
    HDR.wheel_calibration = 0                                                   # Wheel calibration
    HDR.zero_lvl = 0                                                            # Zero level
    HDR.vel = 100                                                               # The soil velocity (Selected in field m/usec). 100 is a random number
    HDR.preprocessing = 'Unknown Preprocessing'                                 # Not in use
    HDR.comment = '----'                                                        # Not in use
    HDR.antenna_FW = '----'                                                     # Receiver firmware version
    HDR.antenna_HW = '----'                                                     # Not in use
    HDR.antenna_FPGA = '----'                                                   # Receiver FPGA version
    HDR.antenna_serial = '----'                                                 # Receiver serial number
    HDR.software_version = '----'                                               # Software version
    HDR.positioning = 0                                                         # Positioning: (0 = no, 1 = TS, 2 = GPS)
    HDR.num_channel = 1                                                         # Number of channels
    HDR.channel_config = 1                                                      # This channel configuration
    HDR.ch_x_offset = 0                                                         # Channel position relative to ext.positioning
    HDR.ch_y_offset = 0                                                         # Channel position relative to ext.positioning
    HDR.meas_direction = 1                                                      # Meas. direction forward or backward
    HDR.relative_direction = 0                                                  # Direction to RL start(clockwise 360)
    HDR.relative_distance = 0                                                   # Distance from RL start to cross section
    HDR.relative_start = 0                                                      # DIstance from profile start to cross section

    with open(HDR.fname + '.iprh', 'w') as fid:
        fid.write('%s%i\n' % ('HEADER VERSION: ', HDR.hdr_version))
        fid.write('%s%i\n' % ('DATA VERSION: ', HDR.data_format))
        fid.write('%s%s\n' % ('DATE: ', HDR.date))
        fid.write('%s%s\n' % ('START TIME: ', HDR.start_time))
        fid.write('%s%s\n' % ('STOP TIME: ', HDR.stop_time))
        fid.write('%s%s\n' % ('ANTENNA: ', HDR.antenna))
        fid.write('%s%0.6f\n' % ('ANTENNA SEPARATION: ', HDR.ant_sep))
        fid.write('%s%i\n' % ('SAMPLES: ', HDR.num_samp))
        fid.write('%s%0.6f\n' % ('SIGNAL POSITION: ', HDR.signal_pos))
        fid.write('%s%i\n' % ('CLIPPED SAMPLES: ', HDR.clipped_samps))
        fid.write('%s%i\n' % ('RUNS: ', HDR.runs))
        fid.write('%s%i\n' % ('MAX STACKS: ', HDR.stacks))
        fid.write('%s%i\n' % ('AUTOSTACKS: ', HDR.auto_stacks))
        fid.write('%s%0.6f\n' % ('FREQUENCY: ', HDR.samp_freq))
        fid.write('%s%0.6f\n' % ('TIMEWINDOW: ', HDR.time_window))
        fid.write('%s%i\n' % ('LAST TRACE: ', HDR.num_trac))
        fid.write('%s%s\n' % ('TRIG SOURCE: ', HDR.trig_source))
        fid.write('%s%0.6f\n' % ('TIME INTERVAL: ', HDR.trac_int_sec))
        fid.write('%s%0.6f\n' % ('DISTANCE INTERVAL: ', HDR.trac_int_met))
        fid.write('%s%0.6f\n' % ('USER DISTANCE INTERVAL: ', HDR.user_trac_int))
        fid.write('%s%0.6f\n' % ('STOP POSITION: ', HDR.stop_pos))
        fid.write('%s%s\n' % ('WHEEL NAME: ', HDR.wheel_name))
        fid.write('%s%0.6f\n' % ('WHEEL CALIBRATION: ', HDR.wheel_calibration))
        fid.write('%s%i\n' % ('ZERO LEVEL: ', HDR.zero_lvl))
        fid.write('%s%i\n' % ('SOIL VELOCITY: ', HDR.vel))
        fid.write('%s%s\n' % ('PREPROCESSING: ', HDR.preprocessing))
        fid.write('%s%s\n' % ('OPERATOR COMMENT: ', HDR.comment))
        fid.write('%s%s\n' % ('ANTENNA F/W: ', HDR.antenna_FW))
        fid.write('%s%s\n' % ('ANTENNA H/W: ', HDR.antenna_HW))
        fid.write('%s%s\n' % ('ANTENNA FPGA: ', HDR.antenna_FPGA))
        fid.write('%s%s\n' % ('ANTENNA SERIAL: ', HDR.antenna_serial))
        fid.write('%s%s\n' % ('SOFTWARE VERSION: ', HDR.software_version))
        fid.write('%s%i\n' % ('POSITIONING: ', HDR.positioning))
        fid.write('%s%i\n' % ('CHANNELS: ', HDR.num_channel))
        fid.write('%s%i\n' % ('CHANNEL CONFIGURATION: ', HDR.channel_config))
        fid.write('%s%0.6f\n' % ('CH_X_OFFSET: ', HDR.ch_x_offset))
        fid.write('%s%0.6f\n' % ('CH_Y_OFFSET: ', HDR.ch_y_offset))
        fid.write('%s%i\n' % ('MEASUREMENT DIRECTION: ', HDR.meas_direction))
        fid.write('%s%i\n' % ('RELATIVE DIRECTION: ', HDR.relative_direction))
        fid.write('%s%0.6f\n' % ('RELATIVE DISTANCE: ', HDR.relative_distance))
        fid.write('%s%0.6f\n' % ('RELATIVE START: ', HDR.relative_start))
    fid.close()

    with open(HDR.fname + '.iprb', 'wb') as fid:
        data.T.astype('short').tofile(fid)
    fid.close()

sys.exit(0)
