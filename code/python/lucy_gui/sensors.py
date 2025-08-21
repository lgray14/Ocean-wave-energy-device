import serial
import numpy as np

from time import process_time,sleep,time

from threading import Thread


class SerialClosedError(Exception):
    pass


class SerialReadError(Exception):
    pass


nbits_analog = {'arduino':10,'bbb':12}


class SensorReader(object):
    average_over = 100
    read_freq = 2.0   # in

    max_values_buffer_size = 1000

    def __init__(self,dev,datatypes,**kwargs):

        self.name = kwargs.pop('name','')
        self.delimiter = kwargs.pop('delimiter',None)
        self.board = kwargs.pop('board','arduino')
        self.ADC_resolution = 2**nbits_analog[self.board]
        import sys 
        # print("~~~~~~~DIRECTORY~~~~~~~~~")
        # print(dir(serial))
        self.serial = serial.Serial(dev,**kwargs)
        self.port = self.serial.port
        self.datatypes = datatypes
        self.buffer = []
        self.values_buffer = []
        self.read_times = []

        # waste a read
        try:
            self.read()
        except Exception:
            pass

        self.stop_threads = True

        # self.start_thread()

        return

    def _read(self):
        while self.serial.is_open:
            try:
                msg = self.serial.readline().decode()
                # msg = ser.read_until().strip().decode()
            except Exception:
                pass
            else:
                break
        return msg

    def threaded_read(self):
        while self.serial.is_open and not self.stop_threads:
            msg = self._read()
            # self.buffer.append(msg)
            # self.read_times.append(time())

            self.buffer.append((time(),msg))


            # sleep(0.01)
            # print('threaded_read: len(self.buffer) = ',len(self.buffer))
            # print(msg,end='')
        return

    def process_buffer(self):
        # print(f'in process_buffer: id(self) = {id(self)}, len(self.buffer) = {len(self.buffer)}')

        # print('process_buffer 0: ',len(self.buffer),len(self.values_buffer),len(self.read_times))
        while len(self.buffer) > 0:
            # print('process_buffer: len(self.buffer) = ',len(self.buffer))

            # vals = self.process_oldest()
            t,vals = self.process_oldest()

            # print('process_buffer 1: ',len(self.buffer),len(self.values_buffer),len(self.read_times))
            self.values_buffer.append(vals)
            self.read_times.append(t)
            # print('process_buffer 2: ',len(self.buffer),len(self.values_buffer),len(self.read_times))

            # s = self.buffer.pop(0)
            # ss = s.split(self.delimiter)
            # try:
            #     vals = [typ(val) for typ,val in zip(self.datatypes,ss)]
            # except Exception:
            #     pass
            # else:
            #     self.values_buffer.append(vals)

            if len(self.values_buffer) > self.max_values_buffer_size:
                self.values_buffer.pop(0)
                self.read_times.pop(0)

            # print('process_buffer 3: ',len(self.buffer),len(self.values_buffer),len(self.read_times))

            # print('process_buffer: vals = ',self.values_buffer)
            # print('process_buffer: time = ',self.read_times)

            if len(self.read_times) != len(self.values_buffer):
                print(self.read_times)
                raise ValueError('Buffer sizes do not match!',
                                 len(self.read_times),len(self.values_buffer))

        return

    def process_oldest(self):
        # print(f'process_oldest: len(self.buffer) = {len(self.buffer)}')
        while True:
            if len(self.buffer) == 0:
                return

            t,s = self.buffer.pop(0)
            # s = self.buffer.pop(0)
            ss = s.split(self.delimiter)
            try:
                # vals = [typ(val) for typ,val in zip(self.datatypes,ss)]
                vals = []
                for typ,strval in zip(self.datatypes,ss):
                    val = typ(strval)
                    if typ == int:
                        val /= self.ADC_resolution - 1
                    vals.append(val)
            except Exception:
                # pass
                print('Error in process_oldest: ',s)
            else:
                # return vals
                return t,vals

    def read_raw(self):
        if not self.serial.is_open:
            print('Port not open!',self.serial)
            raise SerialClosedError
        red = False
        max_iter = 3
        i = 0
        while not red:
            # try:
            #     s = self.serial.readline().decode('utf-8').strip()
            # except Exception:
            #     i += 1
            # else:
            #     red = True
            # if i > max_iter:
            #     raise SerialReadError(f'Could not read serial in {i} attempts')

            try:
                s = self.serial.readline().decode('utf-8').strip()
            except Exception:
                i += 1
                continue

            ss = s.split(self.delimiter)
            try:
                vals = [typ(val) for typ,val in zip(self.datatypes,ss)]
            except Exception:
                i += 1
            else:
                red = True

            if i > max_iter:
                raise SerialReadError(f'Could not read serial in {i} attempts')

        return tuple(vals)

    def read(self,raw=False):
        vals = self.read_raw()
        if raw:
            return vals

        cvals = []
        for typ,val in zip(self.datatypes,vals):
            if typ == int:
                val /= self.ADC_resolution - 1
            cvals.append(val)
        return tuple(cvals)

    def continuous_read(self):
        # self.vals = []
        vals = []
        self.vals = vals
        tstart = process_time()
        while not self.stop_thread:
            v = self.read()
            # vals.append(v)
            vals.append(v)
            if len(vals) > self.average_over:
                vals.pop(0)
            # if len(vals) < self.average_over:
            #     vals.append(v)
            # else:
            #     vals = np.roll(vals,shift=-1)
            #     vals[-1] = v
            t = process_time()

            if t-tstart >= 1/self.read_freq:
                vmean = np.mean(vals,axis=0)
                # print(vmean,' {}   {:7.3f}'.format(len(vals),(process_time()-tstart)*self.read_freq))
                self.current_val = vmean
                vals = []
                self.vals = vals
                tstart = process_time()

        return

    def start_threads(self):
        self.stop_threads = False

        # self.thread = Thread(target=self.continuous_read,name=self.name)
        self.thread_read = Thread(target=self.threaded_read,name=self.name)
        # self.thread_process = Thread(target=self.process_buffer,name=self.name)

        self.thread_read.start()
        # self.thread_process.start()
        # self.thread.join()  # make the main program wait until thread is done
        return

    def stop(self):
        self.stop_threads = True
        return

    def disconnect(self):
        self.stop()
        self.serial.close()
        return

    def connect(self):
        if not self.serial.is_open:
            self.serial.open()
        return

    def start(self):
        self.connect()
        self.start_threads()
        return

    # def stop_thread(self):
    #     self.thread.stop()
    #     return


def threaded_read(reader):
    reader.stop_threads = False
    thread_read = Thread(target=reader.threaded_read,name=reader.name)
    # thread_process = Thread(target=reader.process_buffer,name=reader.name)

    thread_read.start()
    return thread_read
