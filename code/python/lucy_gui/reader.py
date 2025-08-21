import numpy as np

# import sensors as SS

from time import time  #,sleep


class StatusError(Exception):
    pass


class CoreReader(object):
    freq = 20.

    def __init__(self,shape,zero_offset=None,scale=None,**kwargs):
        self.values = np.zeros(shape)
        self.values[:] = np.nan
        self.shape = self.values.shape

        self.n = self.shape[0]
        self.read_times = np.zeros(self.n)
        self.read_times[:] = np.nan

        self._isready = False

        self.stop_thread = False

        self.zero_offset = zero_offset
        self.scale = scale

        # if np.shape(self.zero_offset)[-1] != self.shape[-1]:
        #     raise ValueError('Shapes do not match!')

        # if np.shape(self.scale)[-1] != self.shape[-1]:
        #     raise ValueError('Shapes do not match!')
        return

    def read(self):
        pass

    def continuous_read(self):
        if not self.is_ready:
            raise StatusError('Reader not ready')
        while True:
            # Reads at least once
            self.read()
            # sleep(1/self.freq)
            if self.stop_thread:
                break
        return

    # def start(self):

    @property
    def is_ready(self):
        return self._isready

    def connect(self):
        self._isready = True
        return

    def disconnect(self):
        self._isready = False
        return

    def insert_read(self,vals_raw,time=None):
        self.values[...] = np.roll(self.values,-1,axis=0)

        if self.values.ndim == 1:
            vals = vals_raw[0]
        else:
            m = self.shape[-1]
            vals = vals_raw
            if m == 0:
                vals = vals[0]

        # if self.zero_offset is not None:
        #     vals -= self.zero_offset
        # if self.scale is not None:
        #     vals *= self.scale

        self.values[-1] = vals  # np.random.rand(self.shape[1])[0]
        # print('insert_read',self.values.shape)

        if time is not None:
            self.read_times[:] = np.roll(self.read_times,-1)
            self.read_times[-1] = time
        return

    def set_zero_offsets(self,zero_offsets):
        self.zero_offset = zero_offsets
        return

    def set_zero_offset_i(self,i,zero_offset):
        if self.zero_offset is None:
            self.zero_offset = np.zeros(self.shape[-1])

        self.zero_offset[i] = zero_offset
        # print('set_zero_offset_i',self.zero_offset)
        return

    def set_as_zero_offsets(self):
        # uses the current reading as zero
        zero_offsets = self.values[-1]
        self.set_zero_offsets(zero_offsets)
        return

    def set_as_zero_offset_i(self,i):
        # uses the current reading as zero
        # print(self.values.shape)
        # print('set_as_zero_offset_i: ',self.zero_offset,i,self.shape,self.values.shape)
        # print('set_as_zero_offset_i: ',f'Reader ID: {hex(id(self))}')
        zero_offset = self.values[-1,i]
        # print('set_as_zero_offset_i: ',zero_offset)
        self.set_zero_offset_i(i,zero_offset)
        return

    def set_scales(self,scales):
        self.scale = scales
        return

    def set_scale_i(self,i,scale):
        # print('set_scale_i: ',self.scale,i,scale,self.shape,self.values.shape)
        # print('set_scale_i: ',f'Reader ID: {hex(id(self))}')
        if self.scale is None:
            self.scale = np.ones(self.shape[-1])

        self.scale[i] = scale
        return

    @property
    def scaled_values(self):
        svals = self.values
        if self.zero_offset is not None:
            svals = svals - self.zero_offset  # to ensure the array is a copy
        if self.scale is not None:
            svals = svals*self.scale
        return svals


class FakeReader(CoreReader):
    def read(self):
        if self.values.ndim == 1:
            vals = np.random.rand(1)
        else:
            m = self.shape[-1]
            vals = np.random.rand(m)
        self.insert_read(vals,time=time())
        # self.values[...] = np.roll(self.values,-1,axis=0)

        # if self.values.ndim == 1:
        #     vals = np.random.rand(1)[0]
        # else:
        #     m = self.shape[-1]
        #     vals = np.random.rand(m)
        #     if m == 0:
        #         vals = vals[0]
        # self.values[-1] = vals  # np.random.rand(self.shape[1])[0]
        return


# roller_diameter = 1.5*2.54  # diameter in cm
# PTO_roller_diameter = 1.7*2.54


# pot1_distance_range = 10*roller_diameter*np.pi
# pot2_distance_range = 10*PTO_roller_diameter*np.pi


class DataReader(CoreReader):
    # def __init__(self,shape,reader,devices,zero_offset=None,scale=None):
    def __init__(self,shape,reader,zero_offset=None,scale=None):
        super().__init__(shape,zero_offset=zero_offset,scale=scale)
        self.reader = reader
        self.reader.serial.flush()

    # def read(self):
    #     # self.values[:] = np.roll(self.values,-1)
    #     red = False
    #     while not red:
    #         try:
    #             vals = self.reader.read()  #[0]
    #         except Exception:
    #             pass
    #         else:
    #             red = True
    #             self.insert_read(vals)
    #     # print(vals)
    #     return

    def read(self):
        # print('asdfasdfafdsasdfasdfasdafsd')
        self.reader.process_buffer()
        # vals = np.array(self.reader.values_buffer)
        # for vals in self.reader.values_buffer:
        #     self.insert_read(vals)
        for vals,t in zip(self.reader.values_buffer,self.reader.read_times):
            # print(vals,t)
            self.insert_read(vals,time=t)
        return

    @property
    def is_ready(self):
        return self.reader.serial.is_open

    def connect(self):
        print('IN DataReader CONNECT')
        self.reader.start()
        # self.reader.connect()
        super().connect()
        return

    def disconnect(self):
        self.reader.disconnect()
        super().disconnect()
        return
