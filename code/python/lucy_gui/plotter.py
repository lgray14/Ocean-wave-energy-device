from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as FigureCanvas
# from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib.pyplot import twinx

import numpy as np


# from matplotlib.backend_bases import key_press_handler
# from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
#                                                NavigationToolbar2Tk)
# from matplotlib.figure import Figure

class MovingWindow(object):
    pass


# def plot_d

line_colors = {1: 'tab:blue', 2: 'tab:orange'}


class MplCanvas(FigureCanvas):

    # def __init__(self, parent=None, width=5, height=4, dpi=100):
    def __init__(self, parent=None,**kwargs):
        self.time_window = kwargs.pop('time_window',20.0)

        figkw = {'layout':'constrained'}
        figkw.update(kwargs)
        fig,axs = plt.subplots(1,2,**figkw)
        # fig = Figure(figsize=(width, height), dpi=dpi)
        # self.axes = fig.add_subplot(111)
        self.figure = fig
        self.axes = axs
        super(MplCanvas, self).__init__(fig)

        ax = axs[0]
        ax.set_title('Position')

        ax = axs[1]
        ax.set_title('Current')

        self.line_pot1 = None
        self.line_power = None
        self.line_current = None

        # self.t = np.arange(0, 3, .01)

        return

    @property
    def lines(self):
        return self.line_pot1,self.line_power,self.line_current

    def setup(self):
        self.initialize_position_plot()
        self.initialize_current_plot()

        # self.figure.tight_layout()
        self.update_xaxes()
        return

    def update_xaxes(self):
        for ax in self.axes:
            ax.set_xlim(self.time_window,0.0)
            ax.set_xlabel("Time before now [s]")
        return

    def initialize_position_plot(self,init=True):
        ax = self.axes[0]
        # self.line_pot1, = ax.plot(self.t, 2 * np.sin(2 * np.pi * self.t))
        # self.line_pot1 = plt.Line2D((0,0),(np.nan,)*2,c='tab:blue',label='Pot1')
        # self.line_power = plt.Line2D((0,0),(np.nan,)*2,c='tab:orange',lw=0.75,
        #                             label='Pot2')
        self._initialize_position_plot_line(1,init=init)
        self._initialize_position_plot_line(2,init=init)
        # ax.set_xlabel("time [s]")
        ax.set_ylabel(r"$z(t)$ [cm]")

        ax.set_ylim(-1,1)
        # ax.set_xlim(0,3)

        # ax.legend(handles=(self.line_pot1,self.line_power),ncol=2,
        #           loc='lower center',bbox_to_anchor=(0.5,1.0))
        return

    def _initialize_position_plot_line(self,iline,init=True):
        # line = plt.Line2D((0,0),(np.nan,)*2,c=line_colors[iline],
        #                   label=f'Pot{iline}')
        ax = self.axes[0]
        # print("iline variable:", iline)
        if init:
            line, = ax.plot((0,0),(np.nan,)*2,c=line_colors[iline],
                            label=f'Pot{iline}')
            setattr(self,f'line_pot{iline}',line)
        else:
            print('clearing lines')
            # getattr(self,f'line_pot{iline}').set_data((0,0),(np.nan,)*2)
            self.line_pot1.set_data((0,0),(np.nan,)*2)
            # self.line_power.set_data((0,0),(np.nan,)*2)
        return

    def initialize_current_plot(self):
        ax = self.axes[1]
        # self.line_current = plt.Line2D((0,0),(np.nan,)*2,c='tab:blue')
        self.line_current, = ax.plot((0,0),(np.nan,)*2,c='tab:red', label="Current")
        # self.line_current, = ax.plot(self.t, np.zeros_like(self.t))
        # ax.set_ylim(0,1024)
        ax.set_ylim(-10,10)
        ax.set_ylabel(r'$Current$ [mA]')

        # Create the twin axis and power line only once
        if not hasattr(self, 'ax_current_twin') or self.ax_current_twin is None:
            self.ax_current_twin = ax.twinx()
            self.ax_current_twin.set_ylabel(r'$Power$ [W]')
            self.ax_current_twin.set_ylim(0, 1)  # Set a default y-limits for power
            self.line_power, = self.ax_current_twin.plot((0,0), (np.nan,)*2, c='tab:orange', label='Power')
            # print("Initialized power line and twin axis")
        else:
            # Clear the line if re-initializing
            self.line_power.set_data((0,0), (np.nan,)*2)
            # print("Re-initialized power line and twin axis")
        # Add legend for clarity
        ax.legend(loc='upper left')
        self.ax_current_twin.legend(loc='upper right')
        return

    # def update_frequency(self,new_val):
    #     # retrieve frequency
    #     f = float(new_val)

    #     # update data
    #     y = 2 * np.sin(2 * np.pi * f * self.t)
    #     self.line.set_data(self.t, y)

    #     # required to update canvas and attached toolbar!
    #     self.draw()

    # def update_live(self,values,plot_bools):
    def update_live(self, t_raw, power_data, values, plot_bools, raw_data=False):
        if np.shape(values)[-1] != len(plot_bools):
            raise ValueError('Shapes not matching')

        update_position = plot_bools[0] or plot_bools[1]
        t = -(t_raw - t_raw[-1])

        if update_position:
            ax = self.axes[0]
            if raw_data:
                ax.set_ylim(0, 1)
            for i in [0, 1]:
                line = getattr(self, f'line_pot{i+1}') # position
                if plot_bools[i]:
                    vals = values[:, i]
                    line.set_data(t, vals)
                    adjust_ylim(ax, vals)
                else:
                    line.set_data((0, 0), (np.nan,) * 2)

        ax = self.axes[1]
        line = self.line_current # current
        if plot_bools[2]:
            vals = values[:, -1]
            t = -(t_raw - t_raw[-1])
            line.set_data(t, vals)
            adjust_ylim(ax, vals)
        else:
            line.set_data((0, 0), (np.nan,) * 2)

        # power
        power_times, power_vals = power_data
        # Debug: print power data received by plotter
        # print(f"[PLOTTER DEBUG] Received power_times={power_times[-5:] if len(power_times)>0 else power_times}, power_vals={power_vals[-5:] if len(power_vals)>0 else power_vals}")
        # Use the pre-created twin axis and line
        if hasattr(self, 'line_power') and self.line_power is not None:
            # print("got here")
            if len(power_times) > 0:
                # Align power x-axis to match the current plot's time window
                if len(t_raw) > 0:
                    t0 = t_raw[-1]
                    power_x = -(power_times - t0)
                    self.line_power.set_data(power_x, power_vals)
                    # print("Power times (aligned)", power_x)
                else:
                    self.line_power.set_data(power_times, power_vals)
                self.line_power.set_visible(True)
                # print('set data to t and power_vals')

                # Optionally, keep the power axis y-limits fixed or auto-adjust
                adjust_ylim(self.ax_current_twin, power_vals)
                # # Also, set the x-limits to match the current plot
                # self.ax_current_twin.set_xlim(self.time_window, 0.0)
            else:
                self.line_power.set_data((0, 0), (np.nan,) * 2)
                self.line_power.set_visible(False)
                # print('cleared data')
        # required to update canvas and attached toolbar!
        self.draw()


def adjust_ylim(ax,vals,c=1.1,onlydata=True):

    ylims = ax.get_ylim()

    # print('ylmi: ',ylims,np.nanmin(vals),np.nanmax(vals))
    if np.nanmin(vals) < ylims[0] or np.nanmax(vals) > ylims[1]:
        if onlydata:
            ax.set_ylim(c*np.array([np.nanmin(vals),np.nanmax(vals)])
                        )
        else:
            ax.set_ylim(c*np.array([min(np.nanmin(vals),ylims[0]),
                                    max(np.nanmax(vals),ylims[1])])
                        )
    return
