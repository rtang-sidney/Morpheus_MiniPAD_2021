import numpy as np


START_CURRENT = 0.0
END_CURRENT = 3.0
STEP_SIZE = 0.1
STEPS = (END_CURRENT - START_CURRENT) / STEP_SIZE + 1
STEPS = round(STEPS)

def __init__(self, scan_start, scan_end):
    self.START_INDEX = scan_start
    self.END_INDEX = scan_end

    self.plot_x = np.linspace(start=self.START_CURRENT, stop=self.END_CURRENT, num=self.STEPS)
    self.plot_y = self.plot_x
    self.plot_xx, self.plot_yy = np.meshgrid(self.plot_x, self.plot_y)
    self.plot_zz = np.empty_like(self.plot_xx)

def get_data(self, write_file_name):
    if isinstance(write_file_name, str):
        write_file = open('{:s}.txt'.format(write_file_name), 'w')
    else:
        raise RuntimeError('Write file name has to be a string. Given name {} invalid.'.format(write_file_name))

    sth = -70  # starting point
    for file in range(self.START_INDEX, self.END_INDEX + 1):
        data_file = DataFile(file_index=file)
        data_file.property_sort()
        if data_file.DCT in data_file.scan_name and data_file.x_data.shape[0] > 5:
            pass
        else:
            continue

        if data_file.properties_sorted[data_file.SAMPLE_TABLE_ANGLE] - sth > self.ANGLE_ERROR:
            # self._plot_save_picture(sth=sth)
            sth = round(data_file.properties_sorted[data_file.SAMPLE_TABLE_ANGLE])
        else:
            pass

        # this is to collect the information from the raw data files into the info file
        for i in range(data_file.x_data.shape[0]):
            write_file.write(
                '{:d} {:f} {:f} {:f} {:d} \n'.format(file, sth, data_file.x_data[i],
                                                     data_file.properties_sorted[data_file.DCT6],
                                                     data_file.real_counts[i]))
        # self._plot_points(data_file.x_data, data_file.properties_sorted[data_file.DCT6], data_file.real_counts)

    # self._plot_save_picture(sth=sth)
    write_file.close()

def data_get_fit_plot(self):
    fig, axs = plt.subplots(3, 3, sharex="all", sharey="all", figsize=(10, 10))
    ang1 = 0
    ang2 = 45
    i = 0
    # sth = -70  # starting point
    for file in range(self.START_INDEX, self.END_INDEX + 1):
        data_file = DataFile(file_index=file)
        data_file.property_sort()

        # only continue with the data analysis if there are more than 5 data points
        if data_file.DCT in data_file.scan_name and data_file.x_data.shape[0] > 5:
            pass
        else:
            continue

        # calculate intensity for the plot points
        self._plot_fit_points(data_file.x_data, data_file.properties_sorted[data_file.DCT6], data_file.real_counts)

        # if data_file.properties_sorted[data_file.SAMPLE_TABLE_ANGLE] - sth > self.ANGLE_ERROR:
        #     self._plot_save_picture(sth=sth, fit_type='fit')
        #     sth = round(data_file.properties_sorted[data_file.SAMPLE_TABLE_ANGLE])
        # else:
        #     pass
        sth = data_file.properties_sorted[data_file.SAMPLE_TABLE_ANGLE]
        sth = int(round(sth))
        print(sth)
        if ang1 <= sth <= ang2 and i < 9:
            axi, axj = int(i / 3), i - int(i / 3) * 3
            for j in range(self.plot_zz.shape[0]):
                self.plot_zz[j, :] = counts2polarisation(self.plot_zz[j, :])[0]
            cnt = axs[axi, axj].contourf(self.plot_xx, self.plot_yy, self.plot_zz, cmap='coolwarm')
            cbar = fig.colorbar(cnt)
            axs[axi, axj].set_title(r'Polarisation scan at $\varphi=${:d}°'.format(sth))
            # self.plot_zz = np.empty_like(self.plot_xx)
            print('Picture{}, sth {}.'.format(i, sth))
            i += 1
        else:
            continue

    # self._plot_save_picture(sth=sth, fit_type='fit')
    fig.savefig("FinalScans.png")
