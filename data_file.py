import numpy as np
import re


class DataFile:
    # PATTERN_HEADER = re.compile(r"\**\sDATA\s\**")
    PATTERN_POINTS = re.compile(r"([0-9]*)\sPoints")
    PATTERN_SCAN = re.compile(r"Scanning Variables:\s(tti[1-2]_out[1-2]),\sSteps:\s([1-9]*.*[1-9]*)")

    DATA_HEADER = "**************************** DATA ******************************************\n"
    SCANNING_VARIABLE_KEY = "Scanning Variables"
    NORMED_DETECTOR_EFFICIENCY = 3.24e4
    NORMED_DETECTOR_WAVELENGTH = 1.8e-10  # m
    INCIDENT_WAVE_VECTOR = 1.4e10  # m-1

    FILE_PREFIX = 'morpheus2021n'
    FILE_FORMAT = '.dat'

    PARAMETER_SUFFIX = '_value'
    DATA_START = 'Scan data'
    COUNT_TIME = 'det_preset'
    SCAN_INFO = 'info'

    PLOT_Y_LABEL = 'Normalised counts'

    # Parameter names of the current

    DCT = 'dct'
    COIL_O2 = 'tti1_out1'
    COIL_I1 = 'tti1_out2'
    COIL_I2 = 'tti2_out1'
    COIL_O1 = 'tti2_out2'
    CHANNELS = [COIL_I1, COIL_O1, COIL_I2, COIL_O2]

    I1_POSITION = 'upstream, inner'
    O1_POSITION = 'upstream, outer'
    I2_POSITION = 'downstream, inner'
    O2_POSITION = 'downstream, outer'
    POSITIONS = [I1_POSITION, O1_POSITION, I2_POSITION, O2_POSITION]

    COILs_POSITIONs = dict(zip(CHANNELS, POSITIONS))

    SCATTERING_ANGLE = 'stt'

    # Parameter names of the sample table
    SAMPLE_TABLE_ANGLE = 'sth_st'
    SAMPLE_TABLE_POSITION = ['stx', 'sty', 'stz']
    SAMPLE_TABLE_TILT = ['sgx', 'sgy']
    SAMPLE_TABLE = [SAMPLE_TABLE_ANGLE] + SAMPLE_TABLE_POSITION + SAMPLE_TABLE_TILT

    # Including two pieces of information, i.e. width and height
    SLIT1 = 'ss1'  # (centre_x, centre_y) width x height mm
    SLIT2 = 'ss2'  # but we do not need it at the moment
    SLITs = [SLIT1, SLIT2]

    SLIT_INFO_1 = 'centre_x'
    SLIT_INFO_2 = 'centre_y'
    SLIT_INFO_3 = 'width'
    SLIT_INFO_4 = 'height'
    SLIT_INFOs = [SLIT_INFO_1, SLIT_INFO_2, SLIT_INFO_3, SLIT_INFO_4]

    SLITS_INFOS = []
    for slit in SLITs:
        for slit_info in SLIT_INFOs:
            SLITS_INFOS.append(str(slit) + '_' + str(slit_info))

    # With a form: xx_value : $float unit
    NORMAL_FLOAT_VALUES = CHANNELS + SAMPLE_TABLE + [SCATTERING_ANGLE]

    INFOS_NEEDED = [SCAN_INFO] + NORMAL_FLOAT_VALUES + SLITS_INFOS

    def __init__(self, file_index):
        self.complete_scan = True
        self.filename = self.FILE_PREFIX + "{:0>6d}".format(file_index) + self.FILE_FORMAT
        # print(self.filename)
        f = open(self.filename, 'r')
        lines = f.readlines()

        header_index = lines.index(self.DATA_HEADER)
        header_number = header_index + 3 + 1

        point_number = 0
        self.scan_name = None
        self.scan_step = 0
        for line in lines:
            if re.match(self.PATTERN_POINTS, line):
                point_number = re.search(self.PATTERN_POINTS, line).groups()[0]
                point_number = int(point_number)
                # print(point_number)
            if self.SCANNING_VARIABLE_KEY in line:
                self.scan_name, self.scan_step = re.search(self.PATTERN_SCAN, line).groups()
                self.scan_step = float(self.scan_step)
        if self.scan_name:
            # print(self.scan_name)
            pass
        else:
            raise ValueError("Failed to find the scan name.")

        if point_number > 0:
            try:
                data = np.loadtxt(self.filename, skiprows=header_number, max_rows=point_number)
                self.scan_x, self.scan_count = data[:, 1], data[:, 2]
                data_length = self.scan_x.shape[0]
                for i in range(data_length - 1):
                    if self.scan_x[i + 1] < self.scan_x[i]:
                        self.scan_x[i + 1] = self.scan_x[i] + self.scan_step
            except ValueError:
                self.complete_scan = False
        else:
            raise ValueError("Failed to find the number of scanned points.")

        # print(self.scan_x, self.scan_count)
        # print(lines)
        # print()
        #
        # # reads only the data part of a file
        # # x-axis	timer	mon1	mon2	ctr1    ctr2
        # try:
        #     data = np.loadtxt(self.filename, comments='#')
        #
        #     if data.ndim == 2:
        #         self.data_2d_array = True
        #         self.x_data = data[:, 0]
        #         counts = data[:, -2]
        #         monitor = data[:, -3]
        #         self.real_counts = counts / (self._normalised_counts() * monitor)
        #         self.real_counts = np.int_(self.real_counts)
        #         self.intensity_error = np.sqrt(counts + counts ** 2 / monitor) / (self._normalised_counts() * monitor)
        #         self._raw_data_correction()
        #         self.scan_name = ''
        #         self.x_label = ''
        #         self.y_label = self.PLOT_Y_LABEL
        #
        #         self.__properties = {}
        #         self.properties_sorted = {}
        #     else:
        #         self.data_2d_array = False
        # except UserWarning:
        #     self.data_2d_array = False

    def _raw_data_correction(self):
        error = 1e-3
        data_length = self.x_data.shape[0]
        delete_index = []
        for i in range(data_length):
            # one point was shifted to the next
            if 0 < i < data_length - 2 and self.x_data[i] > self.x_data[i - 1] and abs(
                    self.x_data[i] - self.x_data[i + 1]) < error:
                self.x_data[i] = (self.x_data[i - 1] + self.x_data[i + 1]) / 2.0
            # zero count
            if int(self.real_counts[i]) == 0:
                delete_index.append(i)

        self.x_data = np.delete(self.x_data, delete_index)
        self.real_counts = np.delete(self.real_counts, delete_index)
        # self.monitor = np.delete(self.monitor, delete_index)
        if self.x_data.shape[0] == 0:
            self.data_2d_array = False
        else:
            # first point was shifted to another position
            if self.x_data[0] > self.x_data[1]:
                self.x_data[0] = 0

    @staticmethod
    def _unit_correction(unit):
        if unit == 'A':
            pass
        else:
            unit = 'A'
        return unit

    def _normalised_counts(self):
        wavelength = 2 * np.pi / float(self.INCIDENT_WAVE_VECTOR)
        flux = wavelength / (self.NORMED_DETECTOR_EFFICIENCY * self.NORMED_DETECTOR_WAVELENGTH)
        return flux

    def parse_header(self, data_line):
        # parses a pair of floats: "($float, $float)"
        def slit_matcher(s):
            regex = r'\(\s*([-+]?[0-9]*\.?[0-9]*)\s*,\s*([-+]?[0-9]*\.?[0-9]*\s*)\)\s*' \
                    r'(\s*[-+]?[0-9]*\.?[0-9]*)\s*x\s*([-+]?[0-9]*\.?[0-9]*)\s*'
            try:
                match = re.search(regex, s, re.IGNORECASE)
                return list(map(float, [match.group(1), match.group(2), match.group(3), match.group(4)]))
            except TypeError:
                print("Wrong data type of pair matching")
            except:
                print("Function does not work: slit_matcher")

        # parses a float, returns the first float from the left, anything else is ignored
        def float_matcher(s):
            regex = r'\s*([-+]?[0-9]*\.?[0-9]*)\s*'
            try:
                match = re.search(regex, s, re.IGNORECASE)
                return float(match.group(1))
            except TypeError:
                print("Wrong data type of float matching")
            except OSError:
                pass

        # only the lines with colon in between are needed
        if ':' in data_line:
            # seems like we have a meta-data in the form "key : value"
            # we have to be careful, since raw_data can have multiple entries (more than 2), because
            # somebody thought it is a good idea to have a separator in the data section..
            # for example: "ms2_status : ok: left_idle ...
            key, value = list(map(lambda s: s.strip(), data_line[1:].split(':', maxsplit=1)))
            if self.PARAMETER_SUFFIX in key:
                key = key.replace(self.PARAMETER_SUFFIX, '')
                if key in self.NORMAL_FLOAT_VALUES:
                    value = float_matcher(value)
                    return key, value
                if key in self.SLITs:
                    values = slit_matcher(value)
                    slit_infos = list(map(lambda s: key + '_' + s, self.SLIT_INFOs))
                    return slit_infos, values
            elif self.COUNT_TIME in key or self.SCAN_INFO in key:
                return key, value

    def scan_variable_getter(self, line1, line2):
        self.scan_name, unit = list(map(lambda s: s.split()[1], (line1, line2)))
        if self.scan_name in self.CHANNELS:
            plot_name = self.COILs_POSITIONs[self.scan_name]
            unit = self._unit_correction(unit)
            self.x_label = r'Current $I_{\mathrm{%s}}$ (%s)' % (plot_name, unit)
        else:
            self.x_label = self.scan_name + ' (' + unit + ')'
        return self.scan_name

    # grabs __properties from a file that are encoded in the form of "key: value"
    # Reforming realised by means of parsing funcs
    def _property_getter(self):
        f = open(self.filename, "r")
        lines = f.readlines()
        f.close()

        scan_key = None
        for line in lines:
            if line.startswith('#'):
                if line.startswith('###'):
                    if self.DATA_START in line:
                        data_start_index = lines.index(line)
                        scan_key = self.scan_variable_getter(lines[data_start_index + 1], lines[data_start_index + 2])
                try:
                    key, value = self.parse_header(line)
                except:
                    continue
                if isinstance(key, str):
                    self.__properties.update({key: value})
                elif isinstance(key, list):
                    self.__properties.update(dict(zip(key, value)))
                # except ValueError:
                #     print(self._parse_header(line))
        self.__properties[scan_key] = 'Scanned'

    def property_sort(self):
        self._property_getter()
        for key in self.INFOS_NEEDED:
            try:
                self.properties_sorted.update({key: self.__properties[key]})
            except KeyError:  # for optional infos
                self.properties_sorted.update({key: 'Unused'})
