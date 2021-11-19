import numpy as np
import re


class DataFile:
    # PATTERN_HEADER = re.compile(r"\**\sDATA\s\**")
    PATTERN_POINTS = re.compile(r"([0-9]*)\sPoints")
    PATTERN_SCAN = re.compile(r"Scanning Variables:\s(tti[1-2]_out[1-2]),\sSteps:\s([1-9]*.*[1-9]*)")
    PATTERN_CURRENTS = re.compile(r"(tti[1-2]_out[1-2])\s*=\s*([1-9]*.*[1-9]*)")

    DATA_HEADER = "**************************** DATA ******************************************\n"
    SCANNING_VARIABLE_KEY = "Scanning Variables"
    NORMED_DETECTOR_EFFICIENCY = 3.24e4
    NORMED_DETECTOR_WAVELENGTH = 1.8e-10  # m
    INCIDENT_WAVE_VECTOR = 1.4e10  # m-1

    FOLDER = "../20210907/"
    FILE_PREFIX = 'morpheus2021n'
    FILE_FORMAT = '.dat'

    POWER_SUPPLY_PREFIX = "tti"

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

    CURRENT_O2 = np.NAN
    CURRENT_I1 = np.NAN
    CURRENT_I2 = np.NAN
    CURRENT_O1 = np.NAN
    CURRENTS = [CURRENT_I1, CURRENT_O1, CURRENT_I2, CURRENT_O2]
    COILS_CURRENTS = dict(zip(CHANNELS, CURRENTS))

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
        self.filename = self.FOLDER + self.FILE_PREFIX + "{:0>6d}".format(file_index) + self.FILE_FORMAT
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
            if re.match(self.PATTERN_CURRENTS, line):
                # print(line)
                channel, current = re.search(self.PATTERN_CURRENTS, line).groups()
                if channel in self.CHANNELS:
                    # print(channel, current)
                    self.COILS_CURRENTS[channel] = float(current)

            if self.SCANNING_VARIABLE_KEY in line:
                self.scan_name, self.scan_step = re.search(self.PATTERN_SCAN, line).groups()
                self.scan_step = float(self.scan_step)
        if self.scan_name == self.COIL_O1:
            self.pairing_coil = self.COIL_O2
        elif self.scan_name == self.COIL_O2:
            self.pairing_coil = self.COIL_O1
        elif self.scan_name == self.COIL_I1:
            self.pairing_coil = self.COIL_I2
        elif self.scan_name == self.COIL_I2:
            self.pairing_coil = self.COIL_I1
        else:
            raise RuntimeError("Invalid scan name captured: {}".format(self.scan_name))

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
