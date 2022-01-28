import numpy as np
from data_file import DataFile
import matplotlib.pyplot as plt
from lmfit import Model
import universal_params as univ
from matrix_angles import polarisation2angles

# from lmfit.lineshapes import gaussian2d
plt.rcParams.update({'font.size': 18})
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams['font.sans-serif'] = ['Arial']
SCAN_2D = range(2814, 2854 + 1)


def sine_func_1d(current, amp, para_coil, phase, y_shift):
    # minus sein before c because the two coils are in the opposite direction
    return -amp * np.cos(-para_coil * current + phase) + y_shift


def sine_func_2d(current1, current2, amp, para_coil1, para_coil2, delta_i, delta_f, shift1, shift2, phase_shift,
                 y_shift):
    # minus sein before c because the two coils are in the opposite direction
    # return amp * np.cos(para_coil1 * (current1 + shift1) - para_coil2 * (current2 + shift2) - (
    #         delta_i - delta_f)) + y_shift  # + phase_shift
    return -amp * np.cos(para_coil1 * (current1 + shift1) - para_coil2 * (current2 + shift2) + phase_shift + (
            delta_i - delta_f)) + y_shift


def phase2current(coil_name, precession_phase, shift, para_coil1=None, para_coil2=None, shift1=None, shift2=None):
    if precession_phase == 0:
        return 0
    if coil_name == univ.COILS_POSITIONS[univ.COIL_I1]:
        if para_coil1 is None or shift1 is None:
            raise ValueError("The parameters for Coil I1 are invalid.")
        else:
            current1 = (precession_phase - shift / 2.0) / para_coil1 - shift1
            return current1
    elif coil_name == univ.COILS_POSITIONS[univ.COIL_I2]:
        if para_coil2 is None or shift2 is None:
            raise ValueError("The parameters for Coil I2 are invalid.")
        else:
            current2 = -(precession_phase - shift / 2.0) / para_coil2 - shift2
            return current2
    else:
        raise ValueError("Invalid coil name {:s}".format(coil_name))


def counts2polarisation(count, params):
    # intensity_total = 2 * sine_params[-1]
    pol_max = params[0] / params[-1]
    polarisation = (count - params[-1]) / params[0] * pol_max
    return polarisation


scan_steps = len(SCAN_2D)
coil_scanned = np.empty((scan_steps, scan_steps))
coil_pair = np.empty_like(coil_scanned)
coilname_scanned = None
coilname_pair = None
counts = np.empty_like(coil_scanned)
pol_1dfit = np.empty_like(coil_scanned)

current_pairing_coil = 0.0

fmodel_1d = Model(sine_func_1d, independent_vars=["current"])
# fmodel = Model(sine_function_2d, independent_vars=["current1", "current2"])
fmodel = Model(sine_func_2d, independent_vars=["current1", "current2"])


# create parameters -- these are named from the function arguments --
# giving initial values
def lm_fit_1d(xdata, ydata):
    fmodel_1d.set_param_hint('amp', value=(np.max(ydata) - np.min(ydata)) / 2.0,
                             min=(np.max(ydata) - np.min(ydata)) / 2.5, max=(np.max(ydata) - np.min(ydata)) / 1.5)
    fmodel_1d.set_param_hint('para_coil', value=np.pi / 2.5, min=np.pi / 3.5, max=np.pi / 2)
    fmodel_1d.set_param_hint('phase', value=0, min=- np.pi, max=np.pi)
    fmodel_1d.set_param_hint('y_shift', value=(np.max(ydata) + np.min(ydata)) / 2.0, min=np.min(ydata),
                             max=np.max(ydata))
    fmodel.set_param_hint('delta_i', value=delta_i, vary=False)
    fmodel.set_param_hint('delta_f', value=delta_f, vary=False)
    params = fmodel_1d.make_params()
    result = fmodel_1d.fit(ydata, params, current=xdata)
    return result.params["amp"].value, result.params["para_coil"].value, result.params["phase"].value, result.params[
        "y_shift"].value


def lm_fit_2d(xdata1, xdata2, ydata, delta_i, delta_f):
    # delta_i and delta_f are the scattering angles
    fmodel.set_param_hint('amp', value=(np.max(ydata) - np.min(ydata)) / 2.0, min=(np.max(ydata) - np.min(ydata)) / 3.0,
                          max=np.max(ydata) - np.min(ydata))
    fmodel.set_param_hint('para_coil1', value=np.pi / 2.5, min=np.pi / 3.5, max=np.pi / 2)
    fmodel.set_param_hint('para_coil2', value=np.pi / 2.5, min=np.pi / 3.5, max=np.pi / 2)
    fmodel.set_param_hint('delta_i', value=delta_i, vary=False)
    fmodel.set_param_hint('delta_f', value=delta_f, vary=False)
    fmodel.set_param_hint('shift1', value=0, vary=False, min=-0.05, max=0.05)
    fmodel.set_param_hint('shift2', value=0, vary=False, min=-0.05, max=0.05)
    fmodel.set_param_hint('phase_shift', value=0, min=-0.2 * np.pi, max=0.2 * np.pi)
    fmodel.set_param_hint('y_shift', value=(np.max(ydata) + np.min(ydata)) / 2.0, min=np.min(ydata), max=np.max(ydata))
    params = fmodel.make_params()
    result = fmodel.fit(ydata, params, current1=xdata1, current2=xdata2)
    print(result.fit_report())
    return result.params["amp"].value, result.params["para_coil1"].value, result.params["para_coil2"].value, \
           result.params["delta_i"].value, result.params["delta_f"].value, result.params["shift1"].value, result.params[
               "shift2"].value, result.params["phase_shift"].value, result.params["y_shift"].value


def current_adjust(coil, period):
    if coil < 4 - period:
        coil += period
    elif coil > period:
        coil -= period
    else:
        pass
    return coil


d_spacing = 3.55e-10  # PG 002
# d_spacing = 4.556e-10  #MnSi
wavelength = 4.905921e-10
theta = np.arcsin(wavelength / (2 * d_spacing))
delta_i = theta
delta_f = np.pi - theta

for counter, scan_2d in enumerate(SCAN_2D):
    filename = "Scan{:d}.png".format(scan_2d)
    data_file = DataFile(scan_2d)
    if data_file.complete_scan is False:
        print("Scan No. {} is not complete".format(scan_2d))
        continue
    coilname_scanned = data_file.COILs_POSITIONs[data_file.scan_name]
    coilname_pair = data_file.COILs_POSITIONs[data_file.pair_coil]
    coil_scanned[counter, :] = data_file.scan_x
    coil_pair[counter, :] = data_file.COILS_CURRENTS[data_file.pair_coil]
    counts[counter, :] = data_file.scan_count
    sine_params_1d = lm_fit_1d(data_file.scan_x, data_file.scan_count)
    counts_1dfit = sine_func_1d(data_file.scan_x, *sine_params_1d)
    pol_1dfit[counter, :] = counts2polarisation(counts_1dfit, sine_params_1d)
    if counter == 0:
        phase = sine_params_1d[2] - (delta_i - delta_f)
        print(sine_params_1d, phase)
    # print(sine_params_1d, sine_params_1d[0] / sine_params_1d[-1])

sine_params2d = lm_fit_2d(coil_pair, coil_scanned, counts, delta_i, delta_f)
print(sine_params2d, sine_params2d[0] / sine_params2d[-1], np.sin(sine_params2d[-2]))
counts_2dfit = sine_func_2d(coil_pair, coil_scanned, *sine_params2d)
pol_2dfit = counts2polarisation(counts_2dfit, sine_params2d)

filename1 = "PG_InnerCoils_1DFitted.png"
fig1, ax1 = plt.subplots()
cnt1 = ax1.contourf(coil_scanned, coil_pair, pol_1dfit)
ax1.set_xlabel(r"$I$ (A): {:s}".format(coilname_scanned), fontname='sans-serif')
ax1.set_ylabel(r"$I$ (A): {:s}".format(coilname_pair), fontname="Arial")
cbar1 = fig1.colorbar(cnt1)
ax1.tick_params(axis="both", direction="in")
fig1.savefig(filename1, bbox_inches='tight')
plt.close(fig1)

filename2 = "PG_InnerCoils_2DFitted.png"
fig2, ax2 = plt.subplots()
cnt2 = ax2.contourf(coil_scanned, coil_pair, pol_2dfit)
ax2.set_xlabel(r"$I$ (A): {:s}".format(coilname_scanned), fontname="Arial")
ax2.set_ylabel(r"$I$ (A): {:s}".format(coilname_pair), fontname="Arial")
cbar2 = fig2.colorbar(cnt2)
ax2.tick_params(axis="both", direction="in")
fig2.savefig(filename2, bbox_inches='tight')
plt.close(fig2)

filename3 = "PG_InnerCoils_Counts.png"
fig3, ax3 = plt.subplots()
cnt3 = ax3.contourf(coil_scanned, coil_pair, counts)
ax3.set_xlabel(r"$I$ (A): {:s}".format(coilname_scanned), fontname="Arial")
ax3.set_ylabel(r"$I$ (A): {:s}".format(coilname_pair), fontname="Arial")
cbar3 = fig3.colorbar(cnt3)
ax3.tick_params(axis="both", direction="in")
fig3.savefig(filename3, bbox_inches='tight')
plt.close(fig3)

period1 = 2 * np.pi / sine_params2d[1]
period2 = 2 * np.pi / sine_params2d[2]
phase_shift = sine_params2d[-2]
print(period1, period2, np.rad2deg(phase_shift))
print("theta={:.2f}°, coil scanned: {:s}, coil pair: {:s}".format(np.rad2deg(theta), coilname_scanned, coilname_pair))

# precession_coil1 = np.pi
# current_coil1 = phase2current(coil_name=univ.COILS_POSITIONS[univ.COIL_I1], precession_phase=precession_coil1,
#                               para_coil1=sine_params2d[1], shift1=sine_params2d[5])
AXIS_X = "x"
AXIS_Y = "y"
AXIS_Z = "z"
AXES = [AXIS_X, AXIS_Y, AXIS_Z]
for pi in AXES:
    for pf in AXES:
        alpha_i, beta_i, alpha_f, beta_f = polarisation2angles(pi, pf, 2 * theta)
        coil_i1 = phase2current(coil_name=univ.COILS_POSITIONS[univ.COIL_I1], precession_phase=beta_i,
                                shift=phase_shift, para_coil1=sine_params2d[1], shift1=sine_params2d[5])
        coil_i2 = phase2current(coil_name=univ.COILS_POSITIONS[univ.COIL_I2], precession_phase=beta_f,
                                shift=phase_shift, para_coil2=sine_params2d[2], shift2=sine_params2d[6])
        coil_i1 = current_adjust(coil_i1, period1)
        coil_i2 = current_adjust(coil_i2, period2)

        if pi != AXIS_Z and pf != AXIS_Z:
            count = sine_func_2d(coil_i1, coil_i2, *sine_params2d)
            pol = counts2polarisation(count, sine_params2d)
            print("polarisation = {:.2f}".format(pol))

        print(
            "pi = {:s}, pf = {:s}, coil_o1 = {:.1f}°, coil_i1 = {:.1f}° at {:.2f} A, coil_i2 = {:.1f}° at {:.2f} A, coil_o2 = {:.1f}°".format(
                pi, pf, np.rad2deg(alpha_i), np.rad2deg(beta_i), coil_i1, np.rad2deg(beta_f), coil_i2,
                np.rad2deg(alpha_f)))
