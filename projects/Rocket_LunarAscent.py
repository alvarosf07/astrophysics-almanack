import numpy as np
import matplotlib.pyplot as plt
import openmdao.api as om
from openmdao.utils.assert_utils import assert_near_equal

import dymos as dm
from dymos.examples.polarimplementation.polar_eom import polar_eom


p = om.Problem(model=om.Group())

p.driver = om.pyOptSparseDriver()
p.driver.options['optimizer'] = 'SLSQP'
p.driver.declare_coloring()
#p.driver = om.ScipyOptimizeDriver()

traj = dm.Trajectory()

traj.add_parameter('c', opt=False, val=1.5, units='DU/TU',
                   targets={'phase1': ['c']})

# First Phase
phase1 = dm.Phase(ode_class=FiniteBurnODE,
                 transcription=dm.GaussLobatto(num_segments=5, order=3, compressed=False))

phase1 = traj.add_phase('phase1', burn1)

phase1.set_time_options(fix_initial=True, duration_bounds=(.5, 10), units='TU')
phase1.add_state('r', fix_initial=True, fix_final=False, defect_scaler=100.0,
                rate_source='r_dot', units='DU')
phase1.add_state('theta', fix_initial=True, fix_final=False, defect_scaler=100.0,
                rate_source='theta_dot', units='rad')
phase1.add_state('vr', fix_initial=True, fix_final=False, defect_scaler=100.0,
                rate_source='vr_dot', units='DU/TU')
phase1.add_state('vt', fix_initial=True, fix_final=False, defect_scaler=100.0,
                rate_source='vt_dot', units='DU/TU')
phase1.add_state('accel', fix_initial=True, fix_final=False,
                rate_source='at_dot', units='DU/TU**2')
phase1.add_state('deltav', fix_initial=True, fix_final=False,
                rate_source='deltav_dot', units='DU/TU')
phase1.add_control('u1', rate_continuity=True, rate2_continuity=True, units='deg',
                  scaler=0.01, rate_continuity_scaler=0.001, rate2_continuity_scaler=0.001,
                  lower=-30, upper=30)

phase1.add_timeseries_output('pos_x')
burn1.add_timeseries_output('pos_y')

# Link Phases
#traj.link_phases(phases=['burn1', 'coast', 'burn2'],
                 #vars=['time', 'r', 'theta', 'vr', 'vt', 'deltav'])

#traj.link_phases(phases=['burn1', 'burn2'], vars=['accel'])

p.model.add_subsystem('traj', subsys=traj)

# Finish Problem Setup

# Needed to move the direct solver down into the phases for use with MPI.
#  - After moving down, used fewer iterations (about 30 less)

p.driver.add_recorder(om.SqliteRecorder('two_burn_orbit_raise_example.db'))

p.setup(check=True, mode='fwd')

# Set Initial Guesses
p.set_val('traj.parameters:c', value=1.5, units='DU/TU')

phase1 = p.model.traj.phases.phase1

p.set_val('traj.phase1.t_initial', value=0.0)
p.set_val('traj.phase1.t_duration', value=2.25)
p.set_val('traj.phase1.states:r', value=phase1.interpolate(ys=[1, 1.5],
                                                         nodes='state_input'))
p.set_val('traj.phase1.states:theta', value=phase1.interpolate(ys=[0, 1.7],
                                                             nodes='state_input'))
p.set_val('traj.phase1.states:vr', value=phase1.interpolate(ys=[0, 0],
                                                          nodes='state_input'))
p.set_val('traj.phase1.states:vt', value=phase1.interpolate(ys=[1, 1],
                                                          nodes='state_input'))
p.set_val('traj.phase1.states:accel', value=phase1.interpolate(ys=[0.1, 0],
                                                             nodes='state_input'))
p.set_val('traj.phase1.states:deltav', value=phase1.interpolate(ys=[0, 0.1],
                                                              nodes='state_input'))
p.set_val('traj.phase1.controls:u1',
          value=phase1.interpolate(ys=[-3.5, 13.0], nodes='control_input'))

dm.run_problem(p)


#
# Plot results
#
traj = p.model.traj
exp_out = traj.simulate()

fig = plt.figure(figsize=(8, 4))
fig.suptitle('Two Burn Orbit Raise Solution')
ax_u1 = plt.subplot2grid((2, 2), (0, 0))
ax_deltav = plt.subplot2grid((2, 2), (1, 0))
ax_xy = plt.subplot2grid((2, 2), (0, 1), rowspan=2)

span = np.linspace(0, 2 * np.pi, 100)
ax_xy.plot(np.cos(span), np.sin(span), 'k--', lw=1)
ax_xy.plot(3 * np.cos(span), 3 * np.sin(span), 'k--', lw=1)
ax_xy.set_xlim(-4.5, 4.5)
ax_xy.set_ylim(-4.5, 4.5)

ax_xy.set_xlabel('x ($R_e$)')
ax_xy.set_ylabel('y ($R_e$)')

ax_u1.set_xlabel('time ($TU$)')
ax_u1.set_ylabel('$u_1$ ($deg$)')
ax_u1.grid(True)

ax_deltav.set_xlabel('time ($TU$)')
ax_deltav.set_ylabel('${\Delta}v$ ($DU/TU$)')
ax_deltav.grid(True)

t_sol = dict((phs, p.get_val('traj.{0}.timeseries.time'.format(phs)))
             for phs in ['phase1'])
x_sol = dict((phs, p.get_val('traj.{0}.timeseries.pos_x'.format(phs)))
             for phs in ['phase1'])
y_sol = dict((phs, p.get_val('traj.{0}.timeseries.pos_y'.format(phs)))
             for phs in ['phase1'])
dv_sol = dict((phs, p.get_val('traj.{0}.timeseries.states:deltav'.format(phs)))
              for phs in ['phase1'])
u1_sol = dict((phs, p.get_val('traj.{0}.timeseries.controls:u1'.format(phs), units='deg'))
              for phs in ['phase1'])

t_exp = dict((phs, exp_out.get_val('traj.{0}.timeseries.time'.format(phs)))
             for phs in ['phase1'])
x_exp = dict((phs, exp_out.get_val('traj.{0}.timeseries.pos_x'.format(phs)))
             for phs in ['phase1'])
y_exp = dict((phs, exp_out.get_val('traj.{0}.timeseries.pos_y'.format(phs)))
             for phs in ['phase1'])
dv_exp = dict((phs, exp_out.get_val('traj.{0}.timeseries.states:deltav'.format(phs)))
              for phs in ['phase1'])
u1_exp = dict((phs, exp_out.get_val('traj.{0}.timeseries.controls:u1'.format(phs),
                                    units='deg'))
              for phs in ['phase1'])

for phs in ['phase1']:
    try:
        ax_u1.plot(t_exp[phs], u1_exp[phs], '-', marker=None, color='C0')
        ax_u1.plot(t_sol[phs], u1_sol[phs], 'o', mfc='C1', mec='C1', ms=3)
    except KeyError:
        pass

    ax_deltav.plot(t_exp[phs], dv_exp[phs], '-', marker=None, color='C0')
    ax_deltav.plot(t_sol[phs], dv_sol[phs], 'o', mfc='C1', mec='C1', ms=3)

    ax_xy.plot(x_exp[phs], y_exp[phs], '-', marker=None, color='C0', label='explicit')
    ax_xy.plot(x_sol[phs], y_sol[phs], 'o', mfc='C1', mec='C1', ms=3, label='implicit')

plt.show()
