from stytra import Stytra
import pandas as pd
import numpy as np
from stytra.stimulation.stimuli import (MovingGratingStimulus)
from stytra.stimulation import (Protocol, DynamicStimulus)
from lightparam import Param
import random

class customWrapper(DynamicStimulus):
    def __init__(
        self,
        stim,
        left = 60,
        right = 600,
        **kwargs
    ):
        super().__init__(**kwargs)
        self.name = "conditional"
        self.stim_on = stim
        self._elapsed_difference = 0
        self._elapsed_when_centering_started = 0
        self.reset_phase_shift = 1
        self.left_border  = left
        self.right_border = right
        self.advance = False
        self.dynamic_parameters.append("velx")
        self.velx = 0
        self.duration = self.stim_on.duration
        self.stimulus_dynamic = hasattr(stim, "dynamic_parameters")
        self._dt = 0
        self._past_t = 0
        self._previous_value = False
        self.last_bound = 0
        self.velx_start = self.stim_on.df_param.vel_x

    @property
    def dynamic_parameter_names(self):
        if self.stimulus_dynamic:
            return (
                super().dynamic_parameter_names + self.stim_on.dynamic_parameter_names
            )
        else:
            return super().dynamic_parameter_names

    def get_dynamic_state(self):
        state = super().get_dynamic_state()
        if self.stimulus_dynamic:
            state.update(self.stim_on.get_dynamic_state())
        return state


    def initialise_external(self, experiment):
        super().initialise_external(experiment)
        self.stim_on.initialise_external(experiment)

    def get_state(self):
        state = super().get_state()
        state.update(
            {"On": self.stim_on.get_state()}
        )
        return state

    def start(self):
        super().start()
        self.stim_on.start()

    def check_condition(self):
        x, y, theta = self._experiment.estimator.get_position()
        if y<self.left_border and self.last_bound != 1:
            self.last_bound = 1
            return True
        elif y>self.right_border and self.last_bound != 2:
            self.last_bound = 2
            return True
        else:
            return False

    def update(self):
        self._dt = self._elapsed - self._past_t
        self._past_t = self._elapsed
        if not hasattr(self.stim_on, "vel_x"):
            self.stim_on.vel_x = 0
        self.velx = self.stim_on.vel_x
            
        if self.check_condition():
            # put speed reversal here
            if self.last_bound == 2:
                self.stim_on.df_param.vel_x = self.velx_start*-1
            elif self.last_bound == 1:
                self.stim_on.df_param.vel_x = self.velx_start
            
            if self.velx != 0:
                self.advance = True

                new_phase = min(self.stim_on.current_phase + self.reset_phase_shift, len(self.stim_on.phase_times)-1)
                
                time_added = (
                    self._elapsed
                    - self._elapsed_difference
                    - self.stim_on.phase_times[new_phase]
                )
                self.duration += time_added
                self._elapsed_difference += time_added

        # update the current stimulus
        # if not self.advance:
        self.stim_on._elapsed = self._elapsed - self._elapsed_difference

        self.stim_on.update()

    def paint(self, p, w, h):
        # p.setBrush(QBrush(QColor(0, 0, 0)))
        # p.drawRect(QRect(-1, -1, w + 2, h + 2))
        self.stim_on.paint(p, w, h)
       

class OMR(Protocol):
    name = "free_swimming_omr"

    # setup camera
    stytra_config = dict(
        camera=dict(type="spinnaker", roi=[272, 850, 1504, 350]),
        tracking=dict(method="fish", embedded=False, estimator="position"),
        recording=dict(extension="mp4", kbit_rate=2000))
    
    def __init__(self):
        super().__init__()

        self.t_pre = Param(5.0)  # time of still gratings before they move
        self.t_move = Param(5.0)  # time of gratings movement
        self.grating_vel = Param(10.0)  # gratings velocity
        self.grating_period = Param(50)  # grating spatial period
        self.grating_angle_deg = Param(90.0)  # grating orientation
        self.Red = Param(255, (0,255))
        self.Green = Param(255, (0,255))
        self.Blue  = Param(255, (0,255))
        self.left_border = Param(10, (0, 2000))
        self.right_border = Param(500, (0, 2000))
        self.repeats = Param(10, (1, 100))
        self.pre_protocol_t = Param(60, (0, 1000))
 
    def get_stim_sequence(self):
        # Use six points to specify the velocity step to be interpolated:
        t_base = [0, self.t_pre, self.t_pre, self.t_pre + self.t_move]
        vel_base = [0, 0, self.grating_vel, self.grating_vel]
        vel = []
        t = []
        
        speeds = [1, 2, 3, 4, 5, 6]
        
        
        
        for n in range(self.repeats):
            rand_speeds = random.sample(speeds,len(speeds))
            for k in rand_speeds:
                vel.extend([x*k for x in vel_base])
                if len(t)>0:
                    t_add = t[-1]
                else:
                    t_add = 0
                        
                t.extend([x+t_add for x in t_base])

        df1 = pd.DataFrame(dict(t=t, vel_x=vel))
        
        t_0 = [0,self.pre_protocol_t]
        vel_0  = [0, 0]
        df0 = pd.DataFrame(dict(t=t_0, vel_x=vel_0))
        gratStim = []
        gratStim.append(MovingGratingStimulus(df_param=df0, grating_period = self.grating_period,
                        grating_angle=self.grating_angle_deg* np.pi / 180,
                        grating_col_1=(self.Red, self.Green, self.Blue)))
        
        gratStim.append(customWrapper(stim = MovingGratingStimulus(df_param=df1,
                        grating_period = self.grating_period,
                        grating_angle=self.grating_angle_deg* np.pi / 180,
                        grating_col_1=(self.Red, self.Green, self.Blue)), left = self.left_border, right = self.right_border))

        return gratStim

if __name__ == "__main__":
    s = Stytra(protocol=OMR())
